import pandas
from selenium import webdriver
from selenium.common import exceptions
from pymongo import MongoClient
import time


def get_measures(name: str, table: str, sheet=0):
    # 函数说明：
    # 功能：给定家畜/家禽名称，从给定的excel表格中获取对应的繁殖力指标
    # 参数：name：要获取繁殖力指标的家畜/家禽的名称；table：包含繁殖力指标的excel表格的路径；sheet：要读取的工作簿
    # 返回值：繁殖力指标的列表
    measures_df = pandas.read_excel(table, sheet_name=sheet)
    name_index = ''
    measure_list = []
    for measure_name in measures_df.iloc[0:0, ].columns:
        if measure_name == name or name in measure_name:
            name_index = measure_name
    if name_index == '':
        print('找不到' + name + '对应的指标')
        return
    else:
        for measure in measures_df.loc[:, [name_index]][name_index]:
            if type(measure) is not str:
                continue
            measure_list.append(measure)
    return measure_list


def get_names_dic(table: str, sheet=0):
    # 函数说明：
    # 功能：从给定的excel表格中获取家畜的种类及对应的亚种的名称
    # 参数：table：表格的路径，sheet：要读取的工作簿
    # 返回值：种名：亚种名集合的字典列表
    name_df = pandas.read_excel(table, sheet_name=sheet)
    # 读取表格第一行，即所有种的名称，作为后面读取每一列的索引
    name_dic_list = []
    names_list = name_df.iloc[0:0, ].columns
    # 分别读取每一列的所有亚种名，形成一个集合，并加入字典中
    for name_index in names_list:
        name_list = name_df.loc[:, [name_index]][name_index]
        name_dic = {}
        nset = set()
        for name in name_list:
            if type(name) is not str:
                continue
            nset.add(name)
        nset.add(name_index)
        name_dic[name_index] = nset
        name_dic_list.append(name_dic)
    return name_dic_list


def get_search_phase_list(ndic: dict, flag=1):
    # 函数说明：
    # 功能：通过种名从指标的表格中获取对应的繁殖力指标，并构建用于搜索的语句（针对CNKI）。
    # 参数：ndic：“种名:亚种名集合”集合格式的字典；flag：1.中国知网，2.万方，3.重庆维普
    # 返回值：构建的搜索语句的列表
    if flag not in [1, 2, 3]:
        print("参数错误，flag参数应为1、2、3中的一个")
        return None
    name = list(ndic.keys())[0]  # 种名
    mlist = get_measures(name, '名称与指标.xlsx', 1)
    if isinstance(mlist, type(None)):
        return None
    nlist = ndic[name]
    search_phase_list = []
    i = 0
    ignore_name_list = []
    for n in nlist:
        # 跳过已经搜过的品种
        if n in ignore_name_list:
            continue
        for m in mlist:
            # 跳过当前品种中已经搜过的部分
            if (n == '鸡') & (i < 18):
                i += 1
                continue
            if flag == 1:
                search_phase = 'AB = ‘' + n + '’*‘' + m + '’'
            elif flag == 2:
                search_phase = "摘要:(\"" + n + "\"*\"" + m + "\")"
            else:
                search_phase = 'R=' + n + ' * ' + m
            search_phase_list.append(search_phase)
    return search_phase_list


def get_ref_list(ref_container, search_phase: str):
    # 函数说明
    # 从一个包含所有文献信息的网页元素中获取所有文献（针对中国知网）
    # 参数：ref_container：包含文献信息的网页元素
    # 返回值：一个文献信息的字典列表
    try:
        refs = ref_container.find_elements_by_tag_name('tr')
    except AttributeError:
        print('传入的参数类型错误！')
        return None
    try:
        refs = ref_container.find_elements_by_tag_name('tr')
    except exceptions.NoSuchElementException:
        print('获取文献列表失败！')
        return None
    else:
        ref_list = []
        for ref in refs:
            try:
                # 获取文章的基本信息，包括标题、作者、期刊和发布日期
                title = ref.find_element_by_class_name('name').text
                authors = ref.find_element_by_class_name('author').text
                source = ref.find_element_by_class_name('source').text
                date = ref.find_element_by_class_name('date').text
                # 如果文章支持HTML在线阅读，则获取在线阅读的链接，否则获取摘要恋姬
                operat = ref.find_element_by_class_name('operat')
                try:
                    html_link = operat.find_element_by_class_name('icon-html')
                    link = html_link.get_attribute('href')
                except exceptions.NoSuchElementException:
                    try:
                        read_link = operat.find_element_by_class_name(
                            'icon-read')
                    except exceptions.NoSuchElementException:
                        link = 'Null'
                    else:
                        link = read_link.get_attribute('href')
            except exceptions.StaleElementReferenceException:
                return 'stateelement'
            ref_dic = {
                '标题': title,
                '作者': authors,
                '期刊名': source,
                '出版日期': date,
                '摘要或全文链接': link,
                '搜索语句': search_phase
            }
            ref_list.append(ref_dic)
        return ref_list


def get_all_refs_with_search_pahse(search_phase: str):
    # 函数说明：
    # 功能：根据搜索语句爬取文献(针对中国知网)
    # 参数：search_pahse：用于搜索文献的搜索语句
    # 返回值：包含文献信息的字典所组成的列表
    link = 'https://kns.cnki.net/kns8/AdvSearch?dbprefix=SCDB&&crossDbcodes=CJFQ%2CCDMD%2CCIPD%2CCCND%2CCISD%2CSNAD%2CBDZK%2CGXDB_SECTION%2CCJFN%2CCCJD'
    driver.get(link)
    # 等待网页内容加载完成
    time.sleep(1)
    # 点击“专业检索”按钮进入专业检索模式
    try:
        driver.find_element_by_xpath('/html/body/div[4]/div/div[2]/ul/li[4]')
    except exceptions.NoSuchElementException:
        print('没有找到“专业检索”按钮')
    else:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/ul/li[4]').click()
        time.sleep(0.5)
    # 获取检索语句输入区域
    try:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/div/div[1]/div[1]/div[2]/textarea')
    except exceptions.NoSuchElementException:
        print('没有找到搜索框')
    else:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/div/div[1]/div[1]/div[2]/textarea'
        ).send_keys(search_phase)
    # 取消中英文扩展
    try:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/div/div[1]/div[1]/div[2]/div[1]/div[1]/span[3]/label[1]/input'
        )
    except exceptions.NoSuchElementException:
        print('没有找到中英文扩展选项')
    else:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/div/div[1]/div[1]/div[2]/div[1]/div[1]/span[3]/label[1]/input'
        ).click()
        time.sleep(0.5)
    # 点击搜索按钮
    try:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/div/div[1]/div[1]/div[2]/div[2]/input'
        )
    except exceptions.NoSuchElementException:
        print('没有找到搜索按钮')
    else:
        driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[2]/div/div[1]/div[1]/div[2]/div[2]/input'
        ).click()
    time.sleep(1)
    # 获取包含本页所有文献的对象
    ref_list = []
    # 循环获取所有页的文献信息
    while True:
        try:
            ref_container = driver.find_element_by_xpath(
                '//*[@id="gridTable"]/table/tbody')
        except exceptions.NoSuchElementException:
            print('没有找到文献table！')
            break
        else:
            # 获取本页所有文献，保存为字典
            rlist = get_ref_list(ref_container, search_phase)
            if rlist == 'stateelement':
                time.sleep(2)
                try:
                    ref_container = driver.find_element_by_xpath(
                        '//*[@id="gridTable"]/table/tbody')
                except exceptions.NoSuchElementException:
                    print('没有找到文献table！')
                else:
                    rlist = get_ref_list(ref_container, search_phase)
            if (rlist is None) | (rlist == 'stateelement'):
                print('从传递的参数中获取文献列表失败！')
            else:
                for ref in rlist:
                    ref_list.append(ref)
        # 看是否有下一页，如果有就进入下一页并获取文献信息，否则退出循环
        try:
            next_page_button = driver.find_element_by_xpath(
                '//*[@id="PageNext"]')
        except exceptions.NoSuchElementException:
            break
        else:
            next_page_button.click()
            time.sleep(2)
    return ref_list


def get_all_refs_WF(search_phase: str):
    # 函数说明：
    # 功能：用指定的检索语句在万方检索，获取所有的检索结果
    # 参数：search_phase：检索语句（摘要:("名称"*"繁殖力指标")）
    # 返回值：文献信息字典组成的列表
    link = 'http://new.wanfangdata.com.cn/searchResult/getAdvancedSearch.do?searchType=all'
    driver.get(link)
    # 等待网页内容加载完成
    time.sleep(1)
    # 点击专业搜索按钮
    try:
        pro_button = driver.find_element_by_xpath('//*[@id="expert_search_a"]')
    except exceptions.NoSuchElementException:
        print('没有找到专业搜索按钮')
        return None
    else:
        pro_button.click()
    time.sleep(0.5)
    # 输入检索条件
    try:
        text_area = driver.find_element_by_xpath(
            '//*[@id="expert_search_textarea"]')
    except exceptions.NoSuchElementException:
        print('没有找到输入检索条件的区域')
        return None
    else:
        text_area.send_keys(search_phase)
    time.sleep(0.3)
    try:
        search_button = driver.find_element_by_xpath('//*[@id="ch_button"]')
    except exceptions.NoSuchElementException:
        print('没有找到搜索按钮')
        return None
    else:
        search_button.click()
    time.sleep(2)
    ref_list = []
    while True:
        for i in range(20):
            xpath = '//*[@id="list_div_aa_1"]/div/'
            for j in range(i):
                xpath += 'div[3]/'
            xpath += 'div[2]'
            try:
                ref_type = driver.find_element_by_xpath(xpath +
                                                        '/div[1]/strong').text
                title = driver.find_element_by_xpath(xpath + '/div[1]/a').text
                href = driver.find_element_by_xpath(
                    xpath + '/div[1]/a').get_attribute('href')
            except exceptions.NoSuchElementException:
                continue
            else:
                try:
                    authors = driver.find_element_by_xpath(
                        xpath + '/div[2]/div[1]').text
                except exceptions.NoSuchElementException:
                    authors = ''
                if ref_type == '[期刊论文]' or ref_type == '[会议论文]':
                    source = driver.find_element_by_xpath(
                        xpath + '/div[2]/div[2]').text
                    try:
                        date = driver.find_element_by_xpath(
                            xpath + '/div[2]/div[4]').text
                    except exceptions.NoSuchElementException:
                        date = ''
                elif ref_type == '[学位论文]':
                    source = driver.find_element_by_xpath(
                        xpath + '/div[2]/div[3]').text
                    try:
                        date = driver.find_element_by_xpath(
                            xpath + '/div[2]/div[5]').text
                    except exceptions.NoSuchElementException:
                        date = ''
                else:
                    print("错误的文献类型：" + ref_type)
                    continue
                # 去除特殊字符
                source = source.replace('《', '')
                source = source.replace('》', '')
                date = date.replace('- ', '')
                ref_type = ref_type.replace('[', '')
                ref_type = ref_type.replace(']', '')
                # 构建文献字典
                ref_dic = {
                    '标题': title,
                    '作者': authors,
                    '期刊名': source,
                    '出版日期': date,
                    '摘要或全文链接': href,
                    '文献类型': ref_type,
                    '搜索语句': search_phase
                }
                ref_list.append(ref_dic)
        try:
            page_area = driver.find_element_by_xpath('//*[@id="laypage_1"]')
            next_page_button = page_area.find_element_by_class_name(
                'laypage_next')
        except exceptions.NoSuchElementException:
            print("没有找到下一页按钮")
            break
        else:
            next_page_button.click()
            time.sleep(3)
    return ref_list


def get_all_refs_WP(search_phase: str):
    # 函数说明：
    # 功能：用指定的检索语句在维普检索，获取所有的检索结果
    # 参数：search_phase：检索语句（R=品种名 * 繁殖力指标）
    # 返回值：文献信息字典组成的列表
    link = 'http://qikan.cqvip.com/Qikan/Search/Advance?from=index'
    driver.get(link)
    # 进入检索式检索模式
    try:
        pro_search = driver.find_element_by_xpath(
            '//*[@id="body"]/div/div[1]/div/ul/li[2]')
    except exceptions.NoSuchElementException:
        print('没有找到检索式检索按钮')
        return
    else:
        pro_search.click()
    # 输入检索式
    try:
        search_pahse_box = driver.find_element_by_xpath(
            '//*[@id="re_searchdomainfilter"]/div[2]/textarea')
    except exceptions.NoSuchElementException:
        print('没有找到输入检索式的位置')
        return
    else:
        search_pahse_box.send_keys(search_phase)
    # 点击检索按钮
    try:
        search_button = driver.find_element_by_xpath(
            '//*[@id="re_searchdomainfilter"]/div[4]/button[1]')
    except exceptions.NoSuchElementException:
        print('没有找到检索按钮')
        return
    else:
        search_button.click()
    time.sleep(2)
    ref_list = []
    n = 0
    try:
        last_page = driver.find_element_by_class_name('layui-laypage-last')
    except exceptions.NoSuchElementException:
        n = 1
    else:
        n = int(last_page.text)
    for i in range(n):
        for j in range(20):
            ref_path = '//*[@id="remark"]/dl[' + (j + 1).__str__() + ']'
            try:
                driver.find_element_by_xpath(ref_path)
            except exceptions.NoSuchElementException:
                break
            else:
                try:
                    title_box = driver.find_element_by_xpath(ref_path +
                                                             '/dt/a')
                    title = title_box.text
                    href = title_box.get_attribute('href')
                    authors = driver.find_element_by_xpath(
                        ref_path + '/dd[3]/span[1]').text
                    source = driver.find_element_by_xpath(
                        ref_path + '/dd[3]/span[2]/a').text
                    date = driver.find_element_by_xpath(ref_path +
                                                        '/dd[3]/span[3]').text
                except exceptions.NoSuchElementException:
                    print('没找到文献：' + ref_path)
                else:
                    authors = authors.replace('作者 ', '')
                    source = source.replace('《', '')
                    source = source.replace('》', '')
                    ref_dict = {
                        '标题': title,
                        '作者': authors,
                        '期刊名': source,
                        '出版日期': date,
                        '摘要或全文链接': href,
                        '搜索语句': search_phase
                    }
                    ref_list.append(ref_dict)
        if n == 1:
            break
        try:
            next_page = driver.find_element_by_class_name('layui-laypage-next')
        except exceptions.NoSuchElementException:
            print('没有下一页了')
            break
        else:
            next_page.click()
            time.sleep(2)
    return ref_list


def get_all_refs_ZW(org: str, mea: str, col_name: str):
    # 函数说明：
    # 功能：从知网空间爬取文献的链接
    # 参数：org：物种名，mea：繁殖力指标，集合名
    # 返回值：无
    db_name = "知网链接库"
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]

    if col_name in db.list_collection_names():
        col = db.get_collection(col_name)
    else:
        col = db[col_name]
    print(col.count_documents({}))
    link = "http://search.cnki.com.cn/"
    driver.get(link)
    # 输入物种名和繁殖力指标
    try:
        plus_btn = driver.find_element_by_xpath('/html/body/div[2]/div[2]/div/div[1]/a')
    except exceptions.NoSuchElementException:
        print("没找到plus_btn")
        return
    else:
        plus_btn.click()
        time.sleep(1)
    try:
        select_btn = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div/div[2]/div[1]/p')
    except exceptions.NoSuchElementException:
        print("没找到select_btn")
        return
    else:
        select_btn.click()
        time.sleep(1)
    try:
        sum_btn = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div/div[2]/div[1]/ul/li[7]')
    except exceptions.NoSuchElementException:
        print("没找到sum_btn")
        return
    else:
        sum_btn.click()
    try:
        plus_btn = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div/div[1]/a[2]')
    except exceptions.NoSuchElementException:
        print("没找到plus_btn2")
        return
    else:
        plus_btn.click()
        time.sleep(1)
    try:
        select_btn = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div[2]/div[2]/div[1]/p')
    except exceptions.NoSuchElementException:
        print("没找到select_btn2")
        return
    else:
        select_btn.click()
        time.sleep(1)
    try:
        sum_btn = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div[2]/div[2]/div[1]/ul/li[7]')
    except exceptions.NoSuchElementException:
        print("没找到sum_btn2")
        return
    else:
        sum_btn.click()
    try:
        input_area1 = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div[1]/div[2]/div[1]/input')
        input_area2 = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div[2]/div[2]/div[1]/input')
    except exceptions.NoSuchElementException:
        print("没找到input_area")
        return
    else:
        input_area1.send_keys(org)
        input_area2.send_keys(mea)
    # 点击搜索按钮
    try:
        search_btn = driver.find_element_by_xpath('/html/body/div[2]/div[3]/div[2]/div[2]/div[2]/a[1]')
    except exceptions.NoSuchElementException:
        print("没找到search_btn")
        return
    else:
        search_btn.click()
        time.sleep(3)
    # 循环获取文献信息
    have_result = False
    try:
        no_res = driver.find_element_by_xpath('/html/body/div[2]/div[1]/div[2]/div[11]/div[3]/div[1]/div[1]').text
    except exceptions.NoSuchElementException:
        print("=========have res===========")
        have_result = True
    else:
        if no_res != '抱歉！未检索到相关文献':
            have_result = True
    # title = '//*[@id="article_result"]/div/div[{0}]/p[1]/a[1]'
    # author = '//*[@id="article_result"]/div/div[{0}]/p[3]/span[1]'
    # source = '//*[@id="article_result"]/div/div[{0}]/p[3]/a[1]'
    page_count = 2
    page_path = '//*[@id="PageContent"]/div[1]/div[2]/div[13]/a[{0}]'
    while have_result:
        try:
            result = driver.find_element_by_xpath('//*[@id="article_result"]')
        except exceptions.NoSuchElementException:
            print("===========get result filed=========")
            break
        items = result.find_elements_by_class_name("list-item")
        for item in items:
            try:
                ti = item.find_element_by_class_name("left").text
                li = item.find_element_by_class_name("left").get_attribute("href")
            except exceptions.NoSuchElementException:
                print("==========no element========")
                break
            else:
                ref = {
                    "标题": ti,
                    "链接": li,
                    "指标": mea,
                    "物种": org,
                }
                if col.find_one(ref):
                    continue
                col.insert_one(ref)
        try:
            next_page = driver.find_element_by_xpath(page_path.format(page_count))
        except exceptions.NoSuchElementException:
            break
        else:
            if next_page.text == "下一页>" or next_page.text == "<上一页":
                break
            next_page.click()
            page_count += 1
            time.sleep(3)
    print("==========指标完！指标：{0}======".format(mea))


def get_refs_to_db(search_list: list, col_name: str, flag=1):
    # 函数说明：
    # 功能：获取文献，加入数据库中，如果库已经存在，更新库；如果库不存在，新建库并插入数据
    # 参数：db_name：库的名字；search_list：搜索列表，用于获取文献
    # 返回值：无
    if flag not in [1, 2, 3]:
        print("在get_refs_to_db中参数错误，flag的值应为1、2、3中的一个")
        return
    if flag == 1:
        db_name = '知网连接库'
    elif flag == 2:
        db_name = '万方'
    else:
        db_name = '维普'
    db_list = clint.list_database_names()
    if db_name in db_list:
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]
    col_list = db.list_collection_names()
    if col_name in col_list:
        col = db.get_collection(col_name)
    else:
        col = db[col_name]
    for phase in search_list:
        if flag == 1:
            rlist = get_all_refs_with_search_pahse(phase)
        elif flag == 2:
            rlist = get_all_refs_WF(phase)
        else:
            rlist = get_all_refs_WP(phase)
        if len(rlist) == 0:
            print(phase + " 没有找到文献")
            continue
        for ref in rlist:
            f = col.find_one(ref)
            if f is None:
                col.insert_one(ref)
        print(phase + "已搜索")


def get_refs_db_ZW(measure_list: list, org_name: str, col_name: str):
    # 函数说明：
    # 功能：从知网空间爬取文献的链接
    # 参数：measure_list：物种繁殖力指标名，org_name：物种名
    # 返回值：无
    iglist = ["性发育", "性欲", "交配能力", "繁殖周期", "发情表现", "精子畸形率", "精液品质", "射精量", "精子密度", "精子活力", "精子形态"]
    print("==========亚种开始！亚种名：{0}".format(org_name))
    for m in measure_list:
        if m in iglist:
            continue
        get_all_refs_ZW(org_name, m, col_name)
    print("=========亚种完! org name: {0}========".format(org_name))


# 获取种名和对应的亚种名构成的字典
ndiclist = get_names_dic('名称与指标.xlsx')
# 初始化webdriver
driver = webdriver.Chrome()
driver.implicitly_wait(3)
# 将获取到的文献写入数据库
clint = MongoClient()
# 1.中国知网，2.万方，3.维普
# f = 1
ignore_list = ["猪", "牛", "羊", "马", "猫", "狗", "鸡", "鹅", "鸭", "鸽子"]
for nameDic in ndiclist:
    oname = list(nameDic.keys())[0]
    if oname in ignore_list:
        continue
    mlist = get_measures(oname, '名称与指标.xlsx', 1)
    nlist = nameDic[oname]
    igList = []
    for name in nlist:
        if name in igList:
            continue
        get_refs_db_ZW(mlist, name, oname)
    print("=========物种完! org: {0}=========".format(oname))

""" for name in ndiclist:
    # 跳过已经搜过的物种
    search_phase_list = get_search_phase_list(name, flag=f)
    if isinstance(search_phase_list, type(None)):
        continue
    name_key = list(name.keys())[0]
    if name_key in ignore_list:
        continue
    get_refs_to_db(search_phase_list, name_key, flag=f) """

""" R=鸡 * 性成熟期已搜索
R=鸡 * 发情周期已搜索
R=鸡 * 发情持续期 没有找到文献
R=鸡 * 妊娠期已搜索"""
