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


def get_all_refs_WF(org: str, mea: str, db_name: str):
    # 函数说明：
    # 功能：用指定的检索语句在万方检索，获取所有的检索结果
    # 参数：search_phase：检索语句（摘要:("名称"*"繁殖力指标")）
    # 返回值：文献信息字典组成的列表
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]

    if org in db.list_collection_names():
        col = db.get_collection(org)
    else:
        col = db[org]

    link = 'http://new.wanfangdata.com.cn/searchResult/getAdvancedSearch.do?searchType=all'
    search_phase = '摘要:"{0}" and "{1}"'.format(org, mea)
    driver.get(link)
    # 等待网页内容加载完成
    time.sleep(1)
    # 点击专业搜索按钮
    try:
        pro_button = driver.find_element_by_xpath('/html/body/div[4]/div/div[1]/div[2]/div[1]/span[2]')
    except exceptions.NoSuchElementException:
        print('没有找到专业搜索按钮')
        return None
    else:
        # pro_button.click()
        driver.execute_script("arguments[0].click();", pro_button)
    time.sleep(1)
    # 输入检索条件
    try:
        text_area = driver.find_element_by_xpath(
            '/html/body/div[4]/div/div[1]/div[4]/div[2]/div[2]/div[3]/div[1]/div[1]/textarea')
    except exceptions.NoSuchElementException:
        print('没有找到输入检索条件的区域')
        return None
    else:
        text_area.send_keys(search_phase)
    try:
        search_button = driver.find_element_by_xpath('/html/body/div[4]/div/div[1]/div[6]/span[1]')
    except exceptions.NoSuchElementException:
        print('没有找到搜索按钮')
        return None
    else:
        # search_button.click()
        driver.execute_script("arguments[0].click();", search_button)
    time.sleep(5)
    page_count = 2
    while True:
        for i in range(1, 21):
            title_xpath = '/html/body/div[4]/div/div[2]/div[2]/div[2]/div[2]/div[2]/div/div/div[3]/div[{0}]/div/div[1]/div[2]/span[2]'.format(i)
            ab_xpath = '/html/body/div[4]/div/div[2]/div[2]/div[2]/div[2]/div[2]/div/div/div[3]/div[{0}]/div/div[3]/span[2]'.format(i)
            try:
                title = driver.find_element_by_xpath(title_xpath).text
                ab = driver.find_element_by_xpath(ab_xpath).text
            except exceptions.NoSuchElementException:
                continue
            except exceptions.StaleElementReferenceException:
                title = driver.find_element_by_xpath(title_xpath).text
                ab = driver.find_element_by_xpath(ab_xpath).text
                ref_dict = {
                    "标题": title,
                    "摘要": ab,
                    "指标": mea,
                    "物种": org,
                }
                if col.find_one(ref_dict):
                    continue
                col.insert_one(ref_dict)
            else:
                # 构建文献字典
                ref_dict = {
                    "标题": title,
                    "摘要": ab,
                    "指标": mea,
                    "物种": org,
                }
                if col.find_one(ref_dict):
                    continue
                col.insert_one(ref_dict)
        try:
            next_page_xpath = '/html/body/div[4]/div/div[2]/div[2]/div[2]/div[2]/div[2]/div/div/div[5]/span[{0}]'.format(page_count + 1)
            next_page_button = driver.find_element_by_xpath(next_page_xpath)
        except exceptions.NoSuchElementException:
            print("没有找到下一页按钮")
            break
        else:
            # next_page_button.click()
            if next_page_button.get_attribute('class') == 'next':
                break
            page_count += 1
            driver.execute_script("arguments[0].click();", next_page_button)
            time.sleep(10)


def get_all_refs_WP(org: str, mea: str, db_name: str):
    # 函数说明：
    # 功能：用指定的检索语句在维普检索，获取所有的检索结果
    # 参数：search_phase：检索语句（R=品种名 * 繁殖力指标）
    # 返回值：文献信息字典组成的列表
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]

    if org in db.list_collection_names():
        col = db.get_collection(org)
    else:
        col = db[org]

    link = 'http://qikan.cqvip.com/Qikan/Search/Advance?from=index'
    search_phase = 'R={0}*{1}'.format(org, mea)
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
    time.sleep(10)
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
                except exceptions.NoSuchElementException:
                    print('没找到文献：' + ref_path)
                except exceptions.StaleElementReferenceException:
                    title_box = driver.find_element_by_xpath(ref_path +
                                                             '/dt/a')
                    title = title_box.text
                    href = title_box.get_attribute('href')
                    ref_dict = {
                        "标题": title,
                        "链接": href,
                        "指标": mea,
                        "物种": org,
                    }
                    if col.find_one(ref_dict):
                        continue
                    col.insert_one(ref_dict)
                else:
                    ref_dict = {
                        "标题": title,
                        "链接": href,
                        "指标": mea,
                        "物种": org,
                    }
                    if col.find_one(ref_dict):
                        continue
                    col.insert_one(ref_dict)
        if n == 1:
            break
        try:
            next_page = driver.find_element_by_class_name('layui-laypage-next')
        except exceptions.NoSuchElementException:
            print('没有下一页了')
            break
        else:
            driver.execute_script("arguments[0].click();", next_page)
            time.sleep(5)
    print("==========指标完！指标：{0}======".format(mea))


def get_all_refs_ZW(org: str, mea: str, db_name: str):
    # 函数说明：
    # 功能：从知网空间爬取文献的链接
    # 参数：org：物种名，mea：繁殖力指标，集合名
    # 返回值：无
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]

    if org in db.list_collection_names():
        col = db.get_collection(org)
    else:
        col = db[org]
    print(col.count_documents({}))
    link = "http://search.cnki.com.cn/"
    driver.get(link)
    # 输入物种名和繁殖力指标
    try:
        plus_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[2]/div/div[1]/a')
    except exceptions.NoSuchElementException:
        print("没找到plus_btn")
        return
    else:
        plus_btn.click()
        time.sleep(1)
    try:
        select_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div/div[2]/div[1]/p')
    except exceptions.NoSuchElementException:
        print("没找到select_btn")
        return
    else:
        select_btn.click()
        time.sleep(1)
    try:
        sum_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div/div[2]/div[1]/ul/li[7]')
    except exceptions.NoSuchElementException:
        print("没找到sum_btn")
        return
    else:
        sum_btn.click()
    try:
        plus_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div/div[1]/a[2]')
    except exceptions.NoSuchElementException:
        print("没找到plus_btn2")
        return
    else:
        plus_btn.click()
        time.sleep(1)
    try:
        select_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div[2]/div[2]/div[1]/p')
    except exceptions.NoSuchElementException:
        print("没找到select_btn2")
        return
    else:
        select_btn.click()
        time.sleep(1)
    try:
        sum_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div[2]/div[2]/div[1]/ul/li[7]')
    except exceptions.NoSuchElementException:
        print("没找到sum_btn2")
        return
    else:
        sum_btn.click()
    try:
        input_area1 = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div[1]/div[2]/div[1]/input')
        input_area2 = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div[2]/div[2]/div[1]/input')
    except exceptions.NoSuchElementException:
        print("没找到input_area")
        return
    else:
        input_area1.send_keys(org)
        input_area2.send_keys(mea)
    # 点击搜索按钮
    try:
        search_btn = driver.find_element_by_xpath(
            '/html/body/div[2]/div[3]/div[2]/div[2]/div[2]/a[1]')
    except exceptions.NoSuchElementException:
        print("没找到search_btn")
        return
    else:
        search_btn.click()
        time.sleep(3)
    # 循环获取文献信息
    have_result = False
    try:
        no_res = driver.find_element_by_xpath(
            '/html/body/div[2]/div[1]/div[2]/div[11]/div[3]/div[1]/div[1]'
        ).text
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
                li = item.find_element_by_class_name("left").get_attribute(
                    "href")
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
            next_page = driver.find_element_by_xpath(
                page_path.format(page_count))
        except exceptions.NoSuchElementException:
            break
        else:
            if next_page.text == "下一页>" or next_page.text == "<上一页":
                break
            next_page.click()
            page_count += 1
            time.sleep(3)
    print("==========指标完！指标：{0}======".format(mea))


def get_refs_to_db(measure_list: list, org_name: str, flag=1):
    # 函数说明：
    # 功能：获取文献，加入数据库中，如果库已经存在，更新库；如果库不存在，新建库并插入数据
    # 参数：measure_list：物种繁殖力指标名，org_name：物种名，flag：1.知网，2.万方，3.维普
    # 返回值：无
    if flag not in [1, 2, 3]:
        print("在get_refs_to_db中参数错误，flag的值应为1、2、3中的一个")
        return
    if flag == 1:
        db_name = '新知网连接库'
    elif flag == 2:
        db_name = '新万方摘要库'
    else:
        db_name = '新维普链接库'

    for mea in measure_list:
        if flag == 1:
            get_all_refs_ZW(org_name, mea, db_name)
        elif flag == 2:
            get_all_refs_WF(org_name, mea, db_name)
        else:
            get_all_refs_WP(org_name, mea, db_name)
        print("物种【{0}】的【{1}】结束".format(org_name, mea))


# 获取种名和对应的亚种名构成的字典
ndiclist = {
    "猪": ["猪"],
    "牛": ["奶牛", "水牛", "牦牛"],
    "羊": ["山羊", "绵羊"],
    "马": ["马"],
    "猫": ["猫"],
    "狗": ["犬", "狗"],
    "鸡": ["火鸡", "家鸡"],
    "鸭": ["绿头鸭", "疣鼻栖鸭"]
}
# 初始化webdriver
driver = webdriver.Chrome()
driver.implicitly_wait(3)
driver.set_page_load_timeout(60)
# 将获取到的文献写入数据库
clint = MongoClient()
# 1.中国知网，2.万方，3.维普
for f in [2]:
    print("============第【{0}】个库开始=============".format(f))
    ignore_list = []
    for oname in ndiclist.keys():
        if oname in ignore_list:
            continue
        mlist = get_measures(oname, '名称与指标.xlsx', 1)
        nlist = ndiclist[oname]
        igList = ["奶牛", "水牛"]
        for name in nlist:
            if name in igList:
                continue
            get_refs_to_db(mlist, name, f)
    print("============第【{0}】个库结束=============".format(f))
