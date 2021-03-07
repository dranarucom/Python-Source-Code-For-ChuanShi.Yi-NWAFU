import jieba
from pymongo import MongoClient
from selenium import webdriver
from selenium.common import exceptions
import time
import json


def read_all_ab_in_db(flag: int):
    # 函数说明：
    # 功能：从库中读取链接，获取文章的摘要
    # 参数：flag：1.知网；2.万方；3.维普
    # 返回值：无
    if flag == 1:
        db_name = '中国知网'
    elif flag == 2:
        db_name = '万方'
    elif flag == 3:
        db_name = '维普'
    else:
        print('flag参数值错误！')
        return
    db = clint.get_database(db_name)
    new_db_name = db_name + "摘要库"
    if new_db_name in clint.list_database_names():
        new_db = clint.get_database(new_db_name)
    else:
        new_db = clint[new_db_name]
    col_names = db.list_collection_names()
    from get_db_info import get_name_and_mearsure
    # 跳过已完成的物种
    ignore_list = ["马", "猫", '牛', '鸡', '鸭', '鸽子', '羊', '狗', '猪']
    for cname in col_names:
        if cname in ignore_list:
            continue
        print("物种‘" + cname + "'开始")
        col = db.get_collection(cname)
        if cname in new_db.list_collection_names():
            new_col = new_db.get_collection(cname)
            count = new_col.count_documents({})
        else:
            new_col = new_db[cname]
            count = 0
        i = 0
        with open("not_get_link.json", 'r', encoding='utf-8') as fl:
            jl = json.load(fl)
            if cname in jl.keys():
                count += len(jl[cname])
        # 如果没有爬到内容，保存链接，之后再统一处理
        not_get_link_list = []
        for r in col.find():
            # 跳过已经爬过的内容
            # test
            if i > 20:
                break
            i += 1
            sp = r['搜索语句']
            link = r['摘要或全文链接']
            smd = get_name_and_mearsure(sp, flag)
            if smd:
                s = smd['物种']
                m = smd['指标']
                if flag == 2:
                    ab = get_ab_wf(link)
                elif flag == 1:
                    ab = get_ab_zw(link)
                elif flag == 3:
                    ab = get_ab_wp(link)
                else:
                    return
                if ab:
                    rec = {'物种': s, '指标': m, '摘要': ab}
                    if new_col.find_one({'摘要': ab}):
                        continue
                    new_col.insert_one(rec)
                    continue
                print('"' + link + '",')
                not_get_link_list.append(link)
                continue
            print("没有得到物种和指标，语句：" + sp)
            break
        print("物种‘" + cname + "’已完成！")
        print("库中摘要数：" + str(new_col.count_documents({})))
        with open("not_get_link.json", 'r', encoding='utf-8') as fr:
            nlinks = json.load(fr)
            if cname in nlinks.keys():
                for nl in nlinks[cname]:
                    if nl in not_get_link_list:
                        continue
                    not_get_link_list.append(nl)
        with open("not_get_link.json", "w+", encoding="utf-8") as fw:
            nlinks[cname] = not_get_link_list
            json.dump(nlinks, fw)


def read_not_find_ab(flag: int):
    # 函数说明：
    # 功能：读取之前没有找到摘要的链接，试着再找一次
    # 参数：flag：1.知网；2.万方；3.维普
    # 返回值：无
    if flag == 1:
        db_name = '中国知网摘要库'
        sdb_name = '中国知网'
    elif flag == 2:
        db_name = '万方摘要库'
        sdb_name = '万方'
    elif flag == 3:
        db_name = '维普摘要库'
        sdb_name = '维普'
    else:
        print('flag参数值错误，没有对应的摘要库。')
        return
    db = clint.get_database(db_name)
    sdb = clint.get_database(sdb_name)
    from get_db_info import get_name_and_mearsure
    with open("not_get_link.json", "r", encoding="utf-8") as ngl_jr:
        ngl = json.load(ngl_jr)
        for k in ngl.keys():
            if k not in db.list_collection_names():
                col = db[k]
            else:
                col = db.get_collection(k)
            scol = sdb.get_collection(k)
            for link in ngl[k]:
                sr = scol.find_one({"摘要或全文链接": link})
                sp = sr['搜索语句']
                smd = get_name_and_mearsure(sp, flag)
                if smd:
                    s = smd['物种']
                    m = smd['指标']
                    if flag == 2:
                        ab = get_ab_wf(link)
                    elif flag == 1:
                        ab = get_ab_zw(link)
                    else:
                        ab = get_ab_wp(link)
                    if ab:
                        rec = {'物种': s, '指标': m, '摘要': ab}
                        if col.find_one({'摘要': ab}):
                            ngl[k].remove(link)
                            continue
                        col.insert_one(rec)
                        ngl[k].remove(link)
                        continue
                    print("not find link: " + link)
                    continue
                print("not find name, sp: " + sp)
            print(k + " over! current count: " + str(col.count_documents({})))
    with open("not_get_link.json", "w+", encoding="utf-8") as ngl_jw:
        json.dump(ngl, ngl_jw)
    print("read all not find link over!")


def get_ab_wf(link: str):
    # 函数说明：
    # 功能：根据链接获取文章的摘要，针对万方
    # 参数：link：文章链接
    # 返回值：摘要字符串
    try:
        driver.get(link)
    except exceptions.TimeoutException:
        return
    except exceptions.WebDriverException:
        return
    time.sleep(3)
    try:
        getmore_button = driver.find_element_by_class_name('getMore')
        getmore_button.click()
    except exceptions.ElementNotInteractableException:
        pass
    except exceptions.NoSuchElementException:
        pass
    try:
        ab = driver.find_element_by_xpath(
            '//*[@id="app"]/div[2]/div/div[2]/div[1]/div[3]/div[1]').text
    except exceptions.NoSuchElementException:
        return
    ab = ab.replace('收起∧', '')
    return ab


def get_ab_zw(link: str):
    """ try:
        driver.get(link)
    except exceptions.TimeoutException:
        return
    except exceptions.WebDriverException:
        return
    time.sleep(3) """


def get_ab_wp(link: str):
    # 函数说明：
    # 功能：根据链接获取文章的摘要，针对维普
    # 参数：link：文章链接
    # 返回值：摘要字符串
    try:
        driver.get(link)
    except exceptions.TimeoutException:
        return
    except exceptions.WebDriverException:
        return
    time.sleep(3)
    try:
        ab = driver.find_element_by_xpath('//*[@id="body"]/div/div/div[1]/div[4]/div[1]/span[2]/span').text
    except exceptions.NoSuchElementException:
        return
    else:
        return ab


def get_gene_ab(flag: int):
    # 函数说明：
    # 功能：读取摘要，分词，得到其中的基因信息
    # 参数：flag：1.知网；2.万方；3.维普
    # 返回值：无
    if flag == 1:
        db_name = "知网摘要库"
        new_db_name = "知网基因库"
    elif flag == 2:
        db_name = "万方摘要库"
        new_db_name = "万方基因库"
    elif flag == 3:
        db_name = "维普摘要库"
        new_db_name = "维普基因库"
    else:
        print("输入参数错误，方法：get_gene_ab")
        return
    db = clint.get_database(db_name)
    if new_db_name in clint.list_database_names():
        new_db = clint.get_database(new_db_name)
    else:
        new_db = clint[new_db_name]
    for cname in db.list_collection_names():
        col = db.get_collection(cname)
        if cname in new_db.list_collection_names():
            new_col = new_db.get_collection(cname)
        else:
            new_col = new_db[cname]
        find_gene_list = []
        for abr in col.find():
            ab = jieba.lcut(abr["摘要"])
            for c in ab:
                if ('a' <= c[0] <= 'z') or ('A' <= c[0] <= 'Z'):
                    if search_gene_list(c, abr['物种']):
                        find_gene_list.append({
                            "物种": abr["物种"],
                            "指标": abr["指标"],
                            "基因": c
                        })
        for gr in find_gene_list:
            new_col.insert_one(gr)


def search_gene_list(phase: str, s: str):
    # 函数说明：
    # 功能：将一个英文词在基因列表中比对，返回比对的结果
    # 参数：phase: 一个英文词；s：物种名
    # 返回值：比对成功返回True，否则返回False
    return False


clint = MongoClient()
driver = webdriver.Chrome()
driver.implicitly_wait(3)
f = 1
read_all_ab_in_db(flag=f)
# read_not_find_ab(f)

"""
万方：
鹅 over! current count: 1751
鸡 over! current count: 10026
鸽子 over! current count: 54
猫 over! current count: 275
狗 over! current count: 1513
马 over! current count: 4378
鸭 over! current count: 2416
猪 over! current count: 13993
牛 over! current count: 21131
羊 over! current count: 15535
63行

维普：
over!
"""

""" db = clint.get_database('万方摘要库')
col = db.get_collection('鹅')
print(col.count_documents({})) """
""" print(sr)
print(type(sr))
"""
"""
分词输出：
s1 = 'CPU是E3 1231 V3显卡是1050ti8G内存条我之前是780ti假卡，512MB的显存，卡的要死，换个1050ti应该可以了吧？'
sr = jieba.lcut(s1)
print(sr)
print(type(sr))

结果：
['CPU', '是', 'E3', ' ', '1231', ' ', 'V3', '显卡', '是', '1050ti8G', '内存条', '我', '之前', '是', '780ti', '假卡', '，', '512MB', '的', '显存', '，', '卡', '的',
'要死', '，', '换个', '1050ti', '应该', '可以', '了', '吧', '？']
<class 'list'>
"""
"""
输出英文词：
s1 = 'CPU是E3 1231 V3显卡是1050ti8G内存条我之前是780ti假卡，512MB的显存，卡的要死，换个1050ti应该可以了吧？'
sr = jieba.lcut(s1)
for s in sr:
    if ('a' < s[0] < 'z') or ('A' < s[0] < 'Z'):
        print(s)

结果：
CPU
E3
V3
"""
