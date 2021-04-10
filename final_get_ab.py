from pymongo import MongoClient
from selenium import webdriver
from selenium.common import exceptions
import time
import json


def read_all_ab_in_db(db_name: str, new_db_name: str, flag: int):
    # 函数说明：
    # 功能：从库中读取链接，获取文章的摘要
    # 参数：flag：1.知网；2.万方；3.维普
    # 返回值：无
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        print("没有【{0}】这个库".format(db_name))
        return
    if new_db_name in clint.list_database_names():
        new_db = clint.get_database(new_db_name)
    else:
        new_db = clint[new_db_name]
    if flag not in [1, 2, 3]:
        print("flag的值应为1、2、3中的一个，现在flag的值为：{0}".format(flag))
        return
    # 跳过已完成的物种
    ignore_list = []
    for cname in db.list_collection_names():
        if cname in name_chg_dic.keys() and name_chg_dic[cname] in db.list_collection_names():
            cname = name_chg_dic[cname]
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

        with open("not_get_link.json", 'r', encoding='utf-8') as fl:
            jl = json.load(fl)
            if cname in jl.keys():
                count += len(jl[cname])
        # 如果没有爬到内容，保存链接，之后再统一处理
        not_get_link_list = []
        cursor = col.find(no_cursor_timeout=True)
        for r in cursor:
            link = r['链接']
            s = r["物种"]
            m = r["指标"]
            title = r["标题"]
            if flag == 2:
                # ab = get_ab_wf(link)
                return
            elif flag == 1:
                ab = get_ab_zw(link)
            elif flag == 3:
                ab = get_ab_wp(link)
            else:
                return
            if ab:
                rec = {'标题': title, '物种': s, '指标': m, '摘要': ab}
                if new_col.find_one({'摘要': ab}):
                    continue
                new_col.insert_one(rec)
                continue
            print('"{0}",'.format(link))
            not_get_link_list.append(link)
            continue
        cursor.close()
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


def read_not_find_ab(db_name: str, new_db_name: str, flag: int):
    # 函数说明：
    # 功能：读取之前没有找到摘要的链接，试着再找一次
    # 参数：flag：1.知网；2.万方；3.维普
    # 返回值：无
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        print("没有【{0}】这个库".format(db_name))
        return
    if new_db_name in clint.list_database_names():
        new_db = clint.get_database(new_db_name)
    else:
        new_db = clint[new_db_name]
    if flag not in [1, 2, 3]:
        print("flag的值应为1、2、3中的一个，现在flag的值为：{0}".format(flag))
        return

    with open("not_get_link.json", "r", encoding="utf-8") as ngl_jr:
        ngl = json.load(ngl_jr)
        for k in ngl.keys():
            if k not in new_db.list_collection_names():
                new_col = new_db[k]
            else:
                new_col = new_db.get_collection(k)
            scol = db.get_collection(k)
            for link in ngl[k]:
                sr = scol.find_one({"链接": link})
                if not sr:
                    print("没有找到链接，库：{0}，集合：{1}，链接：{2}".format(db_name, k, link))
                    continue
                s = sr["物种"]
                m = sr["指标"]
                title = sr["标题"]
                if flag == 2:
                    # ab = get_ab_wf(link)
                    return
                elif flag == 1:
                    ab = get_ab_zw(link)
                elif flag == 3:
                    ab = get_ab_wp(link)
                else:
                    return
                if ab:
                    ngl[k].remove(link)
                    rec = {'标题': title, '物种': s, '指标': m, '摘要': ab}
                    if new_col.find_one(rec):
                        continue
                    new_col.insert_one(rec)
                    continue
                print('"{0}",'.format(link))
                continue
            print(k + " over! current count: " + str(new_col.count_documents({})))
    with open("not_get_link.json", "w+", encoding="utf-8") as ngl_jw:
        json.dump(ngl, ngl_jw)
    print("read all not find link over!")


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
    time.sleep(5)
    try:
        ab = driver.find_element_by_xpath('//*[@id="body"]/div/div/div[1]/div[4]/div[1]/span[2]/span').text
    except exceptions.NoSuchElementException:
        return
    else:
        return ab


def get_ab_zw(link: str):
    ab_xpath = '//*[@id="content"]/div[2]/div[4]'
    sub_ab_xpath = '//*[@id="ChDivSummary"]'
    try:
        driver.get(link)
    except exceptions.TimeoutException:
        time.sleep(3)
        return
    time.sleep(5)
    try:
        ab = driver.find_element_by_xpath(ab_xpath).text
    except exceptions.NoSuchElementException:
        try:
            ab = driver.find_element_by_xpath(sub_ab_xpath).text
        except exceptions.NoSuchElementException:
            return
        else:
            return ab
    else:
        return ab


clint = MongoClient()
driver = webdriver.Chrome()
driver.implicitly_wait(3)
name_chg_dic = {"犬": "狗"}
f = 3
dname = "新维普链接库"
nname = "新维普摘要库"
# read_all_ab_in_db(dname, nname, f)
read_not_find_ab(dname, nname, f)
