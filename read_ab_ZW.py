from pymongo import MongoClient
from selenium import webdriver
from selenium.common import exceptions
import time
import json


def read_ab_ZW(link: str):
    ab_xpath = '//*[@id="content"]/div[2]/div[4]'
    try:
        driver.get(link)
    except exceptions.TimeoutException:
        time.sleep(3)
        return
    time.sleep(3)
    try:
        ab = driver.find_element_by_xpath(ab_xpath).text
    except exceptions.NoSuchElementException:
        return
    else:
        return ab


def get_ref_db():
    db_name = "知网链接库"
    new_db_name = "知网摘要库"
    res_db_name = "知网摘要结果暂存库"
    # 链接数据库
    if db_name in client.list_database_names():
        db = client.get_database(db_name)
    else:
        db = client[db_name]
    if new_db_name in client.list_database_names():
        new_db = client.get_database(new_db_name)
    else:
        new_db = client[new_db_name]
    if res_db_name in client.list_database_names():
        res_db = client.get_database(res_db_name)
    else:
        res_db = client[res_db_name]
    # 读取链接，爬取摘要
    ig_col_list = ["马", "羊", "鸡", "鹅", "牛", "鸭", "猫", "狗", "鸽子", "猪"]
    for col_name in db.list_collection_names():
        if col_name in ig_col_list:
            continue
        print("=========物种-{0}-开始===========".format(col_name))
        col = db.get_collection(col_name)
        if col_name in new_db.list_collection_names():
            new_col = new_db.get_collection(col_name)
        else:
            new_col = new_db[col_name]
        if col_name in res_db.list_collection_names():
            res_col = res_db.get_collection(col_name)
        else:
            res_col = res_db[col_name]
        not_get_list = []
        for ref in col.find():
            org = ref["物种"]
            mea = ref["指标"]
            link = ref["链接"]
            if res_col.find_one({"链接": link}):
                continue
            ab = read_ab_ZW(link)
            if ab:
                rec = {
                    "物种": org,
                    "指标": mea,
                    "摘要": ab
                }
                if new_col.find_one(rec):
                    continue
                new_col.insert_one(rec)
                res_col.insert_one({"链接": link})
                continue
            print("not find link: {0},".format(link))
            not_get_list.append(link)
        # 读取现有的失败链接，加入字典
        with open("not_get_link.json", 'r', encoding='utf-8') as fr:
            nlinks = json.load(fr)
            if col_name in nlinks.keys():
                for nl in nlinks[col_name]:
                    if nl in not_get_list:
                        continue
                    not_get_list.append(nl)
        with open("not_get_link.json", "w+", encoding="utf-8") as fw:
            nlinks[col_name] = not_get_list
            json.dump(nlinks, fw)
        print("========物种-{0}-结束，获得数量：{1}===========".format(col_name, new_col.count_documents({})))


def read_not_find_link_ZW():
    file_name = "not_get_link.json"
    db_name = "知网链接库"
    new_db_name = "知网摘要库"
    if db_name in client.list_database_names():
        db = client.get_database(db_name)
    else:
        print("==============没找到库：{0}=============".format(db_name))
        return
    if new_db_name in client.list_database_names():
        new_db = client.get_database(new_db_name)
    else:
        new_db = client[new_db_name]
    with open(file_name, "r", encoding="utf-8") as fr:
        link_dic = json.load(fr)
        for org in link_dic.keys():
            if org in db.list_collection_names():
                col = db.get_collection(org)
            else:
                print("========没找到集合：{0}==============".format(org))
                continue
            if org in new_db.list_collection_names():
                new_col = new_db.get_collection(org)
            else:
                new_col = new_db[org]
            print("============物种：{0}开始==========".format(org))
            link_list = link_dic[org]
            new_link_list = []
            if len(link_list) == 0:
                continue
            for link in link_list:
                rec = col.find_one({"链接": link})
                if rec:
                    oname = rec["物种"]
                    mname = rec["指标"]
                    ab = read_ab_ZW(link)
                    if ab:
                        ref = {
                            "物种": oname,
                            "指标": mname,
                            "摘要": ab
                        }
                        if new_col.find_one(ref):
                            continue
                        new_col.insert_one(ref)
                    new_link_list.append(link)
            link_dic[org] = new_link_list
            print("============物种：{0}结束========".format(org))
    with open(file_name, "w+", encoding="utf-8") as fw:
        json.dump(link_dic, fw)


client = MongoClient()
driver = webdriver.Chrome()
driver.implicitly_wait(3)
# get_ref_db()
read_not_find_link_ZW()


# 36行
"""
{
    "\u9a6c": [],
    "\u7f8a": [
        "https://cdmd.cnki.com.cn/Article/CDMD-10129-1020147061.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-XMYH202003041.htm",
        "https://cdmd.cnki.com.cn/Article/CDMD-10712-1020645880.htm",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001043&dbname=IPFDLAST2017",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-ZYSS201906026.htm",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001012&dbname=IPFDLAST2017",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001013&dbname=IPFDLAST2017"
    ],
    "\u9e21": [
        "https://www.cnki.com.cn/Article/CJFDTOTAL-ZGDB199908047.htm"
    ],
    "\u9e45": [
        "https://www.cnki.com.cn/Article/CJFDTOTAL-FZJH201317012.htm"
    ],
    "\u725b": [
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001010&dbname=IPFDLAST2017",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001020&dbname=IPFDLAST2017",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001189&dbname=IPFDLAST2017",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-KXZH201103050.htm",
        "https://cdmd.cnki.com.cn/Article/CDMD-10504-1020366252.htm",
        "https://cdmd.cnki.com.cn/Article/CDMD-10712-1020907445.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-GATE201506033.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-DYJZ200403021.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-ZGNN201705002.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-HLJX201619031.htm",
        "https://cdmd.cnki.com.cn/Article/CDMD-10759-1011204220.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-HLJD201303009.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-XMZY201411042.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-XMSK201707042.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-XMZY201111033.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-BULL198703027.htm"
    ],
    "\u9e2d": [
        "https://www.cnki.com.cn/Article/CJFDTOTAL-NJKK201008053.htm",
        "https://cdmd.cnki.com.cn/Article/CDMD-10712-1020645767.htm",
        "https://www.cnki.com.cn/Article/CJFDTOTAL-ZXWA200606040.htm"
    ],
    "\u732b": [],
    "\u72d7": [
        "https://www.cnki.com.cn/Article/CJFDTOTAL-YQZZ200401002.htm"
    ],
    "\u9e3d\u5b50": [],
    "\u732a": [
        "https://cdmd.cnki.com.cn/Article/CDMD-10712-1020645742.htm",
        "https://cdmd.cnki.com.cn/Article/CDMD-10712-1020907191.htm",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGXJ201608001171&dbname=IPFDLAST2017",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=IGSQ200809001195&dbname=IPFD9914",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZXMY201505001027&dbname=IPFDLAST2015",
        "https://www.cnki.net/KCMS/detail/detail.aspx?dbcode=IPFD&filename=ZGYL201207001111&dbname=IPFD9914"
    ]
}
"""
