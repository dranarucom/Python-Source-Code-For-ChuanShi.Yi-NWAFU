import json
from Bio import Entrez
from pymongo import MongoClient
import jieba
import time


def get_all_org_name():
    # 函数说明：
    # 功能：从文件中读取物种的所有中文和英文名称
    # 参数：
    # 返回值：以中文名为键，对应英文名为值的字典
    with open("家畜名称.json", "r", encoding="utf-8") as org_name_file:
        names = json.load(org_name_file)
        return names


def get_eng_name(chi_name: str):
    # 函数说明：
    # 功能：输入物种中文名，返回对应的英文名
    # 参数：chi_name：物种中文名
    # 返回值：物种的英文名
    if chi_name in org_names.keys():
        return org_names[chi_name]
    else:
        if chi_name in name_chg_dic.keys():
            return org_names[name_chg_dic[chi_name]]
        else:
            return None


def try_find_in_gene(term_to_find: str, org_name: str):
    # 函数说明：
    # 功能：在NCBI的gene库中检索，看某个物种中是否存在某个基因
    # 参数：term_to_find：待检索的基因名；org_name：物种英文名
    # 返回值：检索结果为存在，返回True；否则返回False
    search_term = str.format("({0}[Gene Name]) AND {1}[Organism]", term_to_find, org_name)
    try:
        h = Entrez.esearch(db="gene", term=search_term)
        hr = Entrez.read(h)
    except Exception:
        print("========search {0} failed=========".format(search_term))
        return False
    else:
        count = int(hr['Count'])
        if count > 0:
            return True
        else:
            return False


def try_search_ab(ab: str, eng_name: str):
    # 函数说明：
    # 功能：搜索一个摘要，找出其中的基因名称
    # 参数：ab：摘要字符串，eng_name：摘要对应物种的英文名
    # 返回值：基因名称的列表
    result_gene_list = []
    search_term_list = []
    # 摘要分词，找出其中包含英文的词
    word_list = jieba.lcut(ab)
    for word in word_list:
        if ('a' < word[0] < 'z') or ('A' < word[0] < 'Z'):
            search_term_list.append(word)
    if len(search_term_list) == 0:
        return None
    # 在NCBI检索每一个英文词，看它是否是物种的一个基因
    for term in search_term_list:
        if try_find_in_gene(term, eng_name):
            result_gene_list.append(term)
            continue
        return None

    return result_gene_list


def search_all_ab_dbs(flag: int):
    # 函数说明：
    # 功能：搜寻一个摘要库
    # 参数：flag：1.知网，2.万方，3.维普
    # 返回值：无
    if flag == 1:
        db_name = "中国知网摘要库"
        new_db_name = "知网基因库"
        res_db_name = "知网摘要结果暂存"
    elif flag == 2:
        db_name = "万方摘要库"
        new_db_name = "万方基因库"
        res_db_name = "万方摘要结果暂存"
    elif flag == 3:
        db_name = "维普摘要库"
        new_db_name = "维普基因库"
        res_db_name = "维普摘要结果暂存"
    else:
        print("=======flag值错误，flag应等于1或2或3，现在flag={0}======".format(flag))
        return

    if db_name not in client.list_database_names():
        print("=====没有 {0} 这个库========".format(db_name))
        return

    db = client.get_database(db_name)
    if new_db_name in client.list_database_names():
        new_db = client.get_database(new_db_name)
    else:
        new_db = client[new_db_name]
    if res_db_name in client.list_database_names():
        res_db = client.get_database(res_db_name)
    else:
        res_db = client[res_db_name]

    ignore_cols = ["鹅", "猫", "鸡", "鸭", "马", "狗"]    # 在gene中搜不到鹅的基因，所以跳过
    for col_name in db.list_collection_names():
        if col_name in ignore_cols:
            continue
        col = db.get_collection(col_name)

        if col_name in new_db.list_collection_names():
            new_col = new_db.get_collection(col_name)
        else:
            new_col = new_db[col_name]
        if col_name in res_db.list_collection_names():
            res_col = res_db.get_collection(col_name)
        else:
            res_col = res_db.get_collection(col_name)

        count = res_col.count_documents({})
        cursor = col.find(no_cursor_timeout=True)
        for ref in cursor:
            if res_col.find_one({"已完成id": ref["_id"]}):
                continue
            ab = ref["摘要"]
            cname = ref["物种"]
            m = ref["指标"]
            ename = get_eng_name(cname)
            if ename is not None:
                glist = try_search_ab(ab, ename)
                if glist:
                    for g in glist:
                        rec = {
                            "物种": cname,
                            "指标": m,
                            "基因": g
                        }
                        if new_col.find_one(rec):
                            continue
                        new_col.insert_one(rec)
                    res_col.insert_one({"已完成id": ref["_id"]})
                    time.sleep(2)
                    count += 1
                    if (count % 500) == 0:
                        print("=====search over! count: {0}======".format(count))
            else:
                print("======没有找到【{0}】对应的英文=====".format(cname))
        print("========search db over! name: {0}，count: {1}, get: {2}=========".format(col_name, count, new_col.count_documents({})))
        cursor.close()


Entrez.email = "1327133922@qq.com"
client = MongoClient()
org_names = get_all_org_name()
name_chg_dic = {"奶牛": "家牛", "鸡": "家鸡", "牛": "家牛", "鸭": "绿头鸭", "羊": "绵羊", "犬": "狗", "黄牛": "家牛"}
f = 3
search_all_ab_dbs(f)

# 鸡：已完成：10180；得到：272
# 鸭：已完成：1691；得到：32
# 马：已完成：614；得到：61
# 狗：已完成：797；得到：71
# 牛：已完成：4571；得到：280

""" cname = "水牛"
# test_terms = ["SPP1", "CSN3", "TNF", "LTF", "fullmetal", "count3", "dome"]
ename = get_eng_name(cname)
test_ab = "发绿色A那你们呢，快速CSN饿啊。发了工资用于血常规，很痛苦年七月，你可能有人fullmetal页面载入。爱人，鹅空气过count3艾玛斯通，有人水稻为dome室内地图。"
if ename is not None:
    re = try_search_ab(test_ab, ename)
    if len(re) > 0:
        print(re) """
