from pymongo import MongoClient
from get_db_info import get_name_and_mearsure
import json


def count_ref_with_org(db_flag=1):
    # 功能：统计某个库中各个物种文献的数量
    # 参数：db_flag：1.知网，2.万方，3.维普
    # 返回值：无
    if db_flag == 1:
        db_name = "中国知网"
    elif db_flag == 2:
        db_name = "万方摘要库"
    elif db_flag == 3:
        db_name = "维普摘要库"
    else:
        print("db flag can't be this value: " + str(db_flag))
        return
    # 连接数据库，使用字典记录统计记录
    db = clint.get_database(db_name)
    count_result = dict()
    count_result["库名"] = db_name
    for col_name in db.list_collection_names():
        col = db.get_collection(col_name)
        for ref in col.find():
            # 获取种名
            if '物种' in ref.keys():
                orga_name = ref['物种']
            elif '搜索语句' in ref.keys():
                orga_name = get_name_and_mearsure(ref['搜索语句'], db_flag)['物种']
            else:
                print("Can not get org name in ref: ")
                print(ref)
                break
            if orga_name not in count_result.keys():
                count_result[orga_name] = 1
            else:
                count_result[orga_name] += 1
    print(count_result)
    file_name = db_name + "物种数量统计.json"
    with open(file_name, "w+", encoding="utf-8") as out_file:
        json.dump(count_result, out_file)
        print("dump over, db name: " + db_name)


def count_all_refs_db(db_name: str):
    # 功能：统计某个库中文献的数量
    # 参数：db_name: 库名
    # 返回值：文献数量
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        print("没有找到名为【{0}】的库".format(db_name))
        return
    col_list = db.list_collection_names()
    ref_count = 0
    for col in col_list:
        collection = db.get_collection(col)
        col_count = collection.count_documents({})
        ref_count += col_count
    print(ref_count)


def get_att_name(col_name: str, att: str, db_flag=1):
    # 功能：获取某个集合中某个文献信息的种类，比如来源期刊或机构的种类
    # 参数：col_name：集合名，att：需要获取的文献信息，db_flag：1.中国知网，2.万方，3.维普
    # 返回值：文献信息种类的列表
    if db_flag == 1:
        db_name = '中国知网'
    elif db_flag == 2:
        db_name = '万方'
    elif db_flag == 3:
        db_name = '维普'
    else:
        print('db_flag参数错误')
        return
    db = clint.get_database(db_name)
    col_list = db.list_collection_names()
    if col_name not in col_list:
        print(db_name + '库中没有' + col_name + '这个集合')
        return
    col = db.get_collection(col_name)
    i = 0
    keys = list()
    for r in col.find():
        if i == 0:
            keys = list(r.keys())
            i += 1
    if att not in keys:
        print(db_name + '库的' + col_name + '集合中没有' + att + '这个属性')
        return
    result_list = []
    for r in col.find():
        if r[att] in result_list:
            continue
        result_list.append(r[att])
    return result_list


def count_att_num(col_name: str, att_name: str, att: str, db_flag=1):
    # 功能：统计某个库某个集合中某类文献的数量，比如来源于某个期刊的文献数量
    # 参数：参数：col_name：集合名，att_name：需要获取的文献信息，db_flag：1.中国知网，2.万方，3.维普
    # 返回值：属性名：数量的字典列表
    if db_flag == 1:
        db_name = '中国知网'
    elif db_flag == 2:
        db_name = '万方'
    elif db_flag == 3:
        db_name = '维普'
    else:
        print('db_flag参数错误')
        return
    db = clint.get_database(db_name)
    col_list = db.list_collection_names()
    if col_name not in col_list:
        print(db_name + '库中没有' + col_name + '这个集合')
        return
    col = db.get_collection(col_name)
    return dict({att: att_name, '数量': col.count_documents({att: att_name})})


def count_repeat_num(src_db_name: str, other_db_name: str):
    # 函数说明：
    # 功能：计算当前库的一个集合中有多少篇文献在目标库中重复
    # 参数：src_db_name：当前库名；other_db_name：目标库名
    # 返回值：无
    if src_db_name == other_db_name:
        print("两个库不能是同一个库！")
        return
    if src_db_name in clint.list_database_names():
        src_db = clint.get_database(src_db_name)
    else:
        print("没有【{0}】这个库".format(src_db_name))
        return
    if other_db_name in clint.list_database_names():
        other_db = clint.get_database(other_db_name)
    else:
        print("没有【{0}】这个库".format(other_db_name))
        return
    print("==========当前库：{0}；目标库：{1}==========".format(src_db_name, other_db_name))
    for cname in src_db.list_collection_names():
        col = src_db.get_collection(cname)
        ocol = other_db.get_collection(cname)
        src_count = col.count_documents({})
        other_count = 0
        cursor = col.find(no_cursor_timeout=True)
        for ref in cursor:
            rec = {
                "摘要": ref["摘要"],
                "指标": ref["指标"],
                "物种": ref["物种"]
            }
            if ocol.find_one(rec):
                other_count += 1
        print("物种：{0}，当前库中数量：{1}，重复数量：{2}".format(cname, src_count, other_count))


def count_all_refs_col(db_name: str):
    # 函数说明：
    # 功能：统计一个库中所有集合中所有文档的数量
    # 参数：db_name：库名
    # 返回值：无
    if db_name not in clint.list_database_names():
        print("没有找到名为【{0}】的库".format(db_name))
        return
    db = clint.get_database(db_name)
    for cname in db.list_collection_names():
        col = db.get_collection(cname)
        ref_count = col.count_documents({})
        print("集合名：{0}，数量：{1}".format(cname, ref_count))


def count_gene_db_with_mea(db_name: str):
    # 函数说明：
    # 功能：统计一个库中各个物种里各个指标记录的数量
    # 参数：db_name：库名
    # 返回值：无
    if db_name not in clint.list_database_names():
        print("没有【{0}】这个库".format(db_name))
        return
    db = clint.get_database(db_name)
    for cname in db.list_collection_names():
        col = db.get_collection(cname)
        mea_count_dic = dict()
        cursor = col.find(no_cursor_timeout=True)
        for ref in cursor:
            if ref["物种"] in mea_count_dic.keys():
                mdic = mea_count_dic[ref["物种"]]
                if ref["指标"] in mdic.keys():
                    mdic[ref["指标"]] += 1
                else:
                    mdic[ref["指标"]] = 1
            else:
                mea_count_dic[ref["物种"]] = {ref["指标"]: 1}
        cursor.close()
        print(mea_count_dic)


clint = MongoClient()
db_list = ["新万方基因库2", "新知网基因库2", "新维普基因库2"]
for src in db_list:
    print("============开始统计库【{0}】=============".format(src))
    count_all_refs_col(src)
""" for dname in ["新万方摘要库", "新知网摘要库", "新维普摘要库"]:
    print("==============统计【{0}】===============".format(dname))
    count_all_refs_col(dname) """
""" for f in [1, 2, 3]:
    count_ref_with_org(f) """
""" db = clint.get_database("中国知网")
f = 3
tdb = "万方"
for cname in db.list_collection_names():
    print(count_repeat_num(cname, tdb, f)) """
""" flag = 1
clist = count_refs(db_flag=flag)
print(clist)
i = 0
att = '期刊名'
att_dict_list = []
for col in clist:
    cname = col['物种']
    att_names = get_att_name(col_name=cname, att=att)
    for n in att_names:
        att_dict_list.append(count_att_num(col_name=cname, att_name=n, att=att))
print(att_dict_list) """
# count_write_refs()
# get_att_and_write('期刊名')

""" {'当前库': '中国知网', '目标库': '万方', '集合': '羊', '总数': 15281, '重复数': 3503}
{'当前库': '中国知网', '目标库': '万方', '集合': '猫', '总数': 220, '重复数': 72}
{'当前库': '中国知网', '目标库': '万方', '集合': '鹅', '总数': 2133, '重复数': 883}
{'当前库': '中国知网', '目标库': '万方', '集合': '鸡', '总数': 4323, '重复数': 1618}
{'当前库': '中国知网', '目标库': '万方', '集合': '牛', '总数': 17280, '重复数': 4936}
{'当前库': '中国知网', '目标库': '万方', '集合': '狗', '总数': 1547, '重复数': 615}
{'当前库': '中国知网', '目标库': '万方', '集合': '鸭', '总数': 3260, '重复数': 1270}
{'当前库': '中国知网', '目标库': '万方', '集合': '鸽子', '总数': 80, '重复数': 23}
{'当前库': '中国知网', '目标库': '万方', '集合': '猪', '总数': 8859, '重复数': 3197}
{'当前库': '中国知网', '目标库': '万方', '集合': '马', '总数': 3454, '重复数': 1017} """

""" {'当前库': '中国知网', '目标库': '维普', '集合': '羊', '总数': 15281, '重复数': 3057}
{'当前库': '中国知网', '目标库': '维普', '集合': '猫', '总数': 220, '重复数': 30}
{'当前库': '中国知网', '目标库': '维普', '集合': '鹅', '总数': 2133, '重复数': 917}
{'当前库': '中国知网', '目标库': '维普', '集合': '鸡', '总数': 4323, '重复数': 1705}
{'当前库': '中国知网', '目标库': '维普', '集合': '牛', '总数': 17280, '重复数': 4755}
{'当前库': '中国知网', '目标库': '维普', '集合': '狗', '总数': 1547, '重复数': 519}
{'当前库': '中国知网', '目标库': '维普', '集合': '鸭', '总数': 3260, '重复数': 1242}
{'当前库': '中国知网', '目标库': '维普', '集合': '鸽子', '总数': 80, '重复数': 29}
{'当前库': '中国知网', '目标库': '维普', '集合': '猪', '总数': 8859, '重复数': 2919}
{'当前库': '中国知网', '目标库': '维普', '集合': '马', '总数': 3454, '重复数': 420} """

""" {'当前库': '维普', '目标库': '万方', '集合': '羊', '总数': 8464, '重复数': 2174}
{'当前库': '维普', '目标库': '万方', '集合': '猫', '总数': 56, '重复数': 31}
{'当前库': '维普', '目标库': '万方', '集合': '鹅', '总数': 1689, '重复数': 896}
{'当前库': '维普', '目标库': '万方', '集合': '鸡', '总数': 11829, '重复数': 4231}
{'当前库': '维普', '目标库': '万方', '集合': '牛', '总数': 20191, '重复数': 5559}
{'当前库': '维普', '目标库': '万方', '集合': '狗', '总数': 991, '重复数': 526}
{'当前库': '维普', '目标库': '万方', '集合': '鸭', '总数': 2002, '重复数': 1052}
{'当前库': '维普', '目标库': '万方', '集合': '鸽子', '总数': 64, '重复数': 29}
{'当前库': '维普', '目标库': '万方', '集合': '猪', '总数': 6787, '重复数': 2875}
{'当前库': '维普', '目标库': '万方', '集合': '马', '总数': 814, '重复数': 372} """
