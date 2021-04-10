from pymongo import MongoClient
import csv
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


def count_refs(db_flag=1):
    # 功能：统计某个库中文献的数量
    # 参数：db_flag：1.知网，2.万方，3.维普；col_name：集合名
    # 返回值：文献数量
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
    ref_num_dict_list = []
    for col in col_list:
        collection = db.get_collection(col)
        ref_num_dict_list.append({
            '物种': col,
            '文献数': collection.count_documents({})
        })
    return ref_num_dict_list


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


def count_write_refs():
    # 功能：统计各个库中各个集合的文献数量，并写入csv文件中
    # 参数：file_name：要写入的csv文件的路径
    # 返回值：没有
    for flag in range(1, 4):
        if flag == 1:
            file_name = '知网文献数量统计.csv'
        elif flag == 2:
            file_name = '万方文献数量统计.csv'
        else:
            file_name = '维普文献数量统计.csv'
        clist = count_refs(db_flag=flag)
        with open(file_name, 'w', encoding='utf-8') as file:
            writer = csv.DictWriter(file, fieldnames={'物种', '文献数'})
            writer.writeheader()
            writer.writerows(clist)


def get_att_and_write(att_name: str):
    # 功能：将所有论文分类并统计各类的论文数量，比如统计不同来源的论文数量
    # 参数：分类的基准，如果对文献来源，则为“期刊名”
    # 返回值：没有
    for flag in range(1, 4):
        clist = count_refs(db_flag=flag)
        for col in clist:
            cname = col['物种']
            if flag == 1:
                file_name = '知网-' + cname + '-' + att_name + '数量统计.csv'
            elif flag == 2:
                file_name = '万方-' + cname + '-' + att_name + '数量统计.csv'
            else:
                file_name = '维普-' + cname + '-' + att_name + '数量统计.csv'
            with open(file_name, 'w', encoding='utf-8') as file:
                writer = csv.DictWriter(file, fieldnames={att_name, '数量'})
                att_names = get_att_name(col_name=cname,
                                         att=att_name,
                                         db_flag=flag)
                for n in att_names:
                    writer.writerow(
                        count_att_num(col_name=cname,
                                      att_name=n,
                                      att=att_name,
                                      db_flag=flag))


def count_repeat_num(col_name: str, target_db_name: str, flag: int):
    # 函数说明：
    # 功能：计算当前库的一个集合中有多少篇文献在目标库中重复
    # 参数：col_name：库名；target_db_name：目标库名，flag：标志要统计的库：1.知网，2.万方，3.维普
    # 返回值：包含统计信息的一个字典
    if flag == 1:
        db_name = '中国知网'
    elif flag == 2:
        db_name = '万方'
    elif flag == 3:
        db_name = '维普'
    else:
        print("参数错误，flag应为1/2/3")
        return
    db = clint.get_database(db_name)
    if col_name not in db.list_collection_names():
        print("没有找到" + col_name + "这个库")
        return
    col = db.get_collection(col_name)
    col_count = col.count_documents({})
    repeat_num = 0
    target_db = clint.get_database(target_db_name)
    target_col = target_db.get_collection(col_name)
    from get_db_info import get_name_and_mearsure
    for r in col.find():
        fr = target_col.find_one({'标题': r['标题']})
        if fr:
            if target_db_name == "中国知网":
                f = 1
            elif target_db_name == "万方":
                f = 2
            elif target_db_name == "维普":
                f = 3
            else:
                print("没有" + target_db_name + "这个库")
                break
            if get_name_and_mearsure(r['搜索语句'], flag) == get_name_and_mearsure(
                    fr['搜索语句'], f):
                repeat_num += 1
    return {
        "当前库": db_name,
        "目标库": target_db_name,
        "集合": col_name,
        "总数": col_count,
        "重复数": repeat_num
    }


def count_gene_db(db_name: str):
    # 函数说明：
    # 功能：统计一个库中各个物种记录的数量
    # 参数：db_name：库名
    # 返回值：无
    if db_name not in clint.list_database_names():
        print("没有找到名为【{0}】的库".format(db_name))
        return
    db = clint.get_database(db_name)
    res_dic = dict()
    ignore_col_name = ["鹅", "鸽子"]
    ignore_org_name = ["乌鸡"]
    for cname in db.list_collection_names():
        if cname in ignore_col_name:
            continue
        col = db.get_collection(cname)
        cursor = col.find(no_cursor_timeout=True)
        for ref in cursor:
            org = ref["物种"]
            if org in ignore_org_name:
                continue
            if org in res_dic.keys():
                res_dic[org] += 1
            else:
                res_dic[org] = 1
        cursor.close()
    print(res_dic)


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
for dname in ["万方基因库", "知网基因库", "维普基因库"]:
    print("==============统计【{0}】===============".format(dname))
    count_gene_db_with_mea(dname)
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
