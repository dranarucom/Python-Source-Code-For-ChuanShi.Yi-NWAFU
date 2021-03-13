from pymongo import MongoClient
import jieba


def search_ab(ab: str, org: str):
    vlist = jieba.lcut(ab)
    res_list = []
    g_db_name = "物种基因库"
    g_db = clint.get_database(g_db_name)
    if org not in g_db.list_collection_names():
        print("物种【{0}】没有对应的基因库".format(org))
        return
    col = g_db.get_collection(org)
    for v in vlist:
        if ('a' < v[0] < 'z') or ('A' < v[0] < 'Z'):
            if col.find_one({"基因": v}):
                res_list.append(v)
    return res_list


def search_all_ab(flag: int):
    if flag == 1:
        db_name = "知网摘要库"
        temp_db_name = "知网摘要暂存库"
        res_db_name = "知网基因库"
    elif flag == 2:
        db_name = "万方摘要库"
        temp_db_name = "万方摘要暂存库"
        res_db_name = "万方基因库"
    elif flag == 3:
        db_name = "维普摘要库"
        temp_db_name = "维普摘要暂存库"
        res_db_name = "维普基因库"
    else:
        print("flag的值应为1、2或3，现在flag的值为：{0}".format(flag))
        return
    db = clint.get_database(db_name)
    if temp_db_name in clint.list_database_names():
        temp_db = clint.get_database(temp_db_name)
    else:
        temp_db = clint[temp_db_name]
    if res_db_name in clint.list_database_names():
        res_db = clint.get_database(res_db_name)
    else:
        res_db = clint[res_db_name]

    ignore_col_name = ["鹅", "鸽子"]
    ignore_org_name = ["乌鸡"]
    for cname in db.list_collection_names():
        if cname in ignore_col_name:
            continue
        col = db.get_collection(cname)
        if cname in temp_db.list_collection_names():
            temp_col = temp_db.get_collection(cname)
        else:
            temp_col = temp_db[cname]
        if cname in res_db.list_collection_names():
            res_col = res_db.get_collection(cname)
        else:
            res_col = res_db[cname]

        cursor = col.find(no_cursor_timeout=True)
        for ref in cursor:
            if temp_col.find_one({"已完成id": ref["_id"]}):
                continue
            ab = ref["摘要"]
            org = ref["物种"]
            mea = ref["指标"]

            if org in ignore_org_name:
                continue
            if org in g_org_list:
                final_org = org
            else:
                if org in org_chg_dic.keys():
                    final_org = org_chg_dic[org]
                else:
                    print("没有【{0}】对应的基因集合".format(org))
                    continue

            glist = search_ab(ab, final_org)
            for g in glist:
                rec = {
                    "基因": g,
                    "指标": mea,
                    "物种": final_org
                }
                if res_col.find_one(rec):
                    continue
                res_col.insert_one(rec)
            temp_col.insert_one({"已完成id": ref["_id"]})
        print("检索摘要库【{0}】结束，共检索数量：{1}，得到结果数量：{2}".format(cname, temp_col.count_documents({}), res_col.count_documents({})))


clint = MongoClient()
g_org_list = ["兔子", "家鸡", "山羊", "水牛", "火鸡", "牛", "牦牛", "狗", "猪", "猫", "疣鼻栖鸭", "绵羊", "马", "鸭"]
org_chg_dic = {"鸡": "家鸡", "家牛": "牛", "奶牛": "牛", "犬": "狗", "羊": "绵羊", "绿头鸭": "鸭", "黄牛": "牛"}
f = 2
search_all_ab(f)
