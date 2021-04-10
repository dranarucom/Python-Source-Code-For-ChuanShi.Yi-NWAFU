from pymongo import MongoClient
import jieba


def search_ab(ab: str, org: str):
    vlist = jieba.lcut(ab)
    res_list = []
    g_db_name = "新物种基因库"
    g_db = clint.get_database(g_db_name)
    if org not in g_db.list_collection_names():
        print("物种【{0}】没有对应的基因库".format(org))
        return
    col = g_db.get_collection(org)
    for v in vlist:
        if ('a' < v[0] < 'z') or ('A' < v[0] < 'Z'):
            if col.find_one({"gene_name": v}):
                res_list.append(v)
    return res_list


def search_all_ab(ab_source_db_name: str, res_db_name: str):
    # 函数说明：
    # 功能：从库中读取摘要，用jieba分词后获取摘要中包含的所有遗传位点
    # 参数：ab_source_db_name：摘要库；res_db_name：获取到的遗传位点储存的库
    # 返回值：无
    if ab_source_db_name in clint.list_database_names():
        db = clint.get_database(ab_source_db_name)
    else:
        print("没有【{0}】这个库".format(ab_source_db_name))
        return

    if res_db_name in clint.list_database_names():
        res_db = clint.get_database(res_db_name)
    else:
        res_db = clint[res_db_name]

    for cname in db.list_collection_names():
        col = db.get_collection(cname)
        if cname in g_org_list:
            final_org = cname
        else:
            if cname in org_chg_dic.keys():
                final_org = org_chg_dic[cname]
            else:
                print("没有【{0}】对应的基因集合".format(cname))
                return
        if final_org in res_db.list_collection_names():
            res_col = res_db.get_collection(final_org)
        else:
            res_col = res_db[final_org]

        count = 0
        cursor = col.find(no_cursor_timeout=True)
        for ref in cursor:
            ab = ref["摘要"]
            mea = ref["指标"]
            title = ref["标题"]
            count += 1
            glist = search_ab(ab, final_org)
            if type(glist) is list and glist:
                for g in glist:
                    rec = {
                        "基因": g,
                        "指标": mea,
                        "物种": final_org,
                        "标题": title
                    }
                    if res_col.find_one(rec):
                        continue
                    res_col.insert_one(rec)
        print("检索摘要库【{0}】结束，得到结果数量：{1}".format(cname, res_col.count_documents({})))
        print("共处理数量：{0}".format(count))
        cursor.close()


clint = MongoClient()
g_org_list = ["绵羊", "山羊", "奶牛", "水牛", "牦牛", "猪", "绿头鸭", "疣鼻栖鸭", "家鸡", "火鸡", "狗", "猫", "马"]
org_chg_dic = {"犬": "狗"}
source_db_list = ["新知网摘要库", "新万方摘要库", "新维普摘要库"]
res_db_list = ["新知网基因库", "新万方基因库", "新维普基因库"]
for i in range(0, 3):
    search_all_ab(source_db_list[i], res_db_list[i])
