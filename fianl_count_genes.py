from pymongo import MongoClient


def count_db_to_db(source_db: str, target_db: str):
    # 函数说明：
    # 功能：统计两个库中有多少相同的内容
    # 参数：source_db：做为对比基础的库；target_db：用于对比的库
    # 返回值：无
    if source_db not in clint.list_database_names():
        print("没有【{0}】这个库".format(source_db))
        return
    if target_db not in clint.list_database_names():
        print("没有【{0}】这个库".format(target_db))
        return

    if source_db == target_db:
        return

    sdb = clint.get_database(source_db)
    tdb = clint.get_database(target_db)

    print("库【{0}】和库【{1}】对比开始".format(source_db, target_db))
    for cname in sdb.list_collection_names():
        if cname not in tdb.list_collection_names():
            print("库【{0}】中没有【{1}】这个集合".format(target_db, cname))
            continue

        scol = sdb.get_collection(cname)
        tcol = tdb.get_collection(cname)
        res_dic = {"物种": cname, "基库总数": scol.count_documents({}), "对象库总数": tcol.count_documents({}), "重复数": 0}
        cursour = scol.find(no_cursor_timeout=True)
        for ref in cursour:
            rec = {
                "物种": ref["物种"],
                "基因": ref["基因"],
                "指标": ref["指标"]
            }
            if tcol.find_one(rec):
                res_dic["重复数"] += 1
        print(res_dic)
        cursour.close


clint = MongoClient()
source_db_list = ["新知网基因库", "新万方基因库", "新维普基因库"]
target_db_list = ["新知网基因库", "新万方基因库", "新维普基因库"]
for sdbname in source_db_list:
    for tdbname in target_db_list:
        count_db_to_db(sdbname, tdbname)
