import json
from pymongo import MongoClient

with open("test_json.json", "w+", encoding="utf-8") as fw:
    j = {}
    k = "k" + str(0/10)
    for i in range(100):
        if i % 10 == 0:
            k = 'k' + str(i/10)
            j[k] = []
        j[k].append("name" + str(i))
    json.dump(j, fw)

""" with open("test_json.json", "r", encoding="utf-8") as fr:
    jr = json.load(fr)
    print(type(jr)) """
db_name = "物种基因名称库"
col_s = "奶牛基因名称库"
col_t = "家牛基因名称库"
clint = MongoClient()
if db_name in clint.list_database_names():
    db = clint.get_database(db_name)
    if col_s in db.list_collection_names() and col_t in db.list_collection_names():
        _s_col = db.get_collection(col_s)
        _t_col = db.get_collection(col_t)
        print("source count: " + str(_s_col.count_documents({})))  # 总数量
        count = 0
        for r in _s_col.find():
            if _t_col.find_one({"基因名": r["基因名"]}):
                count += 1
        print("dep count: " + str(count))  # 重复数量

""" source count: 10747
dep count: 10739 """
