from pymongo import MongoClient


def create_db():
    small_str = "你好世界"
    large_str = "你好世界排水"
    sdb = clint["测试库"]
    scol = sdb["测试集一"]
    ocol = sdb["测试集二"]
    scol.insert_one({"摘要": small_str})
    ocol.insert_one({"摘要": large_str})
    print("创建测试数据成功")


def test_find(src_col: str, other_col: str):
    db = clint.get_database("测试库")
    scol = db.get_collection(src_col)
    ocol = db.get_collection(other_col)
    ocount = 0
    for ref in scol.find():
        rec = {"摘要": ref["摘要"]}
        if ocol.find_one(rec):
            ocount += 1
    print("当前集合：{0}，目标集合：{1}，重复数量：{2}".format(src_col, other_col, ocount))


clint = MongoClient()
create_db()
test_find("测试集一", "测试集二")
test_find("测试集二", "测试集一")
