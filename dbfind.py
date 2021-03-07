from pymongo import MongoClient


def del_repeat(db_name: str):
    clint = MongoClient()
    if db_name not in clint.list_database_names():
        return None
    db = clint.get_database(db_name)
    col_list = db.list_collection_names()
    for col in col_list:
        collection = db.get_collection(col)
        ref_list = collection.find()
        print(len(ref_list))


database_name = "中国知网"
del_repeat(database_name)
