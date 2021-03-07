import pymongo

clint = pymongo.MongoClient()

database = clint['家畜种类']

collection = database['牛']

data = {'id': 123, 'name': 'kingname', 'age': 20, 'scalary': 99999}
collection.insert_one(data)
