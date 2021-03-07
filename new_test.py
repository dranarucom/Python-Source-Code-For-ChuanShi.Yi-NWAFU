# import json
from Bio import Entrez
from pymongo import MongoClient


def try_find_in_gene(term_to_find: str, org_name: str):
    # 函数说明：
    # 功能：在NCBI的gene库中检索，看某个物种中是否存在某个基因
    # 参数：term_to_find：待检索的基因名；org_name：物种名
    # 返回值：检索结果为存在，返回True；否则返回False
    search_term = str.format("({0}[Gene Name]) AND {1}[Organism]", term_to_find, org_name)
    try:
        h = Entrez.esearch(db="gene", term=search_term)
        hr = Entrez.read(h)
    except Exception:
        print("========search {0} failed=========".format(term_to_find))
        return False
    else:
        count = int(hr['Count'])
        if count > 0:
            return True
        else:
            return False


Entrez.email = "1327133922@qq.com"
client = MongoClient()
db_name = "维普基因库"
if db_name in client.list_database_names():
    db = client.get_database(db_name)
    col_name = "猫"
    if col_name in db.list_collection_names():
        col = db.get_collection(col_name)
        for ref in col.find():
            print("========id: {0}======".format(ref["_id"]))
""" with open("家畜名称.json", "r", encoding="utf-8") as name_file:
    names = json.load(name_file)
print(names)
test_terms = ["SPP1", "CSN3", "TNF", "LTF", "fullmetal", "count3", "dome"]
chi_name = "水牛"
eng_name = names[chi_name]
for term in test_terms:
    if try_find_in_gene(term, eng_name):
        print(term) """

# "SPP1", "CSN3", "TNF", "LTE", "fullmetal", "count3", "dome"

# (fullmeatal[Gene Name]) AND horse[Organism]
# 奶牛=家牛
# 乌鸡=0
# 鹅=0

""" ignore_names = ["乌鸡", "奶牛", "鹅"]
for name in names.keys():
    if name in ignore_names:
        continue
    search_term = names[name] + "[Organism]"
    h = Entrez.esearch(
        db='gene',
        term=search_term,
    )
    c = Entrez.read(h)
    count = c["Count"]
    print(search_term)
    print("=======name: " + name + ", count: " + str(count)) """
