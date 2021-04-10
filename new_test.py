# import json
from Bio import Entrez
from pymongo import MongoClient
import pandas
import numpy


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


def txt_to_csv_test(file_name: str, out_name: str):
    t = numpy.loadtxt(file_name)
    tdf = pandas.DataFrame(t)
    tdf.to_excel(out_name, index=False)


def read_gene_list_test(org: str):
    file_base = "{0}基因列表.xlsx"
    file_name = file_base.format(org)
    print(file_name)
    db_name = "物种基因英文键库"
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]
    if org in db.list_collection_names():
        col = db.get_collection(org)
    else:
        col = db[org]
    gene_colunm_name = "Official Symbol"
    chro_colunm_name = "Chromosome"
    annotation_colunm_name = "Annotation"
    mim_id_colunm_name = "MIM_ID"
    gene_id_colunm_name = "Gene_ID"
    gene_list_df = pandas.read_excel(file_name)
    indx_list = [gene_colunm_name, chro_colunm_name, annotation_colunm_name, mim_id_colunm_name, gene_id_colunm_name]
    gene_df = gene_list_df.loc[:, indx_list]
    gene_col = gene_df[gene_colunm_name]
    chro_col = gene_df[chro_colunm_name]
    anno_col = gene_df[annotation_colunm_name]
    mim_col = gene_df[mim_id_colunm_name]
    gene_id_col = gene_df[gene_id_colunm_name]
    count = 0
    print(len(gene_col))
    for i in range(len(gene_col)):
        gene = gene_col[i].replace(" ", "")
        if type(anno_col[i]) is str and "Chromosome" in anno_col[i]:
            anno = anno_col[i]
        elif type(mim_col[i]) is str and "Chromosome" in mim_col[i]:
            anno = mim_col[i]
        elif type(chro_col[i]) is str and "Chromosome" in chro_col[i]:
            anno = chro_col[i]
        else:
            anno = "Null"
        if pandas.isnull(gene_id_col[i]):
            gid = "Null"
        else:
            gid = int(gene_id_col[i])
        print("name: {0}, annotation: {1}, gene_id: {2}".format(gene, anno, gid))
        continue
        g = gene_col[i]
        rec = {"gene_name": g}
        if col.find_one(rec):
            continue
        col.insert_one(rec)
        count += 1
    print(count)


def test_search_term(term: str):
    db = clint.get_database("物种基因英文键库")
    col = db.get_collection("牛")
    search_term = " {0}".format(term)
    if col.find_one({"gene_name": search_term}):
        print("==========【{0}】在集合中".format(term))
    else:
        print("==========【{0}】不在集合中".format(term))


def test_object_id():
    db = clint.get_database("物种基因库")
    col = db.get_collection("牛")
    cursor = col.find(no_cursor_timeout=True)
    import time
    for ref in cursor:
        print(ref["_id"])
        time.sleep(20)


clint = MongoClient()
# test_object_id()
read_gene_list_test("马")
""" test_list = ["RHO", "NOS3", "LTF", "PRNP", "ALB", "IGF1", "rho", "nos3", "ltf", "igf"]
for t in test_list:
    test_search_term(t) """
""" clint = MongoClient()
org_names = ["马", "兔子", "猫", "狗", "火鸡", "家鸡", "疣鼻栖鸭", "鸭", "猪", "牦牛", "水牛", "牛", "山羊", "绵羊"]
for org in org_names:
    read_gene_list_test(org) """
""" fname = "0.txt"
oname = "horse.csv"
txt_to_csv_test(fname, oname) """

""" Entrez.email = "1327133922@qq.com"
client = MongoClient()
db_name = "维普基因库"
if db_name in client.list_database_names():
    db = client.get_database(db_name)
    col_name = "猫"
    if col_name in db.list_collection_names():
        col = db.get_collection(col_name)
        for ref in col.find():
            print("========id: {0}======".format(ref["_id"])) """
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
