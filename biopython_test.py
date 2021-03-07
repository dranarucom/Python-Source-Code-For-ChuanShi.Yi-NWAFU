from Bio import Entrez
import pandas
import csv
from urllib import error
import time
from http import client
import json
from os import path
from os import environ
from pymongo import MongoClient


def get_genes(orga_name_chi: str, orga_name_eng: str):
    # 函数说明：
    # 功能：从NCBI的gene库获取一个物种的基因名称，产生一个csv文件
    # 参数：orga_name_chi：物种的中文名；orga_name_eng：物种的英文名（俗名或学名）
    # 返回值：无
    ret_max = 10000
    if orga_name_chi == '家牛':
        ret_start = 82689
    else:
        ret_start = 0
    file_name = '基因列表/' + orga_name_chi + '基因名称' + '.csv'
    search_term = orga_name_eng + '[Organism]'
    open_mode = 'a'
    with open(file_name, open_mode) as file:
        writer = csv.DictWriter(file,
                                fieldnames={'ID', 'Name', 'OtherAliases'})
        writer.writeheader()
        while True:
            try:
                handle = Entrez.esearch(db='gene',
                                        term=search_term,
                                        retstart=ret_start,
                                        retmax=ret_max)
            except TimeoutError:
                print("Search Time Out!")
                break
            record = Entrez.read(handle)
            id_list = record['IdList']
            if len(id_list) == 0:
                break
            gene_list = []
            for i in range(0, len(id_list)):
                try:
                    h = Entrez.esummary(db='gene', id=id_list[i])
                    r = Entrez.read(h)
                except TimeoutError:
                    print(str(i) + ":Time Out!")
                    break
                except error.HTTPError:
                    print("HTTPError: " + id_list[i])
                    break
                except error.URLError:
                    print("URLError: " + id_list[i])
                    break
                except client.RemoteDisconnected:
                    print("RemoteDisconnected: " + id_list[i])
                    break
                try:
                    name = r['DocumentSummarySet']['DocumentSummary'][0][
                        'Name']
                    other = r['DocumentSummarySet']['DocumentSummary'][0][
                        'OtherAliases']
                except IndexError:
                    print('IndexError: ')
                    print(r)
                    continue
                gene_list.append({
                    'ID': (ret_start + 1),
                    'Name': name,
                    'OtherAliases': other
                })
                ret_start += 1
                time.sleep(0.3)
            writer.writerows(gene_list)


def get_genes_json(orga_name_chi: str, orga_name_eng: str):
    # 函数说明：
    # 功能：从NCBI的gene库获取一个物种的基因名称，产生一个json文件
    # 参数：orga_name_chi：物种的中文名；orga_name_eng：物种的英文名（俗名或学名）
    # 返回值：无
    ret_max = 10000
    file_name = orga_name_chi + '基因名称' + '.json'
    search_term = orga_name_eng + '[Organism]'
    count = 0
    if not path.isfile(file_name):
        with open(file_name, "w+", encoding="utf-8") as nf:
            json.dump({"count": count}, nf)
    else:
        with open(file_name, "r", encoding="utf-8") as fr:
            count = json.load(fr)["count"]
    print("start! name: " + orga_name_chi)
    # 获取物种基因数量
    hcount = Entrez.esearch(
        db='gene',
        term=search_term,
    )
    rcount = Entrez.read(hcount)
    ret_max = rcount['Count']
    # 开始获取所有基因
    try:
        handle = Entrez.esearch(db='gene',
                                term=search_term,
                                retstart=count,
                                retmax=ret_max)
        record = Entrez.read(handle)
        id_list = record['IdList']
    except (TimeoutError, ConnectionAbortedError, error.URLError):
        print("ESearch Net Error!")
    else:
        gene_list = []
        for i in id_list:
            try:
                h = Entrez.esummary(db='gene', id=i)
                r = Entrez.read(h)
            except (TimeoutError, error.HTTPError, error.URLError,
                    client.RemoteDisconnected, ValueError,
                    ConnectionAbortedError, ConnectionError):
                print(str(count) + ":Get Summary Error!")
                break
            try:
                name = r['DocumentSummarySet']['DocumentSummary'][0]['Name']
            except IndexError:
                print(str(count) + ":Get Summary Index Error!")
                print(r)
                continue
            count += 1
            gene_list.append({"id": count, "name": name})
            if (count % 1000 == 0):
                print("Get gene count: " + str(count))
        with open(file_name, "r", encoding="utf-8") as fr:
            jr = json.load(fr)
            for g in gene_list:
                jr["count"] += 1
                key_name = "key" + str(int(g["id"] / ret_max))
                if key_name in jr.keys():
                    if g["name"] in jr[key_name]:
                        continue
                    jr[key_name].append(g["name"])
                    continue
                jr[key_name] = [g["name"]]

        with open(file_name, "w+", encoding="utf-8") as fw:
            json.dump(jr, fw)
            print("dump over, name: " + orga_name_chi + ", count: " +
                  str(jr["count"]))
    print("over! name: " + orga_name_chi)


def get_genes_db(orga_name_chi: str, orga_name_eng_list: list):
    # 函数说明：
    # 功能：从NCBI的gene库获取一个物种的基因名称，保存在本地数据库中
    # 参数：orga_name_chi：物种的中文名；orga_name_eng：物种的英文名（俗名或学名）
    # 返回值：无
    client = MongoClient()
    print("=============start! name: " + orga_name_chi)
    # 获取物种基因数量
    ret_max = 0
    for name in orga_name_eng_list:
        search_term = name + '[Organism]'
        gene_count = get_gene_count_in_ncbi(search_term)
        if gene_count == 0:
            continue
        else:
            ret_max = gene_count
            break
    print(orga_name_chi + "总数: " + str(ret_max))
    if ret_max == 0:
        return

    # 连接数据库
    db_name = "物种基因名称库"
    if db_name in client.list_database_names():
        db = client.get_database(db_name)
    else:
        db = client[db_name]
    col_name = orga_name_chi + "基因名称库"
    if col_name in db.list_collection_names():
        col = db.get_collection(col_name)
    else:
        col = db[col_name]
    count = col.count_documents({})

    # 开始获取所有基因
    try:
        handle = Entrez.esearch(db='gene',
                                term=search_term,
                                retstart=count,
                                retmax=ret_max)
        record = Entrez.read(handle)
        id_list = record['IdList']
    except Exception:
        print("ESearch Net Error!")
    else:
        for gid in id_list:
            if col.find_one({'基因id': gid}):
                continue
            try:
                h = Entrez.esummary(db='gene', id=gid)
                r = Entrez.read(h)
            except Exception:
                print(str(count) + ":Get Summary Error!")
                break
            # time.sleep(0.3)
            try:
                name = r['DocumentSummarySet']['DocumentSummary'][0]['Name']
            except IndexError:
                name = 'null'
            count += 1
            rec = {'基因id': gid, '基因名': name}
            col.insert_one(rec)
            if (count % 1000 == 0):
                print("Get gene count: " + str(count) + "，orga: " +
                      orga_name_chi)


def get_gene_count_in_ncbi(search_term: str):
    # 函数说明：
    # 功能：从NCBI获取一个物种的基因数量
    # 参数：search_term：用于在ncbi中进行检索的语句
    # 返回值：基因数量
    try:
        count_handle = Entrez.esearch(
            db='gene',
            term=search_term,
        )
    except Exception:
        print("get count error")
        return 0
    count_read = Entrez.read(count_handle)
    count = count_read['Count']
    return count


Entrez.email = '1327133922@qq.com'
htt_set = 0  # 1.使用代理；其它：不使用代理
if htt_set == 1:
    # 设置代理ip
    environ["http_proxy"] = "http://61.135.169.121:80"
# https://ip.jiangxianli.com/?page=1&country=%E4%B8%AD%E5%9B%BD
orga_excel = pandas.read_excel('家畜名称.xlsx')
chi_name_list = orga_excel.iloc[0:0, ].columns
ignore_list = ["水牛", "奶牛", "疣鼻栖鸭", "乌鸡"]  # 奶牛与家牛有很多重复的内容，暂时先跳过，以后如果有机会再继续
result_file_mode = 2  # 1.csv，2.数据库, 3.json
for n in range(0, len(chi_name_list)):
    chi_name = chi_name_list[n]
    if chi_name in ignore_list:
        continue
    eng_name_list = []
    eng_names = orga_excel.iloc[:, n]
    for name in eng_names:
        if type(name) is not str:
            continue
        eng_name_list.append(name)
    if result_file_mode == 1:
        get_genes(chi_name, eng_name_list[0])
    elif result_file_mode == 2:
        get_genes_db(chi_name, eng_name_list)
    else:
        get_genes_json(chi_name, eng_name_list[0])
    print(chi_name + '已完成！========================')

# 家牛: 57904
# 奶牛：10747
# 山羊：13103
# 牦牛：20323
# 绵羊：7016
# 家兔：634
# 猫：1184
# 狗：5069
# 火鸡：686
# 家鸡：322
# 鸽子：140
# 乌鸡总数: 0
