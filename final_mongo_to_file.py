from pymongo import MongoClient


def mongo_gene_to_txt(db_name: str, file_name: str):
    # 函数说明：
    # 功能：将mongo基因库中的内容写入一个txt文件中
    # 参数：db_name：库名；file_name：要写入的文件名
    # 返回值：无
    if db_name not in clint.list_database_names():
        print("没有【{0}】这个库".format(db_name))
        return
    sdb = clint.get_database(db_name)
    with open(file_name, "w") as fw:
        fw.write("物种\t指标\t基因\n")
        for cname in sdb.list_collection_names():
            col = sdb.get_collection(cname)
            cur = col.find(no_cursor_timeout=True)
            for ref in cur:
                org = ref["物种"]
                mea = ref["指标"]
                gene = ref["基因"]
                fw.write("{0}\t{1}\t{2}\n".format(org, mea, gene))


clint = MongoClient()
gdb_list = ["新万方基因库2", "新知网基因库2", "新维普基因库2"]
file_list = ["万方基因列表.txt", "知网基因列表.txt", "维普基因列表.txt"]
for i in range(0, 3):
    mongo_gene_to_txt(gdb_list[i], file_list[i])
