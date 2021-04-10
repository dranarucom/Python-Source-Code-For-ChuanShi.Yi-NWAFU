from pymongo import MongoClient
import pandas


def read_gene_list(db_name: str, org_name: str):
    file_base = "{0}基因列表.xlsx"
    file_name = file_base.format(org_name)
    print(file_name)
    if db_name in clint.list_database_names():
        db = clint.get_database(db_name)
    else:
        db = clint[db_name]
    if org_name in db.list_collection_names():
        col = db.get_collection(org_name)
    else:
        col = db[org_name]
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
        if pandas.isnull(gene_col[i]) or type(gene_col[i]) is not str:
            continue
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
        rec = {
            "gene_name": gene,
            "annotation": anno,
            "gene_id": gid
            }
        if col.find_one(rec):
            continue
        col.insert_one(rec)
        count += 1
    print("读取基因数量：{0}".format(count))


clint = MongoClient()
dname = "新物种基因库"
org_list = ["绵羊", "山羊", "奶牛", "水牛", "牦牛", "猪", "绿头鸭", "疣鼻栖鸭", "家鸡", "火鸡", "狗", "猫", "马"]
ig_org_list = []    # 跳过已完成的物种
for oname in org_list:
    if oname in ig_org_list:
        continue
    read_gene_list(dname, oname)
