from pymongo import MongoClient
import pymysql


def common_mongo_to_mysql(ab_db_name: str, gene_db_name: str, sql_table_name: str):
    if ab_db_name not in clint.list_database_names():
        print("没有【{0}】这个库".format(ab_db_name))
        return
    if gene_db_name not in clint.list_database_names():
        print("没有【{0}】这个库".format(gene_db_name))
        return

    ab_db = clint.get_database(ab_db_name)
    gene_db = clint.get_database(gene_db_name)

    sql_cursor = mysql_clint.cursor()
    sql_add_key = """id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
        TITLE VARCHAR(100),
        AB TEXT,
        ORG_NAME VARCHAR(5),
        MEA VARCHAR(10),
        GENE_NAME VARCHAR(25)
        """
    create_sql_table = "CREATE TABLE IF NOT EXISTS {0} ({1})".format(sql_table_name, sql_add_key)
    sql_cursor.execute(create_sql_table)
    # INSERT IGNORE INTO会忽略数据库中已经存在的数据，如果数据库没有数据，就插入新的数据，如果有数据的话就跳过这条数据。这样就可以保留数据库中已经存在数据，达到在间隙中插入数据的目的
    insert_sql = "INSERT IGNORE INTO {0}(TITLE, AB, ORG_NAME, MEA, GENE_NAME) VALUES(%s, %s, %s, %s, %s)".format(sql_table_name)

    for cname in ab_db.list_collection_names():
        ab_col = ab_db.get_collection(cname)
        gene_col = gene_db.get_collection(cname)
        gene_cursor = gene_col.find({})
        for ref in gene_cursor:
            gname = ref['基因']
            m = ref['指标']
            org = ref['物种']
            title = ref['标题']
            r = ab_col.find_one({'标题': title, '指标': m, '物种': org})
            if r:
                ab = r['摘要']
            else:
                ab = "NULL"

            try:
                # 执行插入语句
                sql_cursor.execute(insert_sql, [title, ab, org, m, gname])
                # 提交到数据库执行
                mysql_clint.commit()
            except pymysql.Error:
                # 发生错误时回滚
                mysql_clint.rollback()
                print((title, ab, org, m, gname))
        gene_cursor.close()
    sql_cursor.close()


clint = MongoClient()
sql_host = "localhost"
sql_port = 3306
sql_user = "root"
sql_pwd = "root"
sql_db_name = 'ycs_gene_db'
mysql_clint = pymysql.connect(host=sql_host, port=sql_port, user=sql_user, password=sql_pwd, db=sql_db_name)
ab_db_list = ["新知网摘要库", "新万方摘要库", "新维普摘要库"]
gene_db_list = ["新知网基因库", "新万方基因库", "新维普基因库"]
sql_table_list = ["知网", "万方", "维普"]
ignore_list = ["知网"]
for i in range(0, 3):
    if sql_table_list[i] in ignore_list:
        continue
    common_mongo_to_mysql(ab_db_list[i], gene_db_list[i], sql_table_list[i])
