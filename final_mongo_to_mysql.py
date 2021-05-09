from pymongo import MongoClient
import pymysql


def common_mongo_to_mysql(gene_db_name: str, sql_table_name: str):
    # 函数说明：
    # 功能：将mongoDB中的基因和来源的文献储存到MySQL中
    # 参数：ab_db_name：mongo中的摘要库；gene_db_name：mongo中的基因库；sql_table_name：储存数据的MySQL表
    # 返回值：无
    if gene_db_name not in clint.list_database_names():
        print("没有【{0}】这个库".format(gene_db_name))
        return

    gene_db = clint.get_database(gene_db_name)

    sql_cursor = mysql_clint.cursor()
    sql_add_key = """id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
        TITLE VARCHAR(100),
        AB TEXT,
        ORG_NAME VARCHAR(5),
        MEA VARCHAR(10),
        GENE_NAME VARCHAR(25),
        ANNOTATION VARCHAR(100)
        """
    create_sql_table = "CREATE TABLE IF NOT EXISTS {0} ({1})".format(sql_table_name, sql_add_key)
    sql_cursor.execute(create_sql_table)
    # INSERT IGNORE INTO会忽略数据库中已经存在的数据，如果数据库没有数据，就插入新的数据，如果有数据的话就跳过这条数据。这样就可以保留数据库中已经存在数据，达到在间隙中插入数据的目的
    insert_sql = "INSERT IGNORE INTO {0}(TITLE, AB, ORG_NAME, MEA, GENE_NAME, ANNOTATION) VALUES(%s, %s, %s, %s, %s, %s)".format(sql_table_name)

    for cname in gene_db.list_collection_names():
        gene_col = gene_db.get_collection(cname)
        gene_cursor = gene_col.find({})
        for ref in gene_cursor:
            gname = ref['基因']
            m = ref['指标']
            org = ref['物种']
            title = ref['标题']
            ab = ref["摘要"]
            anno = ref["注释"]

            try:
                # 执行插入语句
                sql_cursor.execute(insert_sql, [title, ab, org, m, gname, anno])
                # 提交到数据库执行
                mysql_clint.commit()
            except pymysql.Error:
                # 发生错误时回滚
                mysql_clint.rollback()
                print((title, ab, org, m, gname))
        gene_cursor.close()
    sql_cursor.close()


def add_column_to_mysql(source_db: str, source_key: str, source_insert_key: str, target_sql_table: str):
    # 函数说明：
    # 功能：从mongo的一个库中提取一系列值插入到MySQL的一个表中作为新的一列
    # 参数：source_db：mongo中的库；source_key：用于在mongo和MySQL之间标识同一个记录的值在mongo中对应的key；target_sql_table：要插入的MySQL表
    # 返回值：无
    if source_db in clint.list_database_names():
        sdb = clint.get_database(source_db)
    else:
        print("mongo中没有【{0}】这个库".format(source_db))
        return

    sql_cursor = mysql_clint.cursor()
    # 在GENE_NAME后面添加GENE_ANNOTATION这一列
    # sql_add_column = "ALTER TABLE {0} ADD COLUMN GENE_ANNOTATION VARCHAR(100) AFTER GENE_NAME".format(target_sql_table)
    # sql_cursor.execute(sql_add_column)
    sql_update = 'UPDATE {0} SET GENE_ANNOTATION = %s WHERE GENE_NAME = %s AND ORG_NAME = %s'.format(target_sql_table)

    for cname in sdb.list_collection_names():
        scol = sdb.get_collection(cname)
        cursor = scol.find({})
        for ref in cursor:
            gname = ref["gene_name"]
            anno = ref["annotation"]
            try:
                # 执行插入语句
                sql_cursor.execute(sql_update, (anno, gname, cname))
                # 提交到数据库执行
                mysql_clint.commit()
            except pymysql.Error:
                # 发生错误时回滚
                mysql_clint.rollback()
                print((anno, gname, cname))
        cursor.close()
    sql_cursor.close()


clint = MongoClient()
sql_host = "localhost"
sql_port = 3306
sql_user = "root"
sql_pwd = "root"
sql_db_name = 'ycs_gene_db'
mysql_clint = pymysql.connect(host=sql_host, port=sql_port, user=sql_user, password=sql_pwd, db=sql_db_name)
""" sql_table_list = ["知网", "万方", "维普"]
for i in range(0, 3):
    add_column_to_mysql("新物种基因库", "annotation", "gene_name", sql_table_list[i]) """

gene_db_list = ["新知网基因库2", "新万方基因库2", "新维普基因库2"]
sql_table_list = ["知网2", "万方2", "维普2"]
ignore_list = []
for i in range(0, 3):
    if sql_table_list[i] in ignore_list:
        continue
    common_mongo_to_mysql(gene_db_list[i], sql_table_list[i])
