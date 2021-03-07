from Bio import Entrez
import os


def get_org_gene_names(org_name_eng: str):
    search_term = org_name_eng + '[orgn]'
    handle = Entrez.esearch(db='gene',
                            term=search_term,
                            retstart=1000
                            )
    record = Entrez.read(handle)
    print(record['Count'])
    print(len(record['IdList']))


Entrez.email = "1327133922@qq.com"
os.environ["http_proxy"] = "112.30.164.18"
get_org_gene_names("goat")
