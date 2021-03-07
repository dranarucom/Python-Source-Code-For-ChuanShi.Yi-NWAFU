import requests
import csv
import json
import time


def get_bangumi_list(link: str):
    # 函数说明
    # 功能：从给的链接中获取所有番剧的信息，保存为一个字典
    # 参数：link：包含番剧信息的链接
    # 返回值：获取到的番剧信息的列表

    # 从链接获取内容并转换为json格式
    html_str = requests.get(link).content.decode()
    json_data = json.loads(html_str)
    # 从json中获取所有番剧的列表
    bangumi_list = json_data['data']['list']
    bangumi_dict_list = []
    # 从番剧列表的每一个元素中提取一个番剧的信息，保存为一个字典，再将字典保存到字典列表中
    for bangumi in bangumi_list:
        bangumi_dict = {
            '番剧名': bangumi['title'],
            '追番人数': bangumi['order'],
            '播放状态': bangumi['index_show'],
            '收费情况': bangumi['badge'],
            '封面链接': bangumi['cover'],
            '番剧链接': bangumi['link']
        }
        bangumi_dict_list.append(bangumi_dict)
    # 返回番剧字典列表
    return bangumi_dict_list


# 番剧的链接，通过两个字符串的组合可以获取任意页的番剧信息（bangumi_link[0]+页数+bangumi_link[1]）
bangumi_link = [
    'https://api.bilibili.com/pgc/season/index/result?season_version=-1&area=-1&is_finish=-1&copyright=-1&season_status=-1&season_month=-1&year=-1&style_id=-1&order=3&st=1&sort=0&page=',
    '&season_type=1&pagesize=20&type=1']
bangumi_dict_list = []
# 获取前20页的番剧列表
for page_count in range(1, 21):
    link = bangumi_link[0] + str(page_count) + bangumi_link[1]
    bangumi_list = get_bangumi_list(link)
    bangumi_dict_list.append(bangumi_list)
    time.sleep(3)

# 将获取到的番剧列表保存为csv文件
with open('番剧列表.csv', 'w', encoding='utf-8') as file:
    writer = csv.DictWriter(file, fieldnames={
        '番剧名',
        '追番人数',
        '播放状态',
        '收费情况',
        '封面链接',
        '番剧链接'
    })
    writer.writeheader()
    # get_bangumi_list函数返回的是一个列表，所以bangumi_dict_list是列表的列表，不是字典的列表，所以不能直接用writerows函数
    # 把它写入文件，而要先获取每一个子列表，再用writerows函数把子列表写入文件，且不能用writerow函数写子列表，因为writer是一
    # 个DictWriter，writerow函数只能写字典，而writerows函数可以写字典列表
    for each_dict in bangumi_dict_list:
        writer.writerows(each_dict)
