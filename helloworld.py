import requests
import re
import time
import os


def get_index_url(murl):
    # 函数说明：
    # 功能：从网页获取目标小说每个章节在线阅读的网址，通过爬取这些网址即可获得每个章节的内容
    # 参数：murl：记录了目录信息的网页网址
    # 返回值：返回以章节名为键、网址为值的字典

    # decode默认网站以utf-8格式编码，实际网站以GBK格式编码，所以将decode参数设为‘GBK'
    html_str = requests.get(murl).content.decode('GBK')
    url_index_dict = {}
    # 如果所需的字符串周围没有特殊的标识字符，可以先获取包含全部所需信息的大块字符，在
    # 从这大块字符中获取所需的信息，这样可以避免获取不必要的信息
    # 这里是获取了所有包含章节名和对应的网址的字符块，之后再从每一个字符块中获取章节名和网址
    url_contain = re.findall('<td class="ccss">(.*?)</td>', html_str, re.S)
    for contian in url_contain:
        # 网站将所有章节保存为一个表格，前面获取了表格中每一格的内容，其中包含没有章节的和章节只有图像而没有文字的，
        # 这里只获取有网址链接的章节和由文字内容组成的章节
        if not ('href' in contian):
            continue
        url = re.search('<a href="(.*?)">', contian, re.S).group(1)
        chapter_name = re.search('>(.*?)</a>', contian, re.S).group(1)
        if chapter_name == '插图':
            continue
        url_index_dict[chapter_name] = url
    return url_index_dict


def create_chapter_text(url, title):
    # 函数说明：
    # 功能：从网站获取章节内容，写入文本文件中
    # 参数：url：在线阅读章节的网址；title：要写入的文件的文件名
    # 返回值：无

    chapter = requests.get(url).content.decode('GBK')
    chapter_content = re.findall('&nbsp;&nbsp;&nbsp;&nbsp;(.*?)<br />', chapter)

    os.makedirs("恋姬无双小说", exist_ok=True)
    with open(os.path.join('恋姬无双小说', title + '.txt'), 'w', encoding='utf-8') as file:
        file.writelines(chapter_content)


main_url = 'http://book.suixw.com/modules/article/reader.php?aid=489'
url_dict = get_index_url(main_url)
for chapter in url_dict.keys():
    create_chapter_text(url_dict[chapter], chapter)
    time.sleep(3)
