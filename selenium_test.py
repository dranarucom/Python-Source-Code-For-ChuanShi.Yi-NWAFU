from selenium import webdriver
from selenium.common import exceptions
import time
import csv

# 从视频网站哔哩哔哩获取所有美国地区的动画番剧列表，写入到一个csv文件中

driver = webdriver.Chrome()
driver.implicitly_wait(3)
bangumi_dict_list = []
# 链接中"area=3"表示地区为美国，其它筛选条件为默认值
link = "https://www.bilibili.com/anime/index/#season_version=-1&area=3&is_finish=-1&copyright=-1&season_status=-1&season_month=-1&year=-1&style_id=-1&order=3&st=1&sort=0&page=1"
driver.get(link)
# 由于番剧分几页显示，所以用循环获取番剧信息，每一次循环获取一页的所有番剧
while True:
    # 先获取包含所有番剧信息的网页元素，再从中获取番剧元素列表，再从列表中取出每一个元素，获取所需的信息，保存为一个字典
    bangumi_list_content = driver.find_element_by_xpath(
        "//*[@id='app']/div[2]/div[1]/ul[2]")
    bangumi_list = bangumi_list_content.find_elements_by_tag_name('li')
    for bangumi in bangumi_list:
        publish_info = bangumi.find_element_by_class_name('pub-info').text
        follower_count = bangumi.find_element_by_class_name('shadow').text
        ban = bangumi.find_element_by_class_name('bangumi-title')
        title = ban.text
        link = ban.get_attribute('href')
        bangumi_dict = {
            '番剧名': title,
            '更新状态': publish_info,
            '追番人数': follower_count,
            '番剧链接': link
        }
        bangumi_dict_list.append(bangumi_dict)
    # 获取了这一页的所有番剧后，试着获取“下一页”的按钮，如果获取成功就模拟浏览器点击按钮以加载下一页内容，否则跳出while循环
    try:
        next_page_button = driver.find_element_by_xpath('//*[@id="app"]/div[2]/div[1]/div/a[@class="p next-page"]')
    except exceptions.NoSuchElementException:
        break
    else:
        # 加载了下一页内容后，要用time.sleep()等内容加载完成，否则容易出现StaleElementReferenceException错误
        next_page_button.click()
        time.sleep(1.5)
# 将获取到的信息写入csv文件
with open('美国番剧列表.csv', 'w', encoding='utf-8') as file:
    writer = csv.DictWriter(file, fieldnames={
        '番剧名',
        '更新状态',
        '追番人数',
        '番剧链接'
    })
    writer.writeheader()
    writer.writerows(bangumi_dict_list)
