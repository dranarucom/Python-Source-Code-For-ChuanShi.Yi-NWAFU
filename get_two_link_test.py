from selenium import webdriver
import time
from selenium.common import exceptions

driver = webdriver.Chrome()
driver.get('http://ziziyy.net/acg/')
time.sleep(2)
# 开启新标签页的JavaScript脚本，开启新标签页同时打开新网页
js = "window.open('http://www.acfun.cn/')"
# 执行JS脚本
driver.execute_script(js)
# 获取所有窗口的句柄，打印结果示例：['CDwindow-5FDB548FDDE42E72D2F0FB2ACCF74691', 'CDwindow-E639D61969C3D004E39062B5F39D6564', 'CDwindow-B888DD6228F077334F219FB37813F428']
n = driver.window_handles
print(n)
# 通过上面的JS脚本开启新标签页时，虽然浏览器显示在新标签页，但如果直接获取新标签页的元素会报错，需要通过driver.swtich_to.window()切换窗口句柄
# 后才能获取新标签页中的元素
try:
    driver.find_element_by_xpath(
        '//*[@id="pagelet_navigation"]/div/ul/li[2]/a')
except exceptions.NoSuchElementException:
    print('没有切换当前标签页')
else:
    print("切换了当前标签页")
# 因为此时还没有切换到新标签页，所以可以获取旧标签页的元素
try:
    t = driver.find_element_by_xpath('/html/body/div[1]/div/ul[1]/a')
except exceptions.NoSuchElementException:
    print("切换了当前标签页")
else:
    print(t.text)
    print("没有切换当前标签页")
time.sleep(1)
# 切换标签页（窗口）
driver.switch_to.window(n[0])
time.sleep(2)
driver.switch_to.window(n[-1])
try:
    driver.find_element_by_xpath(
        '//*[@id="pagelet_navigation"]/div/ul/li[2]/a')
except exceptions.NoSuchElementException:
    print("切换当前标签页不成功")
else:
    print("切换当前标签页成功")
time.sleep(5)
driver.close()
time.sleep(2)
