xpath_list = []
for i in range(20):
    xpath = '//*[@id="list_div_aa_1"]/div/'
    for j in range(i):
        xpath += 'div[3]/'
    xpath += 'div[2]'
    xpath_list.append(xpath)
print(xpath_list)
