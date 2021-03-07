def get_name_and_mearsure(search_pahse: str, flag: int):
    # 函数说明：
    # 功能：输入检索语句，提取其中的物种名和繁殖力指标信息
    # 参数：search_pahse：检索语句；flag：1.中国知网，2.万方，3.维普
    # 返回值：物种名和指标的字典
    import re
    if flag == 1:
        name = re.search("‘(.*?)’", search_pahse).group(1)
        mearsure = re.search(r"\*‘(.*?)’", search_pahse).group(1)
        return {"物种": name, "指标": mearsure}
    elif flag == 2:
        name = re.search('"(.*?)"', search_pahse).group(1)
        mearsure = re.search(r'\*"(.*?)"', search_pahse).group(1)
        return {"物种": name, "指标": mearsure}
    elif flag == 3:
        search_pahse = search_pahse + 'wp'
        name = re.search(r"R=(.*?) \*", search_pahse).group(1)
        mearsure = re.search(r" \* (.*?)wp", search_pahse).group(1)
        return {"物种": name, "指标": mearsure}
    else:
        print("参数错误，应为flag = 1/2/3")
        return None
