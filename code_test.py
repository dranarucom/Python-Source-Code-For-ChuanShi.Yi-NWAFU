def reverseStringII(my_str, size):
    first_str = my_str[size: len(my_str)]
    last_str = my_str[0: size]
    return (first_str + last_str)


clist = reverseStringII('abcdefg', 3)
print(clist)
