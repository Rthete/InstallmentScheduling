'''
FilePath: \InstallmentScheduling\dataGenerator.py
Description:  
Author: rthete
Date: 2023-06-20 19:40:39
LastEditTime: 2023-06-21 14:11:27
'''

import numpy as np

arr_g = np.random.rand(30) * 0.2 + 0.1
arr_g.sort()
arr_g = arr_g[::-1]

arr_o = np.random.rand(30) * 3 + 1

arr_s = np.random.rand(30) * 3 + 1

arr_w = np.random.rand(30) * 10 + 35
arr_w.sort()

with open('data/30-servers/g.txt', mode = 'w') as data:
    for num in arr_g:
        data.write('{:.2f}\n'.format(num))

with open('data/30-servers/o.txt', mode = 'w') as data:
    for num in arr_o:
        data.write('{:.2f}\n'.format(num))

with open('data/30-servers/s.txt', mode = 'w') as data:
    for num in arr_s:
        data.write('{:.2f}\n'.format(num))

with open('data/30-servers/w.txt', mode = 'w') as data:
    for num in arr_w:
        data.write('{:.2f}\n'.format(num))