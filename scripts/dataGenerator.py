'''
FilePath: \InstallmentScheduling\dataGenerator.py
Description:  
Author: rthete
Date: 2023-06-20 19:40:39
LastEditTime: 2023-09-01 00:46:52
'''

import numpy as np
np.random.seed(1)

# 生成处理机参数
def generate_30_servers():
    arr_g = np.random.rand(30) * 0.4 + 0.1
    arr_g.sort()
    arr_g = arr_g[::-1]

    arr_o = np.random.rand(30) * 5 + 0.5

    arr_s = np.random.rand(30) * 5 + 0.5

    arr_w = np.random.rand(30) * 15 + 35
    arr_w.sort()

    with open('data/exp1-30-servers/g.txt', mode = 'w') as data:
        for num in arr_g:
            data.write('{:.2f}\n'.format(num))

    with open('data/exp1-30-servers/o.txt', mode = 'w') as data:
        for num in arr_o:
            data.write('{:.2f}\n'.format(num))

    with open('data/exp1-30-servers/s.txt', mode = 'w') as data:
        for num in arr_s:
            data.write('{:.2f}\n'.format(num))

    with open('data/exp1-30-servers/w.txt', mode = 'w') as data:
        for num in arr_w:
            data.write('{:.2f}\n'.format(num))

# 生成随机整数
def generate_and_save_random_integers(n, error_num, min_value, max_value, output_file):
    # Generate n sets of random integers
    random_integers = np.random.randint(min_value, max_value + 1, size=(n, error_num))
    
    # Write the results to a .txt file
    with open(output_file, "w") as f:
        for integers in random_integers:
            integer_str = ", ".join(str(i) for i in integers)
            f.write(f"{integer_str}\n")
    
    print(f"Random integers saved to {output_file}")

# 生成故障处理机   
def generate_servers_error_place(server_num):
    # Specify the parameters
    num_sets = 30
    error_num = 3
    min_value = 1
    max_value = server_num
    output_file = (f"data/exp3-error-place/error-place-{max_value}-{error_num}.txt")

    # Call the function
    generate_and_save_random_integers(num_sets, error_num, min_value, max_value, output_file)

# 生成故障趟数
def generate_error_installment():
    # Specify the parameters
    num_sets = 30
    error_num = 1
    # 暂不考虑第一趟/最后一趟出错
    min_value = 2
    max_value = 23
    output_file = (f"data/exp3-error-place/error-installment.txt")

    # Call the function
    generate_and_save_random_integers(num_sets, error_num, min_value, max_value, output_file)
    
if __name__ == "__main__":
    generate_servers_error_place(30)
    # generate_error_installment()
    




