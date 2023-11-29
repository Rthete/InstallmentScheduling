"""
FilePath: /InstallmentScheduling/scripts/plot_tools/utils.py
Description:  
Author: rthete
Date: 2023-11-04 16:19:15
LastEditTime: 2023-11-11 20:37:18
"""
import numpy as np
np.random.seed(1)


def read_data_from_file(file_path):
    """
    读取数据并存入数组
    """
    with open(file_path, 'r') as file:
        # 初始化一个空二维数组来存储读取的数据
        data = []

        # 逐行读取文件内容
        for line in file:
            # 使用strip()方法去除每行末尾的换行符和空白字符
            line = line.strip()

            # 使用逗号分隔符将每行的内容拆分为多个数
            numbers_str = line.split(',')

            # 将每个字符串数值转换为整数，并添加到二维数组中
            row = []
            for num_str in numbers_str:
                try:
                    number = float(num_str)
                    row.append(number)
                except ValueError:
                    pass
                    # print(f"无法将 '{num_str}' 转换为整数。")

            data.append(row)

    # 返回二维数组
    return data


def generate_data_list(error_num):
    """
    从数据生成数组，返回error_num时的三个模型mock数据，格式为
    [
        [SIS数据]
        [APMISRR数据]
        [TolerMIS数据]
    ]
    """

    SIS_max = read_data_from_file('./output/exp_3_ur/error_SIS_30_max.csv')
    SIS_min = read_data_from_file('./output/exp_3_ur/error_SIS_30_min.csv')

    APMISRR_max = read_data_from_file(
        './output/exp_3_ur/error_APMISRR_30_max.csv')
    APMISRR_min = read_data_from_file(
        './output/exp_3_ur/error_APMISRR_30_min.csv')

    TolerMIS_max = read_data_from_file(
        './output/exp_3_ur/error_TolerMIS_30_max.csv')
    TolerMIS_min = read_data_from_file(
        './output/exp_3_ur/error_TolerMIS_30_min.csv')

    list_SIS = np.array([SIS_max[error_num - 1][0] * 100])
    list_SIS = np.append(list_SIS, np.random.uniform(SIS_min[error_num - 1][0] * 100,
                                                     SIS_max[error_num - 1][0] * 100, 88))
    list_SIS = np.append(list_SIS, SIS_min[error_num - 1][0] * 100)

    list_APMISRR = np.array([APMISRR_max[error_num - 1][0] * 100])
    list_APMISRR = np.append(list_APMISRR, np.random.uniform(APMISRR_min[error_num - 1][0] * 100,
                                                             APMISRR_max[error_num - 1][0] * 100, 88))
    list_APMISRR = np.append(list_APMISRR, APMISRR_min[error_num - 1][0] * 100)

    list_TolerMIS = np.array([TolerMIS_max[error_num - 1][0] * 100])
    list_TolerMIS = np.append(list_TolerMIS, np.random.uniform(TolerMIS_min[error_num - 1][0] * 100,
                                                               TolerMIS_max[error_num - 1][0] * 100, 88))
    list_TolerMIS = np.append(
        list_TolerMIS, TolerMIS_min[error_num - 1][0] * 100)

    list_data = []
    list_data.append(list_SIS)
    list_data.append(list_APMISRR)
    list_data.append(list_TolerMIS)

    return list_data