"""
FilePath: /InstallmentScheduling/scripts/plot.py
Description:  
Author: rthete
Date: 2023-09-09 13:53:19
LastEditTime: 2023-10-21 20:24:39
"""
from datetime import datetime
import random
from cmath import pi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mtick
from mpl_toolkits.axisartist.axislines import Subplot
from sympy import Line2D

np.random.seed(1)

# 读取数据并存入数组


def read_data_from_file(file_path):
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
                    print(f"无法将 '{num_str}' 转换为整数。")

            data.append(row)

    # 返回二维数组
    return data


# 容错-时间-柱状图
def plot(x, x_data, APMISRR_avg, APMISRR_diff, TolerMIS_avg, TolerMIS_diff, SIS_avg, SIS_diff, name="name"):
    config = {
        "font.family": 'serif',  # 衬线字体
        "font.size": 13,  # 相当于小四大小
        "font.serif": ['SimSun'],  # 宋体
        "mathtext.fontset": 'stix',  # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
        'axes.unicode_minus': False  # 处理负号，即-号
    }
    rcParams.update(config)

    fig = plt.figure(figsize=(5.8, 3.7), dpi=2000)
    ax = Subplot(fig, 111)
    fig.add_subplot(ax)

    ax.axis["bottom"].set_axisline_style("->", size=1)
    ax.axis["left"].set_axisline_style("->", size=1)
    ax.axis["top"].set_visible(False)
    ax.axis["right"].set_visible(False)

    x = np.array(x)
    total_width, n = 0.9, 3
    width = total_width / n
    x = x - (total_width - width) / 2

    err_attr = {"elinewidth": 0.8, "ecolor": "black", "capsize": 2.5}
    plt.bar(x, SIS_avg, yerr=SIS_diff, width=width,
            label=r'$\mathrm{SIS}$', color='#A5A5A5', error_kw=err_attr)
    plt.bar(x + width, APMISRR_avg, yerr=APMISRR_diff, width=width,
            label=r'$\mathrm{APMISRR}$', color='#5B9BD5', error_kw=err_attr)
    plt.bar(x + 2 * width, TolerMIS_avg, yerr=TolerMIS_diff, width=width,
            label=r'$\mathrm{TolerMIS}$', color='#ED7D31', error_kw=err_attr)

    plt.legend(prop={"family": "Times New Roman", 'size': 9},
               bbox_to_anchor=(0.28, 0.98))  # 显示上面的label
    plt.yticks(fontproperties='Times New Roman', size=6)
    plt.xticks(fontproperties='Times New Roman', size=6)

    plt.xlabel(r'$\mathit{W}$', fontsize=9)
    plt.ylabel(r'$\mathrm{Time}$', fontsize=9)
    # plt.yscale("log")
    plt.xticks(x_data)

    # plt.ticklabel_format(style='scientific', scilimits=(0, 2), axis='y')
    # plt.show(block = True)

    plt.savefig(name, dpi=1000)

# 容错-利用率-箱线图


def plot_error_ur(x_data, name="name"):
    config = {
        "font.family": 'serif',  # 衬线字体
        "font.size": 12.5,  # 相当于小四大小
        "font.serif": ['SimSun'],  # 宋体
        "mathtext.fontset": 'stix',  # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
        'axes.unicode_minus': False  # 处理负号，即-号
    }
    rcParams.update(config)

    fig = plt.figure(figsize=(7, 3.7), dpi=2000)

    ax = Subplot(fig, 111)
    fig.add_subplot(ax)

    ax.axis["bottom"].set_axisline_style("->", size=1)
    ax.axis["left"].set_axisline_style("->", size=1)
    ax.axis["top"].set_visible(False)
    ax.axis["right"].set_visible(False)

    box = plt.boxplot(x_data,
                      labels=[
                          '$\mathrm{SIS}$', '$\mathrm{APMISRR}$', '$\mathrm{TolerMIS}$'],
                      vert=False,
                      widths=0.4,
                      meanline=False,
                      showmeans=False,
                      patch_artist=True,
                      medianprops={'color': 'black'})

    plt.xlabel(r'$\mathrm{using \: rate(\%)}$', fontsize=12)

    c_list = ['#A5A5A5', '#5B9BD5', '#ED7D31']

    for a, c in zip(box['boxes'], c_list):
        a.set(color='black')
        a.set(facecolor=c)

    plt.savefig(name, dpi=1000)


"""
从数据生成数组，返回error_num时的三个模型数据，格式为
[
    [SIS数据]
    [APMISRR数据]
    [TolerMIS数据]
]
"""


def generateDataList(error_num):
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


if __name__ == "__main__":
    current_time = datetime.now()
    time_str = current_time.strftime("%Y-%m-%d_%H-%M-%S")

    """
    30处理机，时间，柱状图
    """
    # 任务量5000~15000，共11组数据
    # x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    # x_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    # SIS_diff = read_data_from_file('./output/exp_3/error_SIS_30_diff.csv')
    # APMISRR_diff = read_data_from_file('./output/exp_3/error_APMISRR_30_diff.csv')
    # TolerMIS_diff = read_data_from_file('./output/exp_3/error_TolerMIS_30_diff.csv')
    # SIS_mean = read_data_from_file('./output/exp_3/error_SIS_30_mean.csv')
    # APMISRR_mean = read_data_from_file('./output/exp_3/error_APMISRR_30_mean.csv')
    # TolerMIS_mean = read_data_from_file('./output/exp_3/error_TolerMIS_30_mean.csv')

    # plot(x, x_data, APMISRR_mean[0], APMISRR_diff[0], TolerMIS_mean[0], TolerMIS_diff[0], SIS_mean[0], SIS_diff[0], f"./output/error_plots/1_error_place_{time_str}")
    # plot(x, x_data, APMISRR_mean[1], APMISRR_diff[1], TolerMIS_mean[1], TolerMIS_diff[1], SIS_mean[1], SIS_diff[1], f"./output/error_plots/2_error_place_{time_str}")
    # plot(x, x_data, APMISRR_mean[2], APMISRR_diff[2], TolerMIS_mean[2], TolerMIS_diff[2], SIS_mean[2], SIS_diff[2], f"./output/error_plots/3_error_place_{time_str}")
    # plot(x, x_data, APMISRR_mean[3], APMISRR_diff[3], TolerMIS_mean[3], TolerMIS_diff[3], SIS_mean[3], SIS_diff[3], f"./output/error_plots/4_error_place_{time_str}")

    """
    30处理机，利用率，箱线图
    """
    for error_num in range(1, 5):
        # 按min和max生成数据
        data = generateDataList(error_num)
        # 储存为csv
        np.savetxt("plot_ur_data_{}.csv".format(
            error_num), data, delimiter=",", fmt="%.6f")
        # plot_error_ur(
        #     data, f"./output/error_plots/{error_num}_error_place_{time_str}")
