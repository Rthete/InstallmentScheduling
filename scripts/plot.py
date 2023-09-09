from cmath import pi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mtick
from mpl_toolkits.axisartist.axislines import Subplot
from sympy import Line2D

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


# 画图
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
    plt.bar(x, SIS_avg, yerr=SIS_diff, width=width, label=r'$\mathrm{SIS}$', color='#A5A5A5', error_kw=err_attr)
    plt.bar(x + width, APMISRR_avg, yerr=APMISRR_diff, width=width, label=r'$\mathrm{APMISRR}$', color='#5B9BD5', error_kw=err_attr)
    plt.bar(x + 2 * width, TolerMIS_avg, yerr=TolerMIS_diff, width=width, label=r'$\mathrm{TolerMIS}$', color='#ED7D31', error_kw=err_attr)

    plt.legend(prop={"family": "Times New Roman", 'size': 9}, bbox_to_anchor=(0.28, 0.98))  # 显示上面的label
    plt.yticks(fontproperties='Times New Roman', size=6)
    plt.xticks(fontproperties='Times New Roman', size=6)

    plt.xlabel(r'$\mathit{W}$', fontsize=9)
    plt.ylabel(r'$\mathrm{Time}$', fontsize=9)
    # plt.yscale("log")
    plt.xticks(x_data)

    # plt.ticklabel_format(style='scientific', scilimits=(0, 2), axis='y')
    # plt.show(block = True)

    plt.savefig(name, dpi=1000)

if __name__ == "__main__":
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    x_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    SIS_diff = read_data_from_file('./output/exp_3/error_SIS_30_diff.csv')
    APMISRR_diff = read_data_from_file('./output/exp_3/error_APMISRR_30_diff.csv')
    TolerMIS_diff = read_data_from_file('./output/exp_3/error_TolerMIS_30_diff.csv')
    SIS_mean = read_data_from_file('./output/exp_3/error_SIS_30_mean.csv')
    APMISRR_mean = read_data_from_file('./output/exp_3/error_APMISRR_30_mean.csv')
    TolerMIS_mean = read_data_from_file('./output/exp_3/error_TolerMIS_30_mean.csv')
    plot(x, x_data, APMISRR_mean[0], APMISRR_diff[0], TolerMIS_mean[0], TolerMIS_diff[0], SIS_mean[0], SIS_diff[0], "1_error_place")
    plot(x, x_data, APMISRR_mean[1], APMISRR_diff[1], TolerMIS_mean[1], TolerMIS_diff[1], SIS_mean[1], SIS_diff[1], "2_error_place")
    plot(x, x_data, APMISRR_mean[2], APMISRR_diff[2], TolerMIS_mean[2], TolerMIS_diff[2], SIS_mean[2], SIS_diff[2], "3_error_place")
    plot(x, x_data, APMISRR_mean[3], APMISRR_diff[3], TolerMIS_mean[3], TolerMIS_diff[3], SIS_mean[3], SIS_diff[3], "4_error_place")