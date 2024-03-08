import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axisartist.axislines import Subplot
import numpy as np
import pandas as pd
import seaborn as sns


class Plotter:
    def plot_error_time_bar(
        x,
        x_data,
        APMISRR_avg,
        APMISRR_diff,
        TolerMIS_avg,
        TolerMIS_diff,
        SIS_avg,
        SIS_diff,
        result_fig,
    ):
        """
        容错-时间-柱状图
        """
        config = {
            "font.family": "serif",  # 衬线字体
            "font.size": 9,  # 相当于小四大小
            "font.serif": ["SimSun"],  # 宋体
            "mathtext.fontset": "stix",  # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
            "axes.unicode_minus": False,  # 处理负号，即-号
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
        plt.bar(
            x,
            SIS_avg,
            yerr=SIS_diff,
            width=width,
            label=r"$\mathrm{SIS}$",
            color="#A5A5A5",
            error_kw=err_attr,
        )
        plt.bar(
            x + width,
            APMISRR_avg,
            yerr=APMISRR_diff,
            width=width,
            label=r"$\mathrm{APMISRR}$",
            color="#5B9BD5",
            error_kw=err_attr,
        )
        plt.bar(
            x + 2 * width,
            TolerMIS_avg,
            yerr=TolerMIS_diff,
            width=width,
            label=r"$\mathrm{TolerMIS}$",
            color="#ED7D31",
            error_kw=err_attr,
        )

        plt.legend(
            prop={"family": "Times New Roman", "size": 9}, bbox_to_anchor=(0.28, 0.98)
        ) 
        plt.yticks(fontproperties="Times New Roman", size=6)
        plt.xticks(fontproperties="Times New Roman", size=6)

        plt.xlabel(r"$\mathit{W}$", fontsize=9)
        plt.ylabel(r"$\mathrm{Time}$", fontsize=9)
        # plt.yscale("log")
        x_label = []
        for i in range(5, 16):
            x_label.append(str(i * 1000))
        plt.xticks(x_data, x_label)

        # plt.ticklabel_format(style='scientific', scilimits=(0, 2), axis='y')
        # plt.show(block = True)

        plt.savefig(result_fig, dpi=1000)

    def plot_error_ur(x_data, result_fig):
        """
        deprecated: 容错-利用率-箱线图
        """
        config = {
            "font.family": "serif",  # 衬线字体
            "font.size": 12.5,  # 相当于小四大小
            "font.serif": ["SimSun"],  # 宋体
            "mathtext.fontset": "stix",  # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
            "axes.unicode_minus": False,  # 处理负号，即-号
        }
        rcParams.update(config)

        fig = plt.figure(figsize=(5, 8), dpi=2000)

        ax = Subplot(fig, 111)
        fig.add_subplot(ax)

        ax.axis["bottom"].set_axisline_style("->", size=1)
        ax.axis["left"].set_axisline_style("->", size=1)
        ax.axis["top"].set_visible(False)
        ax.axis["right"].set_visible(False)

        # ax.set_ylim(70, 100)
        ax.set_yticks(np.linspace(70, 100, 11, endpoint=True))

        box = plt.boxplot(
            x_data,
            labels=["$\mathrm{SIS}$", "$\mathrm{APMISRR}$", "$\mathrm{TolerMIS}$"],
            vert=True,
            widths=0.4,
            meanline=False,
            showmeans=False,
            patch_artist=True,
            medianprops={"color": "black"},
        )

        plt.ylabel(r"$\mathrm{using \: rate(\%)}$", fontsize=12)

        c_list = ["#A5A5A5", "#5B9BD5", "#ED7D31"]

        for a, c in zip(box["boxes"], c_list):
            a.set(color="black")
            a.set(facecolor=c)

        plt.savefig(result_fig, dpi=1000)

    def plot_error_time_box(
        x,
        x_data,
        APMISRR_avg,
        APMISRR_diff,
        TolerMIS_avg,
        TolerMIS_diff,
        SIS_avg,
        SIS_diff,
        result_fig,
    ):
        """
        deprecated: 容错-时间-箱线图
        """
        config = {
            "font.family": "serif",  # 衬线字体
            "font.size": 13,  # 相当于小四大小
            "font.serif": ["SimSun"],  # 宋体
            "mathtext.fontset": "stix",  # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
            "axes.unicode_minus": False,  # 处理负号，即-号
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
        plt.bar(
            x,
            SIS_avg,
            yerr=SIS_diff,
            width=width,
            label=r"$\mathrm{SIS}$",
            color="#A5A5A5",
            error_kw=err_attr,
        )
        plt.bar(
            x + width,
            APMISRR_avg,
            yerr=APMISRR_diff,
            width=width,
            label=r"$\mathrm{APMISRR}$",
            color="#5B9BD5",
            error_kw=err_attr,
        )
        plt.bar(
            x + 2 * width,
            TolerMIS_avg,
            yerr=TolerMIS_diff,
            width=width,
            label=r"$\mathrm{TolerMIS}$",
            color="#ED7D31",
            error_kw=err_attr,
        )

        plt.legend(
            prop={"family": "Times New Roman", "size": 9}, bbox_to_anchor=(0.28, 0.98)
        )  # 显示上面的label
        plt.yticks(fontproperties="Times New Roman", size=6)
        plt.xticks(fontproperties="Times New Roman", size=6)

        plt.xlabel(r"$\mathit{W(1e3)}$", fontsize=9)
        plt.ylabel(r"$\mathrm{Time}$", fontsize=9)
        # plt.yscale("log")
        plt.xticks(x_data)

        # plt.ticklabel_format(style='scientific', scilimits=(0, 2), axis='y')
        # plt.show(block = True)

        plt.savefig(result_fig, dpi=1000)


    def plot_error_ur_box_combined(data_path, result_fig):
        """
        容错-利用率-箱线图-组合
        
        seaborn: https://seaborn.pydata.org/generated/seaborn.boxplot.html#seaborn.boxplot
        """
        config = {
            "font.family": "serif",  # 衬线字体
            "font.size": 13,  # 相当于小四大小
            "font.serif": ["SimSun"],  # 宋体
            "mathtext.fontset": "stix",  # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
            "axes.unicode_minus": False,  # 处理负号，即-号
        }
        rcParams.update(config)

        fig = plt.figure(figsize=(6, 6), dpi=2000)

        df = pd.read_csv(data_path)
        color_list = ["#A5A5A5", "#5B9BD5", "#ED7D31"]
        sns.boxplot(
            x="error_num",
            y="value",
            data=df,
            hue="model",
            whis=20,
            palette=color_list,
            saturation=1,
            linecolor="black"
        )

        plt.legend(prop={"family": "Times New Roman", "size": 9}, loc="lower left")
        plt.yticks(fontproperties="Times New Roman", size=9)
        plt.xticks(fontproperties="Times New Roman", size=9)

        plt.xlabel(r"$\mathrm{Number \;\; of\;\; failed \;\; servers}$", fontsize=11)
        plt.ylabel(r"$\mathrm{Utilization \;\; ratio}$", fontsize=11)

        plt.savefig(result_fig, dpi=1000)
