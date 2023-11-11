"""
FilePath: /InstallmentScheduling/scripts/plot_tools/main.py
Description:  
Author: rthete
Date: 2023-09-09 13:53:19
LastEditTime: 2023-11-11 20:07:17
"""
from datetime import datetime

from Plotter import Plotter
from utils import generate_data_list, read_data_from_file


class PlotterLauncher:
    def error_time_bar():
        """
        30处理机，时间，柱状图
        """
        current_time = datetime.now()
        time_str = current_time.strftime("%Y-%m-%d_%H-%M-%S")

        # 任务量5000~15000，共11组数据
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        x_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        SIS_diff = read_data_from_file("../../output/exp_3_time/error_SIS_30_diff.csv")
        APMISRR_diff = read_data_from_file(
            "../../output/exp_3_time/error_APMISRR_30_diff.csv"
        )
        TolerMIS_diff = read_data_from_file(
            "../../output/exp_3_time/error_TolerMIS_30_diff.csv"
        )
        SIS_mean = read_data_from_file("../../output/exp_3_time/error_SIS_30_mean.csv")
        APMISRR_mean = read_data_from_file(
            "../../output/exp_3_time/error_APMISRR_30_mean.csv"
        )
        TolerMIS_mean = read_data_from_file(
            "../../output/exp_3_time/error_TolerMIS_30_mean.csv"
        )

        Plotter.plot_error_time_bar(
            x,
            x_data,
            APMISRR_mean[0],
            APMISRR_diff[0],
            TolerMIS_mean[0],
            TolerMIS_diff[0],
            SIS_mean[0],
            SIS_diff[0],
            f"../../output/error_plots/time_1_error_place_{time_str}",
        )
        Plotter.plot_error_time_bar(
            x,
            x_data,
            APMISRR_mean[1],
            APMISRR_diff[1],
            TolerMIS_mean[1],
            TolerMIS_diff[1],
            SIS_mean[1],
            SIS_diff[1],
            f"../../output/error_plots/time_2_error_place_{time_str}",
        )
        Plotter.plot_error_time_bar(
            x,
            x_data,
            APMISRR_mean[2],
            APMISRR_diff[2],
            TolerMIS_mean[2],
            TolerMIS_diff[2],
            SIS_mean[2],
            SIS_diff[2],
            f"../../output/error_plots/time_3_error_place_{time_str}",
        )
        Plotter.plot_error_time_bar(
            x,
            x_data,
            APMISRR_mean[3],
            APMISRR_diff[3],
            TolerMIS_mean[3],
            TolerMIS_diff[3],
            SIS_mean[3],
            SIS_diff[3],
            f"../../output/error_plots/time_4_error_place_{time_str}",
        )

    def error_time_box():
        pass

    def error_ur_box():
        """
        deprecated: 30处理机，利用率，箱线图
        """
        current_time = datetime.now()
        time_str = current_time.strftime("%Y-%m-%d_%H-%M-%S")

        for error_num in range(1, 5):
            # 按min和max生成数据
            data = generate_data_list(error_num)

            # 储存为csv
            # np.savetxt("./data/error_plot_ur/plot_ur_data_{}.csv".format(
            #     error_num), data, delimiter=",", fmt="%.6f")

            Plotter.plot_error_ur(
                data, f"../../output/error_plots/{error_num}_error_place_{time_str}"
            )
            
    def error_ur_box_combined():
        """
        30处理机，利用率，箱线图，四种情况组合
        """
        current_time = datetime.now()
        time_str = current_time.strftime("%Y-%m-%d_%H-%M-%S")
        
        data_path = "../../output/exp_3_all_exps/error_ur_s30_w5000.csv"
        fig_path = f"../../output/error_plots/ur_error_place_w5000_{time_str}"
        Plotter.plot_error_ur_box_combined(data_path, fig_path)


if __name__ == "__main__":
    # PlotterLauncher.error_time_bar()
    # PlotterLauncher.error_ur_box()
    # PlotterLauncher.error_time_box()
    PlotterLauncher.error_ur_box_combined()
    print("wtf")
