'''
FilePath: \InstallmentScheduling\scripts\imgTypeConvertor.py
Description:  
Author: rthete
Date: 2023-09-01 00:47:08
LastEditTime: 2023-12-03 22:03:51
'''
import os
from PIL import Image
# import cairosvg

# png --> eps 
def png2eps(input_path, output_path):
    try:
        # 打开输入图片
        image = Image.open(input_path)
        # 创建一个新的EPS图像
        eps_image = Image.new("RGB", image.size, (255, 255, 255))
        eps_image.paste(image)
        # 保存EPS图像
        eps_image.save(output_path, "EPS")
        print("图片转换完成！")
    except Exception as e:
        print("转换失败：", str(e))
        
# svg --> pdf
# def svg2pdf(input_svg, output_pdf):
#     cairosvg.svg2pdf(file_obj=open(input_svg, 'rb'), write_to=output_pdf)

def process_images_in_folder(input_path):
    # 获取文件夹中所有文件
    files = os.listdir(input_path)

    # 遍历文件夹中的每个文件
    for file in files:
        # 拼接文件的完整路径
        file_path = os.path.join(input_path, file)

        # 检查文件是否是PNG图片
        if file.lower().endswith(".png"):
            # 调用processor函数处理图片
            output_path = os.path.join("D:\\OneDrive - stu.xidian.edu.cn\\workspace\\任务调度相关\\实验\\TolerMIS实验\\231129_eps格式图例\\eps\\故障处理实验\\using_rate", 
                                       os.path.splitext(os.path.basename(file_path))[0] + ".eps")
            png2eps(file_path, output_path)

if __name__ == "__main__":
    # input_path = "T vs. theta (n=30).svg"
    # png2eps(input_path, "output.eps")
    # svg2pdf(input_path, "output.pdf")
    folder_path = "D:\\OneDrive - stu.xidian.edu.cn\\workspace\\任务调度相关\\实验\\TolerMIS实验\\231129_eps格式图例\\png\\故障处理实验\\using_rate"
    process_images_in_folder(folder_path)