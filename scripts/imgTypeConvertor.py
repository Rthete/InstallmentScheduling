'''
FilePath: \InstallmentScheduling\scripts\imgTypeConvertor.py
Description:  
Author: rthete
Date: 2023-09-01 00:47:08
LastEditTime: 2023-09-09 13:09:37
'''
from PIL import Image
import cairosvg

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
def svg2pdf(input_svg, output_pdf):
    cairosvg.svg2pdf(file_obj=open(input_svg, 'rb'), write_to=output_pdf)

if __name__ == "__main__":
    input_path = "T vs. theta (n=30).svg"
    # png2eps(input_path, "output.eps")
    svg2pdf(input_path, "output.pdf")