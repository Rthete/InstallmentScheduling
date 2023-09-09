<!--
 * @FilePath: \InstallmentScheduling\README.md
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 15:04:04
 * @LastEditTime: 2023-09-09 13:34:07
-->
# README

Implementation of several multi-installment scheduling models, including PMIS, FaDMIS, APMIS-RR.

## Usage

```
InstallmentScheduling
├─ assets   // md文件所用图片
├─ data
│  ├─ 15-servers-w-20           // 初版15个处理机数据，计算能力20~30
│  ├─ 15-servers-w-50           // 初版15个处理机数据，计算能力50~54
│  ├─ APMISRR-homo-no_overhead  // 用于APMISRR模型验证的数据，同构系统无启动开销
│  ├─ exp1-15-servers           // 用于TolerMIS模型对比试验的15个处理机数据
│  ├─ exp1-30-servers           // 用于TolerMIS模型对比试验的30个处理机数据
│  ├─ exp3-error-place          // 用于TolerMIS模型对比试验的故障处理机编号、故障趟编号
│  ├─ MISRR-hete-no_overhead    // 用于MISRR模型验证的数据，异构系统无启动开销
│  └─ MISRR-homo-no_overhead    // 用于MISRR模型验证的数据，同构系统无启动开销
├─ exp
│  ├─ exp_1.cpp // 对比SIS, APMISRR, TolerMIS 
│  ├─ exp_2.cpp // TolerMIS单模型对比
│  └─ exp_3.cpp // 带容错的3个模型实验
├─ include
├─ main.cpp
├─ models   // 模型算法实现
├─ README.md
└─ scripts
   ├─ dataGenerator.py      // 随机生成处理机参数
   └─ imgTypeConvertor.py   // 转换图片格式

```