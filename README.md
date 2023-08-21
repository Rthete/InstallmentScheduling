<!--
 * @FilePath: \InstallmentScheduling\README.md
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 15:04:04
 * @LastEditTime: 2023-08-19 17:37:46
-->
# README

Implementation of several multi-installment scheduling models, including PMIS, FaDMIS, APMIS-RR.

## Usage

```
InstallmentScheduling
├─ data
│  ├─ 15-servers-w-20
│  ├─ 15-servers-w-50
│  ├─ exp1-30-servers
│  ├─ APMISRR-homo-no_overhead
│  ├─ MISRR-hete-no_overhead
│  └─ MISRR-homo-no_overhead
├─ dataGenerator.py               // randomly generate data using numpy
├─ exp                            // experiments
│  └─ exp_1.cpp
│  ├─ exp_2.cpp
│  └─ exp_3.cpp         
├─ include                        // header files
├─ main.cpp
├─ models                         // implementation of models
└─ README.md

```