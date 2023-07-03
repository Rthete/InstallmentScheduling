/*
 * @FilePath: \InstallmentScheduling\include\method.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-05-15 15:54:33
 * @LastEditTime: 2023-06-23 16:38:51
 */
#ifndef _METHOD_H
#define _METHOD_H

#include "MISRR.h"
#include "SIS.h"
#include "PMIS.h"
#include "APMISRR.h"
#include "myAPMISRR.h"
#include "MISRRL.h"
#include "MISRRLL.h"
#include <tuple>

double run_SIS(int server_num=15, int workload=8000, double theta=0.3, string data_path="../data/15-servers-w-20/", vector<int> error_server={});
void run_PMIS();
double run_MISRR(int server_num=15, int m=8, int workload=8000, double theta=0.3, string data_path="../data/15-servers-w-20/", vector<int> error_server={});
void run_MISRR_check();
double run_APMISRR(double lambda);
double run_APMISRR_cost(double lambda, int m);
double run_myAPMISRR(int server_num=15, double lambda=0.5, int m=8, int workload=8000, double theta=0.3, string data_path="../data/15-servers-w-20/", vector<int> error_server={});
double run_MISRRL(double lambda, int m);
double run_MISRRLL(double lambda1, double lambda2, int m);

void test_MISRR_theta();
void test_MISRR_all();
void test_APMISRR_totalTime();
void test_APMISRR_alpha();
void test_APMISRR_beta();
void test_APMISRR_totalTime_cost();
void test_APMISRR_installment();
void test_myAPMISRR_installment();
void test_MISRRL_lambda();
void test_MISRRL_installment();
void test_MISRRL_all();
void test_MISRRLL_lambda2();
void test_MISRRLL_lambda1();
void test_MISRRLL_all();
void test_MISRRLL_2_lambda();

void compare_MISRR_and_MISRRL();

#endif //_METHOD_H
