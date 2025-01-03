/*
 * @FilePath: \InstallmentScheduling\include\method.h
 * @Description:
 * @Author: rthete
 * @Date: 2023-05-15 15:54:33
 * @LastEditTime: 2024-01-05 18:06:05
 */
#ifndef _METHOD_H
#define _METHOD_H

#include <tuple>

#include "APMISRR.h"
#include "MISRR.h"
#include "MISRRL.h"
#include "MISRRLL.h"
#include "PMIS.h"
#include "SIS.h"
#include "myAPMISRR.h"
namespace ModelRunner {
extern bool output_using_rate; // 控制是否输出ur
double run_SIS(int server_num = 15, int workload = 8000, double theta = 0.3,
               string data_path = "../data/15-servers-w-20/",
               vector<int> error_server = {});
void run_PMIS();
double run_MISRR(int server_num = 15, int m = 8, int workload = 8000,
                 double theta = 0.3,
                 string data_path = "../data/15-servers-w-20/",
                 vector<int> error_server = {}, int error_installment = 10,
                 vector<Server> add_servers = {}, int add_installment = 10);
void run_MISRR_conflict(vector<double> &waiting_time, int server_num = 15,
                        int m = 8, int workload = 8000, double theta = 0.3,
                        string data_path = "../data/15-servers-w-20/",
                        vector<int> error_server = {},
                        int error_installment = 10);
void run_MISRR_check();
double run_APMISRR(double lambda);
double run_APMISRR_cost(double lambda, int m);
double run_myAPMISRR(int server_num = 15, double lambda = 0.5, int m = 8,
                     int workload = 8000, double theta = 0.3,
                     string data_path = "../data/15-servers-w-20/",
                     vector<int> error_server = {}, int error_installment = 10,
                     vector<Server> add_servers = {}, int add_installment = 10);
double run_MISRRL(double lambda, int m);
double run_MISRRLL(int server_num = 15, double lambda1 = 0.2,
                   double lambda2 = 0.3, int m = 24, int workload = 8000,
                   double theta = 0.3,
                   string data_path = "../data/15-servers-w-20/",
                   vector<int> error_server = {}, int error_installment = 10,
                   vector<Server> add_servers = {}, int add_installment = 10);
} // namespace ModelRunner

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
void test_MISRRLL_installment();
void test_MISRRLL_all();
void test_MISRRLL_2_lambda();

void compare_MISRR_and_MISRRL();

#endif //_METHOD_H
