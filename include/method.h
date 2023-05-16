/*
 * @FilePath: \InstallmentScheduling\include\method.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-05-15 15:54:33
 * @LastEditTime: 2023-05-16 16:54:50
 */
#ifndef FIRSTMODEL_METHOD_H
#define FIRSTMODEL_METHOD_H

#include "MISRR.h"
#include "PMIS.h"
#include "APMISRR.h"
#include "myAPMISRR.h"
#include "MISRRL.h"
#include <tuple>

void run_PMIS();
tuple<double, double> run_MISRR(double theta, int m);
double run_APMISRR(double lambda);
double run_APMISRR_cost(double lambda, int m);
double run_myAPMISRR(double lambda, int m);
double run_MISRRL(double lambda, int m);

void test_MISRR_theta();
void test_APMISRR_totalTime();
void test_APMISRR_alpha();
void test_APMISRR_beta();
void test_APMISRR_totalTime_cost();
void test_APMISRR_installment();
void test_myAPMISRR_installment();
void test_MISRRL_lambda();
void test_MISRRL_installment();

void compare_MISRR_and_MISRRL();

#endif //FIRSTMODEL_METHOD_H
