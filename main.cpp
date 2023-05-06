/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-05-06 21:06:43
 */

#include "include/MISRR.h"
#include "include/PMIS.h"
#include "include/APMISRR.h"

void runPMIS() {
    auto workload = 1000;   // total workload
    auto serverN = 15;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    auto m = 30;            // installment size
    
    auto fPMIS = fopen("../output/PMIS.txt", "w");
    fprintf(fPMIS, "----------m: %d----------\n", m);
    fprintf(fPMIS, "workload\tPMIS\n");

    PMIS pmis(serverN, theta);
    pmis.getDataFromFile();
    pmis.setW((double)workload);
    pmis.setM(m);
    pmis.initValue();
    pmis.getOptimalModel();

    fprintf(fPMIS, "%d\t\t%.2lf\n", workload, pmis.getOptimalTime());
    pmis.printResult();

    cout << "done" << endl;
}

void runMISRR() {
    auto workload = 1;   // total workload
    auto serverN = 15;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    auto m = 8;            // installment size
    
    auto fMISRR = fopen("../output/MISRR.txt", "w");
    fprintf(fMISRR, "----------m: %d----------\n", m);
    fprintf(fMISRR, "workload\tMISRR\n");

    MISRR misrr(serverN, theta, m);
    misrr.getDataFromFile();
    misrr.setW((double)workload);
    misrr.initValue();
    misrr.getOptimalModel();
    misrr.theLastInstallmentGap();

    fprintf(fMISRR, "%d\t\t%.2lf\n", workload, misrr.getOptimalTime());

    cout << "done" << endl;
}

/**
 * 运行不带启动开销的APMISRR算法
*/
double runAPMISRR(double lambda) {
    auto serverN = 14;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    auto m = 9;            // installment size

    cout << "**********************run APMISRR**********************" << endl;
    cout << serverN + 1 << " servers, m = " << m << ", theta = " << theta << ", lambda = " << lambda << endl;
    cout << "last installment load = " << lambda / m << endl;
    cout << "each internal installment load = " << (m - lambda) / (m * (m - 1)) << endl;

    APMISRR apmisrr(serverN, theta);
    apmisrr.getDataFromFile();
    apmisrr.setM((int)m);
    apmisrr.setLambda((double)lambda);
    apmisrr.initValue();
    apmisrr.isSchedulable();
    return apmisrr.getOptimalTime();
}

/**
 * 运行带启动开销的APMISRR算法
*/
double runAPMISRR_cost(double lambda) {
    auto serverN = 14;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    auto m = 8;            // installment size

    cout << "**********************run APMISRR**********************" << endl;
    cout << serverN << " servers, m = " << m << ", theta = " << theta << ", lambda = " << lambda << endl;
    cout << "last installment load = " << lambda / m << endl;
    cout << "each internal installment load = " << (m - lambda) / (m * (m - 1)) << endl;

    APMISRR apmisrr(serverN, theta);
    apmisrr.getDataFromFile();
    apmisrr.setM((int)m);
    apmisrr.setLambda((double)lambda);
    apmisrr.initValue_cost();
    if(apmisrr.isSchedulable_cost())
        return apmisrr.getOptimalTime_cost(); // test APMISRR total time
    // return apmisrr.getAlpha(); // test APMISRR alpha_1
    // return apmisrr.getBeta(); // test APMISRR beta_1
    return 0;
}

/**
 * 测试total time与\lambda的关系
*/
void test_APMISRR_totalTime() {
    double lambda[10];
    double time[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        time[i] = runAPMISRR(lambda[i]);
    }

    // 实验：总时间与\lambda线性正相关
    cout << "***********" << "(time[i] - time[i-1]) / 0.1" << "***********" << endl;
    for(int i = 1; i < 10; i++) {
        cout << (time[i] - time[i-1]) / 0.1 << endl;
    }
    cout << "Total makespan is linearly dependant on \\lambda." << endl;
}

/**
 * 测试alpha[1]与\lambda的关系
*/
void test_APMISRR_alpha() {
    double lambda[10];
    double time[10];
    double alpha1[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        alpha1[i] = runAPMISRR_cost(lambda[i]);
    }

    // 实验：alpha[1]与\lambda负线性相关
    cout << "***********" << "(alpha1[i] - alpha1[i-1]) / 0.1" << "***********" << endl;
    for(int i = 1; i < 10; i++) {
        cout << (alpha1[i] - alpha1[i-1]) / 0.1 << endl;
    }
    cout << "alpha_1 is linearly dependant on \\lambda." << endl;
}

/**
 * 测试beta[1]与\lambda的关系
*/
void test_APMISRR_beta() {
    double lambda[10];
    double time[10];
    double beta1[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        beta1[i] = runAPMISRR_cost(lambda[i]);
    }

    // 实验：beta[1]与\lambda负线性相关
    cout << "***********" << "(beta1[i] - beta1[i-1]) / 0.1" << "***********" << endl;
    for(int i = 1; i <= 10; i++) {
        cout << "lambda = " << lambda[i] << ", beta_1 = " << beta1[i] << endl;
        cout << (beta1[i] - beta1[i-1]) / 0.1 << endl;
    }
    cout << "beta_1 is linearly dependant on \\lambda." << endl;
}

/**
 * 测试total time与lambda的关系
*/
void test_APMISRR_totalTime_cost() {
    double lambda[10];
    double time[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        time[i] = runAPMISRR_cost(lambda[i]);
    }

    // // 实验：总时间与\lambda线性正相关
    // cout << "***********" << "(time[i] - time[i-1]) / 0.1" << "***********" << endl;
    // for(int i = 1; i < 10; i++) {
    //     cout << (time[i] - time[i-1]) / 0.1 << endl;
    // }
    // cout << "Total makespan is linearly dependant on \\lambda." << endl;
}

int main() {
    // runMISRR();
    // runPMIS();
    // runAPMISRR(0.6);
    // runAPMISRR_cost(0.6);
    // test_APMISRR_totalTime();
    // test_APMISRR_alpha();
    // test_APMISRR_beta();
    test_APMISRR_totalTime_cost();
    return 0;
}