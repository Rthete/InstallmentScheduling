/*
 * @FilePath: \InstallmentScheduling\models\method.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-05-15 15:51:07
 * @LastEditTime: 2023-05-24 13:25:45
 */

#include "method.h"
void run_PMIS() {
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

/**
 * 运行MISRR算法
*/
 tuple<double, double> run_MISRR(double theta, int m) {
    auto workload = 8000;   // total workload
    auto serverN = 15;      // number of servers
    // auto theta = 10;       // Ratio of the output load size to input load size
    // auto m = 8;            // installment size

    cout << "**********************run MISRR**********************" << endl;
    cout << serverN << " servers, m = " << m << ", theta = " << theta << endl;
    cout << "total load = " << 8000 << endl;

    MISRR misrr(serverN, theta, m);
    misrr.getDataFromFile();
    misrr.setW((double)workload);
    misrr.initValue();
    misrr.getOptimalModel();
    cout << "misrr.getUsingRate(): " << misrr.getUsingRate() << endl;
    cout << "misrr.getOptimalTime(): " << misrr.getOptimalTime() << endl;
    misrr.theLastInstallmentGap();

    // fprintf(fMISRR, "%d\t\t%.2lf\n", workload, misrr.getOptimalTime());
    // return make_tuple(misrr.getOptimalM(), misrr.getUsingRate());
    return make_tuple(misrr.getOptimalTime(), misrr.getUsingRate());
    // printf("%d\t\t%.2lf\n", workload, misrr.getOptimalTime());
}

double run_MISRR(int m) {
    auto workload = 8000;   // total workload
    auto serverN = 15;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    // auto m = 8;            // installment size

    cout << "**********************run MISRR**********************" << endl;
    cout << serverN << " servers, m = " << m << ", theta = " << theta << endl;
    cout << "total load = " << 8000 << endl;

    MISRR misrr(serverN, theta, m);
    misrr.getDataFromFile();
    misrr.setW((double)workload);
    misrr.initValue();
    misrr.getOptimalModel();
    cout << "misrr.getUsingRate(): " << misrr.getUsingRate() << endl;
    cout << "misrr.getOptimalTime(): " << misrr.getOptimalTime() << endl;
    misrr.theLastInstallmentGap();
    // return misrr.getUsingRate();
    return misrr.getOptimalTime();
}


/**
 * 运行不带启动开销的APMISRR算法
*/
double run_APMISRR(double lambda) {
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
double run_APMISRR_cost(double lambda, int m) {
    auto serverN = 14;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    // auto m = 8;            // installment size

    cout << "**********************run APMISRR**********************" << endl;
    cout << serverN + 1 << " servers, m = " << m << ", theta = " << theta << ", lambda = " << lambda << endl;
    cout << "total load = " << 8000 << endl;
    cout << "last installment load = " << lambda / m * 8000 << endl;
    cout << "each internal installment load = " << (m - lambda) * 8000 / (m * (m - 1)) << endl;

    APMISRR apmisrr(serverN, theta);
    apmisrr.getDataFromFile();
    apmisrr.setM((int)m);
    apmisrr.setLambda((double)lambda);
    apmisrr.initValue_cost();
    apmisrr.isSchedulable_cost();
    return apmisrr.getOptimalTime_cost(); // test APMISRR total time
    // return apmisrr.getAlpha(); // test APMISRR alpha_1
    // return apmisrr.getBeta(); // test APMISRR beta_1
    return 0;
}

/**
 * 运行带启动开销，非阻塞，去掉P0的APMISRR算法
*/
double run_myAPMISRR(double lambda, int m) {

    auto serverN = 15;
    auto theta = 0.3;

    cout << "**********************run myAPMISRR**********************" << endl;
    cout << serverN << " servers, m = " << m << ", theta = " << theta << ", lambda = " << lambda << endl;
    cout << "total load = " << 8000 << endl;
    cout << "last installment load = " << lambda / m * 8000 << endl;
    cout << "each internal installment load = " << (m - lambda) * 8000 / (m * (m - 1)) << endl;
    
    myAPMISRR myapmisrr(serverN, theta);
    myapmisrr.getDataFromFile();
    myapmisrr.setM((int)m);
    myapmisrr.setLambda((double)lambda);
    myapmisrr.initValue();
    // myapmisrr.getOptimalTime();
    int isSchedulable = -1;
    isSchedulable = myapmisrr.isSchedulable();
    if (isSchedulable != 1)
        return 0;
    // return myapmisrr.getUsingRate();
    return myapmisrr.getOptimalTime();
}

double run_MISRRL(double lambda, int m) {
    auto serverN = 15;
    auto theta = 0.3;
    auto workload = 8000;

    cout << "**********************run MISRRL**********************" << endl;
    cout << serverN << " servers, m = " << m << ", theta = " << theta << ", lambda = " << lambda << endl;
    cout << "total load = " << 8000 << endl;
    cout << "last installment load = " << lambda / m * 8000 << endl;
    cout << "each internal installment load = " << (m - lambda) * 8000 / (m * (m - 1)) << endl;
    
    MISRRL misrrl(serverN, theta, m);
    misrrl.getDataFromFile();
    misrrl.setW((double)workload);
    misrrl.setLambda((double)lambda);
    misrrl.initValue();
    misrrl.getOptimalModel();
    // misrrl.setM((int)m);
    // misrrl.isSchedulable();
    misrrl.theLastInstallmentGap(to_string(lambda));
    if(misrrl.isSchedulable == 0) {
        cout << "not schedulable" << endl;
        return 0;
    }
    cout << "misrrl.getUsingRate(): " << misrrl.getUsingRate() << endl;
    cout << "misrrl.getOptimalTime(): " << misrrl.getOptimalTime() << endl;
    // return misrrl.getUsingRate();
    return misrrl.getOptimalTime();
}

/**
 * MISRR: 测试theta对可行性的影响
*/
void test_MISRR_theta() {
    // FILE * fpResult;
    // fpResult = fopen("../output/test_MISRR_theta_calM.csv", "w");
    for(int i = 0; i < 20; i+=1) {
        run_MISRR(i, 8);
        // cout << "theta = " << i;
        // cout<< get<0>(run_MISRR(i, 8)) << endl;
        // cout << "\ttime = " << get<0>(run_MISRR(i, 0)) << "\tusing rate = " << get<1>(run_MISRR(i, 0)) << endl;
        // fprintf(fpResult, "%d,%lf,%lf\n", i, get<0>(run_MISRR(i, 0)), get<1>(run_MISRR(i, 0)));
    }
    // fclose(fpResult);
}

/**
 * MISRR:测试趟数m的影响
*/
void test_MISRR_all() {
    FILE * fpResult;
    fpResult = fopen("../output/MISRR/test_MISRR_time.csv", "w");
    for(int m = 3; m <= 40; m++) {
        fprintf(fpResult, "%lf\n", run_MISRR(m));
    }
    fclose(fpResult);
}

/**
 * 测试total time与\lambda的关系
*/
void test_APMISRR_totalTime() {
    double lambda[10];
    double time[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        time[i] = run_APMISRR(lambda[i]);
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
        alpha1[i] = run_APMISRR_cost(lambda[i], 8);
    }

    // 实验：alpha[1]与\lambda负线性相关
    cout << "***********" << "(alpha1[i] - alpha1[i-1]) / 0.1" << "***********" << endl;
    for(int i = 1; i < 10; i++) {
        cout << (alpha1[i] - alpha1[i-1]) / 0.1 << endl;
    }
    // cout << "alpha_1 is linearly dependant on \\lambda." << endl;
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
        beta1[i] = run_APMISRR_cost(lambda[i], 8);
    }

    // 实验：beta[1]与\lambda负线性相关
    cout << "***********" << "(beta1[i] - beta1[i-1]) / 0.1" << "***********" << endl;
    for(int i = 1; i <= 10; i++) {
        cout << "lambda = " << lambda[i] << ", beta_1 = " << beta1[i] << endl;
        cout << (beta1[i] - beta1[i-1]) / 0.1 << endl;
    }
    // cout << "beta_1 is linearly dependant on \\lambda." << endl;
}

/**
 * 测试total time与lambda的关系
*/
void test_APMISRR_totalTime_cost() {
    double lambda[10];
    double time[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        time[i] = run_APMISRR_cost(lambda[i], 8);
    }

    // 实验：总时间与\lambda线性正相关
    cout << "***********" << "(time[i] - time[i-1]) / 0.1" << "***********" << endl;
    for(int i = 1; i < 10; i++) {
        cout << (time[i] - time[i-1]) / 0.1 << endl;
    }
    // cout << "Total makespan is linearly dependant on \\lambda." << endl;
}

/**
 * APMISRR: 测试总时间与趟数m的关系
*/
void test_APMISRR_installment() {
    FILE * fpResult;
    fpResult = fopen("../output/test.csv", "w");
    for(int m = 3; m <= 20; m++) {
        fprintf(fpResult, "%lf\n", run_APMISRR_cost(0.2, m));
    }
    fclose(fpResult);
}

/**
 * myAPMISRR: 测试总时间与趟数m的关系
 * lambda < 0.2 || > 0.8不可行
*/
void test_myAPMISRR_installment() {
    FILE * fpResult;
    fpResult = fopen("../output/myAPMISRR/test_myAPMISRR_time.csv", "w");
    for(int m = 3; m <= 40; m++) {
        for(double lambda = 0.1; lambda < 1; lambda+=0.1) {
            fprintf(fpResult, "%lf,", run_myAPMISRR(lambda, m));
        }
        fprintf(fpResult, "\n");
    }
    fclose(fpResult);
}

/**
 * MISRRL: 测试总时间与lambda的关系
 * 结论：lambda越小，用时越短
*/
void test_MISRRL_lambda() {
    for(double lambda = 0; lambda < 1; lambda+=0.1) {
        run_MISRRL(lambda, 8);
        // run_MISRRL(lambda, 20);
    }
}

/**
 * MISRRL: 测试总时间与m的关系
 * 结论: 有最优趟数m
*/
void test_MISRRL_installment() {
    for(int m = 3; m < 40; m++) {
        run_MISRRL(0.2, m);
    }
}

void test_MISRRL_all() {
    FILE * fpResult;
    fpResult = fopen("../output/MISRRL/test_MISRRL_time.csv", "w");
    for(int m = 3; m <= 40; m++) {
        for(double lambda = 0.1; lambda < 1; lambda+=0.1) {
            fprintf(fpResult, "%lf,", run_MISRRL(lambda, m));
        }
        fprintf(fpResult, "\n");
    }
    fclose(fpResult);
}


void compare_MISRR_and_MISRRL() {
    for(int m = 3; m < 20; m++) {
        run_MISRR(0.3, m);
        run_MISRRL(0.2, m);
    }
}