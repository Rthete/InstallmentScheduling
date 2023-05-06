/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-05-04 23:24:22
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
    auto m = 5;            // installment size
    
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

double runAPMISRR(double lambda) {
    auto serverN = 14;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    auto m = 9;            // installment size

    cout << "**********************run APMISRR**********************" << endl;
    cout << serverN << " servers, m = " << m << ", theta = " << theta << ", lambda = " << lambda << endl;
    cout << "last installment load = " << lambda / m << endl;
    cout << "each internal installment load = " << (m - lambda) / (m * (m - 1)) << endl;
    
    // auto fAPMISRR = fopen("../output/APMISRR.txt", "w");
    // fprintf(fAPMISRR, "----------m: %d----------\n", m);
    // fprintf(fAPMISRR, "workload\tAPMISRR\n");

    APMISRR apmisrr(serverN, theta);
    apmisrr.getDataFromFile();
    apmisrr.setM((int)m);
    apmisrr.setLambda((double)lambda);
    apmisrr.initValue_cost();
    // apmisrr.initValue();
    // apmisrr.isSchedulable();
    // return apmisrr.getOptimalTime();
    return 0;
}

void testAPMISRR() {
    double lambda[10];
    double time[10];
    for(int i = 0; i < 11; i++) {
        lambda[i] = 0.1 * i;
        time[i] = runAPMISRR(lambda[i]);
    }

    // 实验：总时间与\lambda线性正相关
    // cout << "***********" << "(time[i] - time[i-1]) / 0.1" << "***********" << endl;
    // for(int i = 1; i < 10; i++) {
    //     cout << (time[i] - time[i-1]) / 0.1 << endl;
    // }
    // cout << "Total makespan is linearly dependant on \\lambda." << endl;
}

int main() {
    // runMISRR();
    // runPMIS();
    testAPMISRR();
    // runAPMISRR(0.6);
    return 0;
}