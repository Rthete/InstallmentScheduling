/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-03-15 16:43:42
 */

#include "include/MISRR.h"

int main() {
    auto workload = 1000;   // total workload
    auto serverN = 15;      // number of servers
    auto theta = 0.3;       // Ratio of the output load size to input load size
    auto m = 30;            // installment size
    
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
    misrr.printResult();

    cout << "done" << endl;
    return 0;
}