/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-03-20 19:54:19
 */

#include "include/MISRR.h"

int modelMISRR() {
    auto workload = 2000;   // total workload
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

    return 1;
}

int main() {

    if(modelMISRR()) {
        cout << "done" << endl;
    }

    return 0;
}