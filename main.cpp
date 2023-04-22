/*
 * @FilePath: \InstallmentScheduling\main.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 16:19:26
 * @LastEditTime: 2023-04-22 16:32:17
 */

#include "include/MISRR.h"
#include "include/PMIS.h"

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
}

int main() {
    // runMISRR();
    runPMIS();
    return 0;
}