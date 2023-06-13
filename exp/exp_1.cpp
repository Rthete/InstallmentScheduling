/*
 * @FilePath: \InstallmentScheduling\exp\exp_1.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-06-12 18:59:30
 * @LastEditTime: 2023-06-13 19:38:47
 */
#include "exp_1.h"

void without_error_15_m() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_15_m.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int m = 3; m <= 40; m+=1) {
        fprintf(fpResult, "%.2f,", run_SIS(8000, 0.3));
        fprintf(fpResult, "%.2f,", run_myAPMISRR(0.5, m, 8000, 0.3));
        fprintf(fpResult, "%.2f\n", run_MISRR(m, 8000, 0.3));
    }
    fclose(fpResult);
}

void ur_without_error_15_m() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_15_m.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int m = 3; m <= 40; m+=1) {
        fprintf(fpResult, "%.4f,", run_SIS(8000, 0.3));
        fprintf(fpResult, "%.4f,", run_myAPMISRR(0.5, m, 8000, 0.3));
        fprintf(fpResult, "%.4f\n", run_MISRR(m, 8000, 0.3));
    }
    fclose(fpResult);
}

void without_error_15_W() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_15_W.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 25000; W+=1000) {
        fprintf(fpResult, "%.2f,", run_SIS(W, 0.3));
        fprintf(fpResult, "%.2f,", run_myAPMISRR(0.5, 24, W, 0.3));
        fprintf(fpResult, "%.2f\n", run_MISRR(24, W, 0.3));
    }
    fclose(fpResult);
}

void ur_without_error_15_W() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_15_W.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 25000; W+=1000) {
        fprintf(fpResult, "%.4f,", run_SIS(W, 0.3));
        fprintf(fpResult, "%.4f,", run_myAPMISRR(0.5, 24, W, 0.3));
        fprintf(fpResult, "%.4f\n", run_MISRR(24, W, 0.3));
    }
    fclose(fpResult);
}

void without_error_15_theta() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_15_theta.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(double theta = 0.1; theta <= 1; theta+=0.1) {
        fprintf(fpResult, "%.2f,", run_SIS(8000, theta));
        fprintf(fpResult, "%.2f,", run_myAPMISRR(0.5, 24, 8000, theta));
        fprintf(fpResult, "%.2f\n", run_MISRR(24, 8000, theta));
    }
    fclose(fpResult);
}

void ur_without_error_15_theta() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_15_theta.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(double theta = 0.1; theta <= 1; theta+=0.1) {
        fprintf(fpResult, "%.4f,", run_SIS(8000, theta));
        fprintf(fpResult, "%.4f,", run_myAPMISRR(0.5, 24, 8000, theta));
        fprintf(fpResult, "%.4f\n", run_MISRR(24, 8000, theta));
    }
    fclose(fpResult);
}