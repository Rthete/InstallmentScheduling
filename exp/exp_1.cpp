/*
 * @FilePath: \InstallmentScheduling\exp\exp_1.cpp
 * @Description: 对比SIS, APMISRR, TolerMIS (不含容错)
 * @Author: rthete
 * @Date: 2023-06-12 18:59:30
 * @LastEditTime: 2023-09-13 13:08:02
 */
#include "exp_1.h"

void without_error_15_m() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_15_m.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int m = 3; m <= 40; m+=1) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(15, 8000, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(15, 0.5, m, 8000, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.2f\n", ModelRunner::run_MISRR(15, m, 8000, 0.3, "../data/exp1-15-servers/"));
    }
    fclose(fpResult);
}

void ur_without_error_15_m() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_15_m.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int m = 3; m <= 40; m+=1) {
        fprintf(fpResult, "%.4f,", ModelRunner::run_SIS(15, 8000, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.4f,", ModelRunner::run_myAPMISRR(15, 0.5, m, 8000, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.4f\n", ModelRunner::run_MISRR(15, m, 8000, 0.3, "../data/exp1-15-servers/"));
    }
    fclose(fpResult);
}

void without_error_15_W() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_15_W.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 25000; W+=1000) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(15, W, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(15, 0.5, 24, W, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.2f\n", ModelRunner::run_MISRR(15, 24, W, 0.3, "../data/exp1-15-servers/"));
    }
    fclose(fpResult);
}

void ur_without_error_15_W() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_15_W.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 25000; W+=1000) {
        fprintf(fpResult, "%.4f,", ModelRunner::run_SIS(15, W, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.4f,", ModelRunner::run_myAPMISRR(15, 0.5, 24, W, 0.3, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.4f\n", ModelRunner::run_MISRR(15, 24, W, 0.3, "../data/exp1-15-servers/"));
    }
    fclose(fpResult);
}

void without_error_15_theta() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_15_theta.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(double theta = 0.1; theta <= 1; theta+=0.1) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(15, 8000, theta, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(15, 0.5, 24, 8000, theta, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.2f\n", ModelRunner::run_MISRR(15, 24, 8000, theta, "../data/exp1-15-servers/"));
    }
    fclose(fpResult);
}

void ur_without_error_15_theta() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_15_theta.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(double theta = 0.1; theta <= 1; theta+=0.1) {
        fprintf(fpResult, "%.4f,", ModelRunner::run_SIS(15, 8000, theta, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.4f,", ModelRunner::run_myAPMISRR(15, 0.5, 24, 8000, theta, "../data/exp1-15-servers/"));
        fprintf(fpResult, "%.4f\n", ModelRunner::run_MISRR(15, 24, 8000, theta, "../data/exp1-15-servers/"));
    }
    fclose(fpResult);
}

void without_error_30_m() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_30_m.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int m = 3; m <= 40; m+=1) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(30, 8000, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(30, 0.5, m, 8000, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.2f\n", ModelRunner::run_MISRR(30, m, 8000, 0.3, "../data/exp1-30-servers/"));
    }
    fclose(fpResult);
}

void ur_without_error_30_m() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_30_m.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int m = 3; m <= 40; m+=1) {
        fprintf(fpResult, "%.4f,", ModelRunner::run_SIS(30, 8000, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.4f,", ModelRunner::run_myAPMISRR(30, 0.5, m, 8000, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.4f\n", ModelRunner::run_MISRR(30, m, 8000, 0.3, "../data/exp1-30-servers/"));
    }
    fclose(fpResult);
}

void without_error_30_W() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_30_W.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 15000; W+=1000) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(30, W, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(30, 0.5, 24, W, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.2f\n", ModelRunner::run_MISRR(30, 0, W, 0.3, "../data/exp1-30-servers/"));
    }
    fclose(fpResult);
}

void ur_without_error_30_W() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_30_W.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 25000; W+=1000) {
        fprintf(fpResult, "%.4f,", ModelRunner::run_SIS(30, W, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.4f,", ModelRunner::run_myAPMISRR(30, 0.5, 24, W, 0.3, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.4f\n", ModelRunner::run_MISRR(30, 24, W, 0.3, "../data/exp1-30-servers/"));
    }
    fclose(fpResult);
}

void without_error_30_theta() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/without_error_30_theta.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(double theta = 0.1; theta <= 1; theta+=0.1) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(30, 8000, theta, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(30, 0.5, 24, 8000, theta, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.2f\n", ModelRunner::run_MISRR(30, 24, 8000, theta, "../data/exp1-30-servers/"));
    }
    fclose(fpResult);
}

void ur_without_error_30_theta() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/ur_without_error_30_theta.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(double theta = 0.1; theta <= 1; theta+=0.1) {
        fprintf(fpResult, "%.4f,", ModelRunner::run_SIS(30, 8000, theta, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.4f,", ModelRunner::run_myAPMISRR(30, 0.5, 24, 8000, theta, "../data/exp1-30-servers/"));
        fprintf(fpResult, "%.4f\n", ModelRunner::run_MISRR(30, 24, 8000, theta, "../data/exp1-30-servers/"));
    }
    fclose(fpResult);
}

void with_error_15() {
    FILE * fpResult;
    fpResult = fopen("../output/exp_1/with_error_15.csv", "w");
    fprintf(fpResult, "SIS,APMISRR,toler-MIS\n");
    for(int W = 5000; W <= 25000; W+=1000) {
        fprintf(fpResult, "%.2f,", ModelRunner::run_SIS(15, W, 0.3, "../data/15-servers-w-20/", {5}));
        fprintf(fpResult, "%.2f,", ModelRunner::run_myAPMISRR(15, 0.5, 24, W, 0.3, "../data/15-servers-w-20/", {5}));
        fprintf(fpResult, "%.2f", ModelRunner::run_MISRR(15, 24, W, 0.3, "../data/15-servers-w-20/", {5}));
        fprintf(fpResult, "\n");
    }
    fclose(fpResult);
}