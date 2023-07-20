/*
 * @FilePath: \InstallmentScheduling\exp\exp_2.cpp
 * @Description: TolerMIS单模型对比
 * @Author: rthete
 * @Date: 2023-07-09 15:01:22
 * @LastEditTime: 2023-07-19 17:38:24
 */
#include "exp_2.h"

namespace exp_2{
    void W_vs_m_15() {
        FILE * fpResult;
        fpResult = fopen("../output/exp_2/W_vs_m_15.csv", "w");
        fprintf(fpResult, "m,");
        for(int W = 1000; W <= 5000; W+=1000) {
            fprintf(fpResult, "%d,", W);
        }
        fprintf(fpResult, "\n");
        for(int m = 3; m <= 40; m+=1) {
            fprintf(fpResult, "%d,", m);
            for(int W = 1000; W <= 5000; W+=1000) {
                fprintf(fpResult, "%.2f,", run_MISRR(15, m, W, 0.3, "../data/exp1-15-servers/"));
            }
            fprintf(fpResult, "\n");
        }
    }

    void theta_vs_m_15() {
        FILE * fpResult;
        fpResult = fopen("../output/exp_2/theta_vs_m_15.csv", "w");
        fprintf(fpResult, "m,");
        for(double theta = 0.1; theta <= 1; theta+=0.1) {
            fprintf(fpResult, "%.1f,", theta);
        }
        fprintf(fpResult, "\n");
        for(int m = 3; m <= 40; m+=1) {
            fprintf(fpResult, "%d,", m);
            for(double theta = 0.1; theta <= 1; theta+=0.1) {
                fprintf(fpResult, "%.2f,", run_MISRR(15, m, 5000, theta, "../data/exp1-15-servers/"));
            }
            fprintf(fpResult, "\n");
        }
    }

    void W_vs_m_30() {
        FILE * fpResult;
        fpResult = fopen("../output/exp_2/W_vs_m_30.csv", "w");
        fprintf(fpResult, "m,");
        for(int W = 1000; W <= 5000; W+=1000) {
            fprintf(fpResult, "%d,", W);
        }
        fprintf(fpResult, "\n");
        for(int m = 3; m <= 40; m+=1) {
            fprintf(fpResult, "%d,", m);
            for(int W = 1000; W <= 5000; W+=1000) {
                fprintf(fpResult, "%.2f,", run_MISRR(30, m, W, 0.3, "../data/exp1-30-servers/"));
            }
            fprintf(fpResult, "\n");
        }
    }

    void theta_vs_m_30() {
        FILE * fpResult;
        fpResult = fopen("../output/exp_2/theta_vs_m_30.csv", "w");
        fprintf(fpResult, "m,");
        for(double theta = 0.1; theta <= 1; theta+=0.1) {
            fprintf(fpResult, "%.1f,", theta);
        }
        fprintf(fpResult, "\n");
        for(int m = 3; m <= 40; m+=1) {
            fprintf(fpResult, "%d,", m);
            for(double theta = 0.1; theta <= 1; theta+=0.1) {
                fprintf(fpResult, "%.2f,", run_MISRR(30, m, 5000, theta, "../data/exp1-30-servers/"));
            }
            fprintf(fpResult, "\n");
        }
    }
} // namspace exp_2
