#include "exp_4.h"
#include "exp_3.h"

#include "method.h"

namespace exp_4 {
void without_error_15_m() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_without_error_15_m.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_myAPMISRR(15, 0.5, m, 8000, 0.3,
                                       "../data/exp1-15-servers/"));
    fprintf(
        fpResult, "%.2f,",
        ModelRunner::run_MISRR(15, m, 8000, 0.3, "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.2f\n",
            ModelRunner::run_MISRRLL(15, 0.2, 0.3, m, 8000, 0.3,
                                     "../data/exp1-15-servers/"));
  }
  fclose(fpResult);
}

void ur_without_error_15_m() {
  ModelRunner::output_using_rate = 1;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_ur_without_error_15_m.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(15, 0.5, m, 8000, 0.3,
                                       "../data/exp1-15-servers/"));
    fprintf(
        fpResult, "%.4f,",
        ModelRunner::run_MISRR(15, m, 8000, 0.3, "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(15, 0.2, 0.3, m, 8000, 0.3,
                                     "../data/exp1-15-servers/"));
  }
  fclose(fpResult);
}

void without_error_30_m() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_without_error_30_m.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_myAPMISRR(30, 0.5, m, 8000, 0.3,
                                       "../data/exp1-30-servers/"));
    fprintf(
        fpResult, "%.2f,",
        ModelRunner::run_MISRR(30, m, 8000, 0.3, "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.2f\n",
            ModelRunner::run_MISRRLL(30, 0.2, 0.3, m, 8000, 0.3,
                                     "../data/exp1-30-servers/"));
  }
  fclose(fpResult);
}

void ur_without_error_30_m() {
  ModelRunner::output_using_rate = 1;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_ur_without_error_30_m.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(30, 0.5, m, 8000, 0.3,
                                       "../data/exp1-30-servers/"));
    fprintf(
        fpResult, "%.4f,",
        ModelRunner::run_MISRR(30, m, 8000, 0.3, "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(30, 0.2, 0.3, m, 8000, 0.3,
                                     "../data/exp1-30-servers/"));
  }
  fclose(fpResult);
}

void without_error_15_W() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_without_error_15_W.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  // for(int W = 5000; W <= 25000; W+=1000) { // 画图差值太小
  for (int W = 5000; W <= 5100; W += 10) {
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_myAPMISRR(15, 0.5, 24, W, 0.3,
                                       "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_MISRR(15, 24, W, 0.3, "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.2f\n",
            ModelRunner::run_MISRRLL(15, 0.2, 0.3, 24, W, 0.3,
                                     "../data/exp1-15-servers/"));
  }
  fclose(fpResult);
}

void ur_without_error_15_W() {
  ModelRunner::output_using_rate = 1;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_ur_without_error_15_W.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int W = 5000; W <= 25000; W += 1000) {
    // for(int W = 5000; W <= 6000; W+=100) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(15, 0.5, 24, W, 0.3,
                                       "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_MISRR(15, 24, W, 0.3, "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(15, 0.2, 0.3, 24, W, 0.3,
                                     "../data/exp1-15-servers/"));
  }
  fclose(fpResult);
}

void without_error_30_W() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_without_error_30_W.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int W = 5000; W <= 25000; W += 1000) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(30, 0.5, 24, W, 0.3,
                                       "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_MISRR(30, 24, W, 0.3, "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(30, 0.2, 0.3, 24, W, 0.3,
                                     "../data/exp1-30-servers/"));
  }
  fclose(fpResult);
}

void ur_without_error_30_W() {
  ModelRunner::output_using_rate = 1;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_ur_without_error_30_W.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (int W = 5000; W <= 25000; W += 1000) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(30, 0.5, 24, W, 0.3,
                                       "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_MISRR(30, 24, W, 0.3, "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(30, 0.2, 0.3, 24, W, 0.3,
                                     "../data/exp1-30-servers/"));
  }
  fclose(fpResult);
}

void without_error_15_theta() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_without_error_15_theta.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (double theta = 0.1; theta <= 1; theta += 0.1) {
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_myAPMISRR(15, 0.5, 24, 8000, theta,
                                       "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_MISRR(15, 24, 8000, theta,
                                   "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.2f\n",
            ModelRunner::run_MISRRLL(15, 0.2, 0.3, 24, 8000, theta,
                                     "../data/exp1-15-servers/"));
  }
  fclose(fpResult);
}

void ur_without_error_15_theta() {
  ModelRunner::output_using_rate = 1;
  FILE *fpResult;
  fpResult =
      fopen("../output/exp_4/MISRRLL_ur_without_error_15_theta.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (double theta = 0.1; theta <= 1; theta += 0.1) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(15, 0.5, 24, 8000, theta,
                                       "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_MISRR(15, 24, 8000, theta,
                                   "../data/exp1-15-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(15, 0.2, 0.3, 24, 8000, theta,
                                     "../data/exp1-15-servers/"));
  }
  fclose(fpResult);
}

void without_error_30_theta() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_without_error_30_theta.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (double theta = 0.1; theta <= 1; theta += 0.1) {
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_myAPMISRR(30, 0.5, 24, 8000, theta,
                                       "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.2f,",
            ModelRunner::run_MISRR(30, 24, 8000, theta,
                                   "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.2f\n",
            ModelRunner::run_MISRRLL(30, 0.2, 0.3, 24, 8000, theta,
                                     "../data/exp1-30-servers/"));
  }
  fclose(fpResult);
}

void ur_without_error_30_theta() {
  ModelRunner::output_using_rate = 1;
  FILE *fpResult;
  fpResult =
      fopen("../output/exp_4/MISRRLL_ur_without_error_30_theta.csv", "w");
  fprintf(fpResult, "APMISRR,tolerMIS,MISRRLL\n");
  for (double theta = 0.1; theta <= 1; theta += 0.1) {
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_myAPMISRR(30, 0.5, 24, 8000, theta,
                                       "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f,",
            ModelRunner::run_MISRR(30, 24, 8000, theta,
                                   "../data/exp1-30-servers/"));
    fprintf(fpResult, "%.4f\n",
            ModelRunner::run_MISRRLL(30, 0.2, 0.3, 24, 8000, theta,
                                     "../data/exp1-30-servers/"));
  }
  fclose(fpResult);
}

void W_vs_m_15() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_W_vs_m_15.csv", "w");
  fprintf(fpResult, "m,");
  for (int W = 1000; W <= 5000; W += 1000) {
    fprintf(fpResult, "%d,", W);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%d,", m);
    for (int W = 1000; W <= 5000; W += 1000) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(15, 0.2, 0.3, m, W, 0.3,
                                       "../data/exp1-15-servers/"));
    }
    fprintf(fpResult, "\n");
  }
}

void theta_vs_m_15() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_theta_vs_m_15.csv", "w");
  fprintf(fpResult, "m,");
  for (double theta = 0.1; theta <= 1; theta += 0.1) {
    fprintf(fpResult, "%.1f,", theta);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%d,", m);
    for (double theta = 0.1; theta <= 1; theta += 0.1) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(15, 0.2, 0.3, m, 5000, theta,
                                       "../data/exp1-15-servers/"));
    }
    fprintf(fpResult, "\n");
  }
}

void W_vs_m_30() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_W_vs_m_30.csv", "w");
  fprintf(fpResult, "m,");
  for (int W = 1000; W <= 5000; W += 1000) {
    fprintf(fpResult, "%d,", W);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%d,", m);
    for (int W = 1000; W <= 5000; W += 1000) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(30, 0.2, 0.3, m, W, 0.3,
                                       "../data/exp1-30-servers/"));
    }
    fprintf(fpResult, "\n");
  }
}

void theta_vs_m_30() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_theta_vs_m_30.csv", "w");
  fprintf(fpResult, "m,");
  for (double theta = 0.1; theta <= 1; theta += 0.1) {
    fprintf(fpResult, "%.1f,", theta);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m += 1) {
    fprintf(fpResult, "%d,", m);
    for (double theta = 0.1; theta <= 1; theta += 0.1) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(30, 0.2, 0.3, m, 5000, theta,
                                       "../data/exp1-30-servers/"));
    }
    fprintf(fpResult, "\n");
  }
}

void labmda1_T_vs_m_15() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_labmda1_T_vs_m_15.csv", "w");
  fprintf(fpResult, "m,");
  for (double lambda1 = 0.1; lambda1 < 1; lambda1 += 0.1) {
    fprintf(fpResult, "%.1f,", lambda1);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m++) {
    fprintf(fpResult, "%d,", m);
    for (double lambda1 = 0.1; lambda1 < 1; lambda1 += 0.1) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(15, lambda1, 0.3, m, 8000, 0.3,
                                       "../data/exp1-15-servers/"));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void labmda2_T_vs_m_15() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_labmda2_T_vs_m_15.csv", "w");
  fprintf(fpResult, "m,");
  for (double lambda2 = 0.1; lambda2 < 1; lambda2 += 0.1) {
    fprintf(fpResult, "%.1f,", lambda2);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m++) {
    fprintf(fpResult, "%d,", m);
    for (double lambda2 = 0.1; lambda2 < 1; lambda2 += 0.1) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(15, 0.2, lambda2, m, 8000, 0.3,
                                       "../data/exp1-15-servers/"));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void labmda1_T_vs_m_30() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_labmda1_T_vs_m_30.csv", "w");
  fprintf(fpResult, "m,");
  for (double lambda1 = 0.1; lambda1 < 1; lambda1 += 0.1) {
    fprintf(fpResult, "%.1f,", lambda1);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m++) {
    fprintf(fpResult, "%d,", m);
    for (double lambda1 = 0.1; lambda1 < 1; lambda1 += 0.1) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(30, lambda1, 0.3, m, 8000, 0.3,
                                       "../data/exp1-30-servers/"));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void labmda2_T_vs_m_30() {
  ModelRunner::output_using_rate = 0;
  FILE *fpResult;
  fpResult = fopen("../output/exp_4/MISRRLL_labmda2_T_vs_m_30.csv", "w");
  fprintf(fpResult, "m,");
  for (double lambda2 = 0.1; lambda2 < 1; lambda2 += 0.1) {
    fprintf(fpResult, "%.1f,", lambda2);
  }
  fprintf(fpResult, "\n");
  for (int m = 3; m <= 40; m++) {
    fprintf(fpResult, "%d,", m);
    for (double lambda2 = 0.1; lambda2 < 1; lambda2 += 0.1) {
      fprintf(fpResult, "%.2f,",
              ModelRunner::run_MISRRLL(30, 0.2, lambda2, m, 8000, 0.3,
                                       "../data/exp1-30-servers/"));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void error_MISRRLL_30() {
  ModelRunner::output_using_rate = 0;

  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;
  std::ofstream meanFile("../output/exp_4_error/error_MISRRLL_30_mean.csv");
  std::ofstream diffFile("../output/exp_4_error/error_MISRRLL_30_diff.csv");
  std::ofstream maxFile("../output/exp_4_error/error_MISRRLL_30_max.csv");
  std::ofstream minFile("../output/exp_4_error/error_MISRRLL_30_min.csv");

  /* 故障i个处理机 */
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    exp_3::readTXTFile("../data/exp3-error-place/error-place-30-" +
                           std::to_string(i) + ".txt",
                       error_places);
    exp_3::readTXTFile("../data/exp3-error-place/error-installment-2.txt",
                       error_installment);

    /* 任务量w */
    int index_w = 0;
    for (int w = 5000; w < 16000; w += 1000) {
      int index_30 = 0;
      double result_sum = 0;
      double result_min = 100000, result_max = 0;

      /* 每种30次取平均 */
      int not_schedulable_cnt = 0;
      for (const auto &row : error_places) {
        double result = ModelRunner::run_MISRRLL(
            30, 0.8, 0.7, 0, w, 0.3, "../data/exp1-30-servers/", row,
            error_installment[index_w][index_30++]);
        if (result == 0)
          not_schedulable_cnt++;
        result_sum += result;
        if (result != 0)
          result_min = std::min(result_min, result);
        result_max = std::max(result_max, result);
        // return;
      }
      meanFile << result_sum / (30 - not_schedulable_cnt) << ",";
      diffFile << result_max - result_min << ",";
      maxFile << result_max << ",";
      minFile << result_min << ",";

      index_w++;
    }
    meanFile << std::endl;
    diffFile << std::endl;
    maxFile << std::endl;
    minFile << std::endl;
  }
  meanFile.close();
  diffFile.close();
  maxFile.close();
  minFile.close();
  //   transposeCSV("../output/exp_3/error_TolerMIS_30_mean.csv",
  //                "../output/exp_3/t_error_TolerMIS_30_mean.csv");
  //   transposeCSV("../output/exp_3/error_TolerMIS_30_diff.csv",
  //                "../output/exp_3/t_error_TolerMIS_30_diff.csv");
}

void ur_error_MISRRLL_30() {
  ModelRunner::output_using_rate = 1;

  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;
  std::ofstream meanFile("../output/exp_4_error/ur_error_MISRRLL_30_mean.csv");
  std::ofstream diffFile("../output/exp_4_error/ur_error_MISRRLL_30_diff.csv");
  std::ofstream maxFile("../output/exp_4_error/ur_error_MISRRLL_30_max.csv");
  std::ofstream minFile("../output/exp_4_error/ur_error_MISRRLL_30_min.csv");

  /* 故障i个处理机 */
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    exp_3::readTXTFile("../data/exp3-error-place/error-place-30-" +
                           std::to_string(i) + ".txt",
                       error_places);
    exp_3::readTXTFile("../data/exp3-error-place/error-installment-2.txt",
                       error_installment);

    /* 任务量w */
    int index_w = 0;
    for (int w = 5000; w < 16000; w += 1000) {
      int index_30 = 0;
      double result_sum = 0;
      double result_min = 100000, result_max = 0;

      /* 每种30次取平均 */
      for (const auto &row : error_places) { // 这里直接用20趟
        double result = ModelRunner::run_MISRRLL(
            30, 0.2, 0.3, 34, w, 0.3, "../data/exp1-30-servers/", row,
            error_installment[index_w][index_30++]);
        result_sum += result;
        result_min = std::min(result_min, result);
        result_max = std::max(result_max, result);
      }
      meanFile << result_sum / 30 << ",";
      diffFile << result_max - result_min << ",";
      maxFile << result_max << ",";
      minFile << result_min << ",";

      index_w++;
    }
    meanFile << std::endl;
    diffFile << std::endl;
    maxFile << std::endl;
    minFile << std::endl;
  }
  meanFile.close();
  diffFile.close();
  maxFile.close();
  minFile.close();
  //   transposeCSV("../output/exp_3/error_TolerMIS_30_mean.csv",
  //                "../output/exp_3/t_error_TolerMIS_30_mean.csv");
  //   transposeCSV("../output/exp_3/error_TolerMIS_30_diff.csv",
  //                "../output/exp_3/t_error_TolerMIS_30_diff.csv");
}
} // namespace exp_4