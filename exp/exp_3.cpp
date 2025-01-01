/*
 * @FilePath: /InstallmentScheduling/exp/exp_3.cpp
 * @Description: 带容错的3个模型实验
 * @Author: rthete
 * @Date: 2023-08-19 17:41:41
 * @LastEditTime: 2023-11-18 18:33:53
 */
#include "exp_3.h"

namespace exp_3 {
// 从txt文件中读取数据
void readTXTFile(const std::string &file_path,
                 std::vector<std::vector<int>> &error_place) {
  std::ifstream file(file_path);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      std::vector<int> row;
      std::stringstream ss(line);
      std::string num;
      while (std::getline(ss, num, ',')) {
        row.push_back(std::stoi(num));
      }
      error_place.push_back(row);
    }
    file.close();
  } else {
    std::cerr << "Error opening file for reading: " << file_path << std::endl;
  }
}

// 表格数据调换xy轴
void transposeCSV(const std::string &inputFilePath,
                  const std::string &outputFilePath) {
  std::ifstream inputFile(inputFilePath);
  std::ofstream outputFile(outputFilePath);

  if (!inputFile.is_open()) {
    std::cerr << "Error opening input file: " << inputFilePath << std::endl;
    return;
  }

  if (!outputFile.is_open()) {
    std::cerr << "Error opening output file: " << outputFilePath << std::endl;
    return;
  }

  std::vector<std::vector<std::string>> matrix;
  std::string line;

  // Read input CSV and store it as a matrix
  while (std::getline(inputFile, line)) {
    std::vector<std::string> row;
    std::istringstream ss(line);
    std::string cell;

    while (std::getline(ss, cell, ',')) {
      row.push_back(cell);
    }

    matrix.push_back(row);
  }

  // Transpose the matrix
  if (!matrix.empty()) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    for (size_t col = 0; col < cols; ++col) {
      for (size_t row = 0; row < rows; ++row) {
        outputFile << matrix[row][col];
        if (row < rows - 1) {
          outputFile << ",";
        }
      }
      outputFile << std::endl;
    }
  }

  inputFile.close();
  outputFile.close();
}

/*冲突时间实验v2 横坐标任务量*/
void error_TolerMIS_30_conflict_workload() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;

  error_places.clear();
  readTXTFile("../data/exp3-error-conflict/error-place-30-conflict.txt",
              error_places);
  readTXTFile("../data/exp3-error-conflict/error-installment-conflict.txt",
              error_installment);

  auto i = 1;
  for (const auto &row : error_places) {
    std::ofstream conflictFile(
        "../output/exp_3_conflict/error_TolerMIS_30_conflict_workload_n_" +
        std::to_string(i++) + ".csv");

    // 5000任务量，总趟数为19
    for (auto workload = 5000; workload < 16000; workload += 1000) {
      vector<double> result;
      // 获取每个处理机的冲突时间
      ModelRunner::run_MISRR_conflict(result, 30, 0, workload, 0.3,
                                      "../data/exp1-30-servers/", row, 18);
      for (auto it = result.begin(); it != result.end(); ++it) {
        conflictFile << *it;
        // 如果不是最后一个元素，添加逗号
        if (std::next(it) != result.end()) {
          conflictFile << ",";
        }
      }
      conflictFile << std::endl;
    }
    conflictFile.close();
  }
}

/*冲突时间实验v2 横坐标趟数*/
void error_TolerMIS_30_conflict_inst() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;

  error_places.clear();
  readTXTFile("../data/exp3-error-conflict/error-place-30-conflict.txt",
              error_places);
  readTXTFile("../data/exp3-error-conflict/error-installment-conflict.txt",
              error_installment);

  auto i = 1;
  for (const auto &row : error_places) {
    std::ofstream conflictFile(
        "../output/exp_3_conflict/error_TolerMIS_30_conflict_n_" +
        std::to_string(i++) + ".csv");

    // 5000任务量，总趟数为19
    for (auto inst = 2; inst < 19; inst++) {
      vector<double> result;
      // 获取每个处理机的冲突时间
      ModelRunner::run_MISRR_conflict(result, 30, 0, 5000, 0.3,
                                      "../data/exp1-30-servers/", row, inst);
      for (auto it = result.begin(); it != result.end(); ++it) {
        conflictFile << *it;
        // 如果不是最后一个元素，添加逗号
        if (std::next(it) != result.end()) {
          conflictFile << ",";
        }
      }
      conflictFile << std::endl;
    }
    conflictFile.close();
  }
}

/*冲突时间实验*/
void error_TolerMIS_30_conflict() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;
  std::ofstream conflictFile(
      "../output/exp_3_conflict/error_TolerMIS_30_conflict_cmp.csv");
  error_places.clear();
  readTXTFile("../data/exp3-error-conflict/error-place-30-conflict.txt",
              error_places);
  readTXTFile("../data/exp3-error-conflict/error-installment-conflict.txt",
              error_installment);

  /* 每种30次取平均 */
  /* 任务量5000/10000 */
  for (const auto &row : error_places) {
    vector<double> result;
    ModelRunner::run_MISRR_conflict(result, 30, 0, 5000, 0.3,
                                    "../data/exp1-30-servers/", row,
                                    error_installment[0][0]);
    // ModelRunner::run_MISRR_conflict(result, 30, 0, 10000, 0.3,
    // "../data/exp1-30-servers/",
    //                    row, error_installment[1][0]);
    // 对比没有冲突情况
    ModelRunner::run_MISRR_conflict(result, 30, 0, 5000, 0.3,
                                    "../data/exp1-30-servers/", row,
                                    error_installment[2][0]);
    // ModelRunner::run_MISRR_conflict(result, 30, 0, 10000, 0.3,
    // "../data/exp1-30-servers/",
    //                    row, error_installment[3][0]);
    for (auto it = result.begin(); it != result.end(); ++it) {
      conflictFile << *it;
      // 如果不是最后一个元素，添加逗号
      if (std::next(it) != result.end()) {
        conflictFile << ",";
      }
    }
    conflictFile << std::endl;
  }

  conflictFile.close();
}

/* SIS模型，共30个处理机，故障1/2/3/4个处理机，每种实验30次 */
void error_SIS_30() {
  std::vector<std::vector<int>> error_places;
  std::ofstream meanFile("../output/exp_3_time/error_SIS_30_mean.csv");
  std::ofstream diffFile("../output/exp_3_time/error_SIS_30_diff.csv");
  std::ofstream maxFile("../output/exp_3_time/error_SIS_30_max.csv");
  std::ofstream minFile("../output/exp_3_time/error_SIS_30_min.csv");

  /* 故障i个处理机 */
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) +
                    ".txt",
                error_places);
    /* 任务量w */
    for (int w = 5000; w < 16000; w += 1000) {
      double result_sum = 0;
      double result_min = 100000, result_max = 0;

      /* 每种30次取平均 */
      for (const auto &row : error_places) {
        double result =
            ModelRunner::run_SIS(30, w, 0.3, "../data/exp1-30-servers/", row);
        result_sum += result;
        result_min = std::min(result_min, result);
        result_max = std::max(result_max, result);
      }
      meanFile << result_sum / 30 << ",";
      diffFile << result_max - result_min << ",";
      maxFile << result_max << ",";
      minFile << result_min << ",";
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
  //   transposeCSV("../output/exp_3/error_SIS_30_mean.csv",
  //                "../output/exp_3/t_error_SIS_30_mean.csv");
  //   transposeCSV("../output/exp_3/error_SIS_30_diff.csv",
  //                "../output/exp_3/t_error_SIS_30_diff.csv");
}

/* TolerMIS模型，共30个处理机，故障1/2/3/4个处理机，每种实验30次 */
void error_TolerMIS_30() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;
  std::ofstream meanFile("../output/exp_3_time/error_TolerMIS_30_mean.csv");
  std::ofstream diffFile("../output/exp_3_time/error_TolerMIS_30_diff.csv");
  std::ofstream maxFile("../output/exp_3_time/error_TolerMIS_30_max.csv");
  std::ofstream minFile("../output/exp_3_time/error_TolerMIS_30_min.csv");

  /* 故障i个处理机 */
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) +
                    ".txt",
                error_places);
    readTXTFile("../data/exp3-error-place/error-installment-2.txt",
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
        double result =
            ModelRunner::run_MISRR(30, 0, w, 0.3, "../data/exp1-30-servers/",
                                   row, error_installment[index_w][index_30++]);
        if (result == 0) {
          not_schedulable_cnt++;
          cout << "!!!!!!!!" << endl;
        }
        result_sum += result;
        result_min = std::min(result_min, result);
        result_max = std::max(result_max, result);
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

/* APMISRR模型，共30个处理机，故障1/2/3/4个处理机，每种实验30次 */
void error_APMISRR_30() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;

  std::ofstream meanFile("../output/exp_3_time/error_APMISRR_30_mean.csv");
  std::ofstream diffFile("../output/exp_3_time/error_APMISRR_30_diff.csv");
  std::ofstream maxFile("../output/exp_3_time/error_APMISRR_30_max.csv");
  std::ofstream minFile("../output/exp_3_time/error_APMISRR_30_min.csv");

  /* 故障i个处理机 */
  vector<int> m = {19, 21, 23, 25, 26, 28, 29, 30, 31, 33, 34};
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) +
                    ".txt",
                error_places);
    readTXTFile("../data/exp3-error-place/error-installment-2.txt",
                error_installment);

    /* 任务量w */
    int index_w = 0;
    for (int w = 5000; w < 16000; w += 1000) {
      int index_30 = 0;
      double result_sum = 0;
      double result_min = 100000, result_max = 0;
      /* 每种30次取平均 */
      for (const auto &row : error_places) {
        double result = ModelRunner::run_myAPMISRR(
            30, 0.5, m[index_w], w, 0.3, "../data/exp1-30-servers/", row,
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
  //   transposeCSV("../output/exp_3/error_APMISRR_30_mean.csv",
  //                "../output/exp_3/t_error_APMISRR_30_mean.csv");
  //   transposeCSV("../output/exp_3/error_APMISRR_30_diff.csv",
  //                "../output/exp_3/t_error_APMISRR_30_diff.csv");
}

/*利用率，3个模型，共30个处理机，故障1/2/3/4个处理机，任务量5000，每种实验30次
 */
void error_cmp_3_models_ur() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;
  vector<int> m = {19, 21, 23, 25, 26, 28, 29, 30, 31, 33, 34};

  std::ofstream allFile("../output/exp_3_all_exps/error_ur_s30_w5000.csv");
  allFile << "model,error_num,value" << std::endl;

  /* 故障i个处理机 */
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) +
                    ".txt",
                error_places);
    readTXTFile("../data/exp3-error-place/error-installment-2.txt",
                error_installment);
    // std::ofstream singleFile("../output/exp_3_all_exps/error_ur_30_" +
    //                          std::to_string(i) + ".csv");
    /* 任务量w */
    int index_w = 0;
    for (int w = 5000; w < 5500; w += 1000) {
      int index_30 = 0;

      /* 每种30次 */
      for (const auto &row : error_places) {
        double result =
            ModelRunner::run_SIS(30, w, 0.3, "../data/exp1-30-servers/", row);
        // singleFile << result << ",";
        allFile << "SIS," << i << "," << result << std::endl;

        result = ModelRunner::run_myAPMISRR(
            30, 0.5, m[index_w], w, 0.3, "../data/exp1-30-servers/", row,
            error_installment[index_w][index_30]);
        // singleFile << result << "\n";
        allFile << "APMISRR," << i << "," << result << std::endl;

        result =
            ModelRunner::run_MISRR(30, 0, w, 0.3, "../data/exp1-30-servers/",
                                   row, error_installment[index_w][index_30]);
        // singleFile << result << ",";
        allFile << "TolerMIS," << i << "," << result << std::endl;

        index_30++;
      }
      index_w++;
    }
    // singleFile.close();
  }
}

/* 时间，3个模型，共30个处理机，故障1/2/3/4个处理机，任务量5000~15000，每种实验30次
 */
void error_save_each_time() {
  std::vector<std::vector<int>> error_places;
  std::vector<std::vector<int>> error_installment;
  vector<int> m = {19, 21, 23, 25, 26, 28, 29, 30, 31, 33, 34};

  /* 故障i个处理机 */
  for (int i = 1; i <= 4; ++i) {
    error_places.clear();
    readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) +
                    ".txt",
                error_places);
    readTXTFile("../data/exp3-error-place/error-installment-2.txt",
                error_installment);

    std::ofstream allFile("../output/exp_3_all_exps/error_time_s30_n" +
                          std::to_string(i) + ".csv");
    allFile << "model,workload,value" << std::endl;

    /* 任务量w */
    int index_w = 0;
    for (int w = 5000; w < 16000; w += 1000) {
      int index_30 = 0;
      double result_sum = 0;
      double result_min = 100000, result_max = 0;
      // std::ofstream singleFile("../output/exp_3_all_exps/error_time_s30_n" +
      //                          std::to_string(i) + "_W" + std::to_string(w) +
      //                          ".csv");
      // singleFile << "SIS,APMISRR,TolerMIS\n";

      /* 每种30次取平均 */
      for (const auto &row : error_places) {
        double result =
            ModelRunner::run_SIS(30, w, 0.3, "../data/exp1-30-servers/", row);
        allFile << "SIS," << w << "," << result << std::endl;
        // singleFile << result << ",";
        result = ModelRunner::run_myAPMISRR(
            30, 0.5, m[index_w], w, 0.3, "../data/exp1-30-servers/", row,
            error_installment[index_w][index_30]);
        allFile << "APMISRR," << w << "," << result << std::endl;
        // singleFile << result << ",";
        result =
            ModelRunner::run_MISRR(30, 0, w, 0.3, "../data/exp1-30-servers/",
                                   row, error_installment[index_w][index_30++]);
        allFile << "TolerMIS," << w << "," << result << std::endl;
        // singleFile << result << "\n";
        result_sum += result;
        result_min = std::min(result_min, result);
        result_max = std::max(result_max, result);
      }

      index_w++;
      // singleFile.close();
    }
  }
}
} // namespace exp_3