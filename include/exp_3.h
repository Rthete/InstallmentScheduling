/*
 * @FilePath: /InstallmentScheduling/include/exp_3.h
 * @Description:
 * @Author: rthete
 * @Date: 2023-08-19 17:44:04
 * @LastEditTime: 2023-11-18 18:34:33
 */
#ifndef _EXP_3_H
#define _EXP_3_H

#include "method.h"
#include <algorithm>
#include <fstream>
#include <sstream>

namespace exp_3 {
void transposeCSV(const std::string &inputFilePath,
                  const std::string &outputFilePath);
void error_SIS_30();
void error_TolerMIS_30();
void error_APMISRR_30();
void error_TolerMIS_30_conflict();
void error_cmp_3_models_ur();
void error_save_each_time();
void error_TolerMIS_30_conflict_inst();
void error_TolerMIS_30_conflict_workload();

// 从txt文件中读取数据
template <typename T>
void readTXTFile(const std::string &file_path,
                 std::vector<std::vector<T>> &file_data) {
  std::ifstream file(file_path);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      std::vector<T> row;
      std::stringstream ss(line);
      std::string num;
      while (std::getline(ss, num, ',')) {
        if constexpr (std::is_same<T, double>::value) {
          row.push_back(std::stod(num));
        } else {
          row.push_back(std::stoi(num));
        }
      }
      file_data.push_back(row);
    }
    file.close();
  } else {
    std::cerr << "Error opening file for reading: " << file_path << std::endl;
  }
}
} // namespace exp_3

#endif //_EXP_3_H