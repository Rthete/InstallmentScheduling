/*
 * @FilePath: \InstallmentScheduling\exp\exp_3.cpp
 * @Description: 带容错的3个模型实验
 * @Author: rthete
 * @Date: 2023-08-19 17:41:41
 * @LastEditTime: 2023-09-09 13:30:25
 */
#include "exp_3.h"

namespace exp_3{
    // 从txt文件中读取数据
    void readTXTFile(const std::string& file_path, std::vector<std::vector<int>>& error_place) {
        std::ifstream file(file_path);
        if(file.is_open()) {
            std::string line;
            while(std::getline(file, line)) {
                std::vector<int> row;
                std::stringstream ss(line);
                std::string num;
                while(std::getline(ss, num, ',')) {
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
    void transposeCSV(const std::string& inputFilePath, const std::string& outputFilePath) {
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

    // SIS模型，共30个处理机，故障1/2/3/4个处理机，每种实验30次
    void error_SIS_30() {
        std::vector<std::vector<int>> error_places;
        std::ofstream outputFile("../output/exp_3/error_SIS_30.csv");
        for(int i = 1; i <= 4; ++i) {
            error_places.clear();
            readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) + ".txt", error_places);
            
            for (const auto& row : error_places) {
                double result = run_SIS(30, 8000, 0.3, "../data/exp1-30-servers/", row);
                outputFile << result << ",";
            }
            outputFile << std::endl;
            
        }
        outputFile.close();
        transposeCSV("../output/exp_3/error_SIS_30.csv", "../output/exp_3/error_SIS_30_transpose.csv");
        
    }

    void error_TolerMIS_30() {
        std::vector<std::vector<int>> error_places;
        std::vector<std::vector<int>> error_installment;

        std::ofstream outputFile("../output/exp_3/error_TolerMIS_30.csv");
        for(int i = 1; i <= 4; ++i) {
            error_places.clear();
            readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) + ".txt", error_places);
            readTXTFile("../data/exp3-error-place/error-installment.txt", error_installment);
            int index = 0;
            for (const auto& row : error_places) {
                double result = run_MISRR(30, 24, 8000, 0.3, "../data/exp1-30-servers/", row, error_installment[index++][0]);
                outputFile << result << ",";
            }
            outputFile << std::endl;
        }
        outputFile.close();
        transposeCSV("../output/exp_3/error_TolerMIS_30.csv", "../output/exp_3/error_TolerMIS_30_transpose.csv");
    }

    void error_APMISRR_30() {
        std::vector<std::vector<int>> error_places;
        std::vector<std::vector<int>> error_installment;

        std::ofstream outputFile("../output/exp_3/error_APMISRR_30.csv");
        for(int i = 1; i <= 4; ++i) {
            error_places.clear();
            readTXTFile("../data/exp3-error-place/error-place-30-" + std::to_string(i) + ".txt", error_places);
            readTXTFile("../data/exp3-error-place/error-installment.txt", error_installment);
            int index = 0;
            for (const auto& row : error_places) {
                double result = run_myAPMISRR(30, 0.5, 24, 8000, 0.3, "../data/exp1-30-servers/", row, error_installment[index++][0]);
                outputFile << result << ",";
            }
            outputFile << std::endl;
        }
        outputFile.close();
        transposeCSV("../output/exp_3/error_APMISRR_30.csv", "../output/exp_3/error_APMISRR_30_transpose.csv");
    }
}