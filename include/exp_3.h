/*
 * @FilePath: /InstallmentScheduling/include/exp_3.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-08-19 17:44:04
 * @LastEditTime: 2023-10-21 21:25:46
 */
#ifndef _EXP_3_H
#define _EXP_3_H

#include "method.h"
#include <algorithm>
#include <fstream>
#include <sstream>

namespace exp_3
{
    void readTXTFile(const std::string& filePath, std::vector<std::vector<int>>& error_place);
    void transposeCSV(const std::string& inputFilePath, const std::string& outputFilePath);
    void error_SIS_30();
    void error_TolerMIS_30();
    void error_APMISRR_30();
    void error_TolerMIS_30_conflict();
} // namespace exp_3

#endif //_EXP_3_H