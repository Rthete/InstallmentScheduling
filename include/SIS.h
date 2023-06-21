/*
 * @FilePath: \InstallmentScheduling\include\SIS.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-06-12 19:55:07
 * @LastEditTime: 2023-06-20 21:12:19
 */
//
// Created by kqzhang on 2022/5/13.
//

#ifndef FIRSTMODEL_SIS_H
#define FIRSTMODEL_SIS_H

#include "header.h"
#include "Server.h"

class SIS {
public:
    SIS(int valueN, double valueTheta);

    /**
     * assign model
     */
    void setW(double value);
    void initValue();
    void getDataFromFile(string data_path="../data/w-20/");
    void getOptimalModel();

    // error
    void error(vector<int> &errorPlace);

    // print result
    void printResult();

    // result
    void calUsingRate();
    double getOptimalTime();
    double getUsingRate();

private:
    int n;                                          // 处理机个数
    int m = 0;                                      // 计算趟数
    double theta;                                   // 结果回传比例
    double W;                                       // 总任务量
    double optimalTime;                             // 调度最短时间
    double usingRate;                               // 处理机使用率
    vector<Server> servers;                         // 所有的处理器

    // error
    double WWithoutError;                           // 未发生错误时的总任务量
    int numberWithoutError;                         // 未发生错误时的处理机个数
    double leftW = 0;                               // 发生错误剩下的任务量
    double startTime;                               // 发生错误时的开始时间

    // alpha
    vector<double> alpha;                       // 第一趟调度中分配量

    // time
    vector<double> usingTimes;                  // 每个处理器的处理时间总和

    // print
    string outputName;
};

#endif //MINE_SIS_H

