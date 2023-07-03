/*
 * @FilePath: \InstallmentScheduling\include\PMIS.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-14 15:47:46
 * @LastEditTime: 2023-03-15 15:39:13
 */
//
// Created by xd_wxl on 2022/12/7.
//

#ifndef _PMIS_H
#define _PMIS_H

#include "header.h"
#include "Server.h"

class PMIS {
public:
    PMIS(int valueN, double valueTheta);

    /**
     * assign model
     */
    void setM(int value);
    void setW(double value);
    void initValue();
    void getDataFromFile();
    void getOptimalModel();

    // error
    void error(vector<int> &errorPlace, int errorInstallment);

    // print result
    void printResult();

    // get best installment
    double getTime();
    int getBestInstallment();

    // result
    double getOptimalTime();
    double getUsingRate();

private:
    int n;                                          // 处理机个数
    int m;                                          // 计算趟数
    double theta;                                   // 结果回传比例
    double W = 0;                                   // 总任务量
    double V = 0;                                   // 每一趟调度的任务量
    double optimalTime = 0;                         // 调度最短时间
    double usingRate;                               // 处理机使用率
    vector<Server> servers;                         // 所有的处理器

    // error
    double WWithoutError = 0;                       // 未发生错误时的总任务量
    int numberWithoutError = 0;                     // 未发生错误时的处理机个数
    int beforeInstallment = 0;                      // 未发生错误时的趟数
    double leftW = 0.0;                             // 发生错误剩下的任务量
    double startTime = 0.0;                         // 发生错误时的开始时间

    // alpha && beta
    vector<double> alpha;                           // 内部调度中分配量
    vector<double> beta;                            // 最后一趟调度中分配量

    // time
    vector<double> usingTimes;                      // 每个处理器的处理时间总和

    // print
    string outputName;

    // alpha && beta function
    void initAlpha();
    void initBeta();
};


#endif //_PMIS_H
