/*
 * @FilePath: \InstallmentScheduling\include\APMISRR.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-04-22 16:49:36
 * @LastEditTime: 2023-05-06 21:06:35
 */

#ifndef FIRSTMODEL_APMISRR_H
#define FIRSTMODEL_APMISRR_H

#include "header.h"
#include "Server.h"

class APMISRR {
public:
    APMISRR(int valueN, double valueTheta);

    /**
     * assign model
     */
    void setM(int value);
    void setW(double value);
    void setLambda(double value);

    double getAlpha();
    double getBeta();
    void initValue_cost();
    void initValue();
    void getDataFromFile();
    // void getOptimalModel();

    // error
    // void error(vector<int> &errorPlace, int errorInstallment);

    // print result
    // void printResult();

    // get best installment
    // double getTime();
    // int getBestInstallment();

    // result
    double getOptimalTime();
    double getOptimalTime_cost();
    // double getUsingRate();

    // check schedulable
    void isSchedulable();
    int isSchedulable_cost();

private:
    int n;                                          // 处理机个数
    int m;                                          // 计算趟数
    double theta;                                   // 结果回传比例
    double W = 0;                                   // 总任务量
    double V = 0;                                   // 每一趟调度的任务量
    double optimalTime = 0;                         // 调度最短时间
    double usingRate;                               // 处理机使用率
    vector<Server> servers;                         // 所有的处理器

    double lambda = 0;                              // 最后一趟调节参数
    double P = 0;

    // error
    double WWithoutError = 0;                       // 未发生错误时的总任务量
    int numberWithoutError = 0;                     // 未发生错误时的处理机个数
    int beforeInstallment = 0;                      // 未发生错误时的趟数
    double leftW = 0.0;                             // 发生错误剩下的任务量
    double startTime = 0.0;                         // 发生错误时的开始时间

    // alpha && beta
    vector<double> alpha;                           // 内部调度中每趟各处理器的负载分量
    vector<double> beta;                            // 最后一趟调度中各处理器的负载分量

    // time
    vector<double> usingTimes;                      // 每个处理器的处理时间总和

    // print
    string outputName;

    // alpha && beta function
    // void initAlpha();
    // void initBeta();
};


#endif //FIRSTMODEL_APMISRR_H