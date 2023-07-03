/*
 * @FilePath: \InstallmentScheduling\include\MISRRL.h
 * @Description:  
 * @Author: rthete
 * @Date: 2023-03-13 22:37:48
 * @LastEditTime: 2023-05-18 13:48:54
 */

#ifndef _MISRRL_H
#define _MISRRL_H

#include "header.h"
#include "Server.h"

class MISRRL {
public:
    MISRRL(int valueN, double valueTheta);
    MISRRL(int valueN, double valueTheta, int installment);

    void initValue();
    void getDataFromFile();
    void getOptimalModel();
    void setW(double value);
    void error(vector<int> &errorPlace, int errorInstallment);
    void printResult();
    double getOptimalTime() const;
    double getUsingRate() const;
    int getOptimalM() const;
    void theLastInstallmentGap(string title);

    void setLambda(double value);
    int isSchedulable = 1;

private:
    int n;                                      // 处理机个数
    int m = 0;                                  // 调度趟数
    double theta;                               // 结果回传比例
    double W = 0;                               // 总任务量
    double V = 0;                               // 每趟内部调度的任务量大小
    double Vb = 0;                              // 最后一趟任务量
    double usingRate;                           // 处理机使用率
    double optimalTime = 0;                     // 调度最短时间
    vector<Server> servers;                     // 所有的处理器

    double l = 0;

    // error
    double WWithoutError = 0;                   // 未发生错误时的总任务量
    int serversNumberWithoutError = 0;          // 未发生错误时的处理机个数
    int beforeInstallment = 0;                  // 未发生错误时的趟数
    double timeGap;                             // 冲突时间
    double leftW = 0.0;                         // 发生错误剩下的任务量
    double startTime = 0.0;                     // 发生错误时的开始时间

    // beta
    vector<double> beta;                        // 内部调度分配量
    vector<double> mu;                          // 内部调度 mu
    vector<double> eta;                         // 内部调度 eta

    // alpha
    vector<double> alpha;                       // 第一趟调度中分配量
    vector<double> Delta;                       // 第一趟调度 Delta
    vector<double> Phi;                         // 第一趟调度 Phi
    vector<double> Epsilon;                     // 第一趟调度 Epsilon

    // gamma
    vector<double> gamma;                       // 最后一次调度中分配量
    vector<double> Lambda;                      // 最后一次调度 Lambda
    vector<double> Psi;                         // 最后一次调度 Psi
    vector<double> P;                           // 最后一次调度 Rho

    // time gap
    double old_V;
    vector<double> old_beta;

    // time
    vector<double> usingTime;

    // print
    string outputName;

    void calAlpha();
    void calBeta();
    void calGamma();
    void calOptimalM();

};


#endif //_MISRRL_H
