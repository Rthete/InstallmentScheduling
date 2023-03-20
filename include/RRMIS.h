//
// Created by xd_wxl on 2023/1/17.
//

#ifndef FIRSTMODEL_RRMIS_H
#define FIRSTMODEL_RRMIS_H

#include "header.h"
#include "Server.h"

class RRMIS {
public:
    RRMIS(int valueN, double valueTheta);

    /**
     * assign model
     */
    void initValue();
    void getDataFromFile();
    void getOptimalModel();
    void setW(double value);

    // error
    void error(vector<int> &errorPlace, int errorInstallment);

    // print result
    void printResult();

    // init value
    void setM(int value);
    int getBestM();
    int getM();

    // result
    double getOptimalTime() const;
    double getUsingRate() const;

    // other
    void theLastInstallmentGap();
private:
    int n;                                      // 处理机个数
    int m;                                      // 调度趟数
    double theta;                               // 结果回传比例
    double W = 0;                               // 总任务量
    double V = 0;                               // 每趟调度的任务量大小
    double usingRate;                           // 处理机使用率
    double optimalTime = 0;                     // 调度最短时间
    vector<Server> servers;                     // 所有的处理器

    // error
    double WWithoutError = 0;                   // 未发生错误时的总任务量
    int serversNumberWithoutError = 0;          // 未发生错误时的处理机个数
    int leftM = 0;                              // 如果有错误、增加，之前的趟数 | 所有处理机
    double leftW = 0.0;                         // 发生错误剩下的任务量
    double startTime = 0.0;                     // 发生错误时的开始时间

    // beta
    vector<double> beta;                        // 内部调度分配量
    vector<double> Mu;                          // 内部调度 mu
    vector<double> Eta;                         // 内部调度 eta

    // alpha
    vector<double> alpha;                       // 第一趟调度中分配量
    vector<double> Delta;                       // 第一趟调度 delta
    vector<double> Phi;                         // 第一趟调度 phi
    vector<double> Epsilon;                     // 第一趟调度 Epsilon

    // gamma
    vector<double> gamma;                       // 最后一次调度中分配量
    vector<double> Lambda;                      // 最后一次调度 Lambda
    vector<double> Psi;                         // 最后一次调度 Psi
    vector<double> P;                           // 最后一次调度 P

    double old_V;
    vector<double> old_beta;

    // time
    vector<double> usingTime;

    // print
    string outputName = "RRMIS";

    void calAlpha();
    void calBeta();
    void calGamma();
    void calOptimalM();
};

#endif //FIRSTMODEL_RRMIS_H
