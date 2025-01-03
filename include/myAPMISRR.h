/*
 * @FilePath: \InstallmentScheduling\include\myAPMISRR.h
 * @Description:  APMISRR add cost, non-block, remove P0
 * @Author: rthete
 * @Date: 2023-05-12 15:51:41
 * @LastEditTime: 2023-09-17 14:36:57
 */
#ifndef _myAPMISRR_H
#define _myAPMISRR_H

#include "Server.h"
#include "header.h"

class myAPMISRR {
public:
  myAPMISRR(int valueN, double valueTheta);

  /**
   * assign model
   */
  void setM(int value);
  void setW(double value);
  void setLambda(double value);

  double getAlpha();
  double getBeta();
  void initValue();
  void getDataFromFile(string data_path = "../data/15-servers-w-20/");

  // result
  void calOptimalTime();
  void calUsingRate();
  double getOptimalTime();
  double getUsingRate();

  // check schedulable
  int isSchedulable();

  void SISinitValue();
  void error(vector<int> &errorPlace, int errorInstallment);
  void addServer(int installment, vector<Server> &new_servers);

private:
  int n;                  // 处理机个数
  int m;                  // 计算趟数
  double theta;           // 结果回传比例
  double W = 0;           // 总任务量
  double V = 0;           // 内部调度每一趟调度的任务量
  double Vb = 0;          // 最后一趟调度的任务量
  double optimalTime = 0; // 调度最短时间
  double usingRate;       // 处理机使用率
  vector<Server> servers; // 所有的处理器

  double lambda = 0; // 最后一趟调节参数
  double P = 0;

  // error
  double WWithoutError = 0;   // 未发生错误时的总任务量
  int numberWithoutError = 0; // 未发生错误时的处理机个数
  int beforeInstallment = 0;  // 未发生错误时的趟数
  double leftW = 0.0;         // 发生错误剩下的任务量
  double startTime = 0.0;     // 发生错误时的开始时间

  // alpha && beta
  vector<double> alpha; // 内部调度中每趟各处理器的负载分量
  vector<double> beta;  // 最后一趟调度中各处理器的负载分量

  // time
  vector<double> usingTime; // 每个处理器的处理时间总和

  // print
  string outputName;

  // alpha && beta function
  // void initAlpha();
  // void initBeta();
};

#endif //_myAPMISRR_H