/*
 * @FilePath: \InstallmentScheduling\models\myAPMISRR.cpp
 * @Description: APMISRR add cost, non-block, remove P0
 * @Author: rthete
 * @Date: 2023-05-12 15:55:34
 * @LastEditTime: 2023-09-17 15:04:01
 */

#include "myAPMISRR.h"

myAPMISRR::myAPMISRR(int valueN, double valueTheta)
    : n(valueN), theta(valueTheta), usingRate(0),
      servers(vector<Server>(MAX_N)),

      // error
      numberWithoutError(valueN),

      // alpha & beta
      alpha(vector<double>(MAX_N, 0)), beta(vector<double>(MAX_N, 0)),

      // time
      usingTime(vector<double>(MAX_N, 0)) {}

// init value
void myAPMISRR::setW(double value) { this->W = value; }

void myAPMISRR::setM(int value) { this->m = value; }

void myAPMISRR::setLambda(double value) { this->lambda = value; }

double myAPMISRR::getAlpha() { return alpha[0]; }

double myAPMISRR::getBeta() { return beta[0]; }

void myAPMISRR::initValue() {
  // 每趟调度任务量
  this->V = (this->m - this->lambda) * this->W / (m * (m - 1));
  cout << "V = " << V << endl;

  // 式(11)
  double temp = 0;
  for (int i = 0; i < this->n; ++i) {
    temp += 1 / servers[i].getW();
  }
  double temp1 = 0;
  for (int i = 0; i < this->n; ++i) {
    temp1 += (servers[i].getS() / servers[i].getW());
  }
  this->P = (this->V + temp1) / temp;
  cout << "P = " << P << endl;

  // 式(12)
  double sum = 0;
  for (int i = 0; i < this->n; ++i) {
    alpha[i] = (P - servers[i].getS()) / (servers[i].getW() * this->V);
    // cout << "alpha[" << i << "] = " << alpha[i] << endl;
    sum += alpha[i];
    // cout << "start cost for computation: " << servers[i].getS() << endl;
    // cout << "compute time in first installment: " << alpha[i] *
    // servers[i].getW() * this->V << endl; cout << "start cost + compute time =
    // " << servers[i].getS() + alpha[i] * servers[i].getW() * this->V << endl;
    // // ==P
  }
  // cout << "sum = " << sum << endl; // alpha[i]之和==1

  temp = 0;
  for (int i = 1; i < this->n; ++i) {
    temp += servers[i].getO() + alpha[i] * servers[i].getG() * this->V +
            servers[i].getO() + alpha[i] * theta * servers[i].getG() * this->V;
  }

  // 内部调度间空闲时间
  double d = P - temp;
  cout << "d = " << d << endl;

  this->Vb = lambda / m * this->W;
  // 式(15)
  vector<double> a(this->n, 0);
  for (int i = 0; i < this->n; ++i) {
    a[i] =
        (servers[i].getW() + theta * servers[i].getG()) / servers[i + 1].getW();
  }
  vector<double> c(this->n, 0);
  for (int i = 0; i < this->n; ++i) {
    c[i] = (-alpha[i] * (1 + theta) * servers[i].getG() * this->V +
            servers[i].getS() - servers[i].getO() - servers[i + 1].getS() -
            servers[i + 1].getO()) /
           (servers[i + 1].getW() * this->Vb);
  }

  // 式(17)
  vector<double> D(this->n, 0);
  for (int i = 1; i < this->n; ++i) {
    for (int j = 0; j < i - 1; ++j) {
      double multi_a = 1;
      for (int k = j + 1; k < i - 1; ++k) {
        multi_a *= a[k];
      }
      D[i] += c[j] * multi_a;
    }
  }

  // 式(20)
  double A = 1;
  for (int i = 1; i < this->n; ++i) {
    double multi_a = 1;
    for (int j = 0; j < i - 1; ++j) {
      multi_a *= a[j];
    }
    A += multi_a;
  }

  // 式(21)
  double B = 0;
  for (int i = 1; i < this->n; ++i) {
    B += D[i];
  }

  // 式(19)
  beta[0] = (1 - B) / A;

  // 式(18)
  sum = beta[0];
  for (int i = 1; i < this->n; ++i) {
    double multi_a = 1;
    for (int j = 0; j < i - 1; ++j) {
      multi_a *= a[j];
    }
    beta[i] = multi_a * beta[0] + D[i];
    // cout << "beta[i] = " << beta[i] << endl;
    // sum += beta[i];
  }
  // cout << sum << endl;
}

void myAPMISRR::calOptimalTime() {
  // 式(29)
  this->optimalTime = (m - 1) * P + servers[0].getO() + servers[0].getS() +
                      beta[0] * servers[0].getW() * this->Vb;
  for (int i = 0; i < this->n; ++i) {
    this->optimalTime +=
        servers[i].getO() + beta[i] * servers[i].getG() * theta * this->Vb;
  }
  this->optimalTime -= servers[0].getO();
  // cout << "optimalTime = " << optimalTime << endl;

  this->optimalTime += this->startTime;
}

void myAPMISRR::calUsingRate() {

  for (int i = 0; i < this->n; i++) {
    int id = servers[i].getId();
    usingTime[id] +=
        ((m - 1) * (alpha[i] * V * servers[i].getW() + servers[i].getS()) +
         beta[i] * Vb * servers[i].getW() + servers[i].getS() +
         beta[i] * Vb * servers[i].getG() * theta);
  }

  // cal using rate
  usingRate = 0.0;
  for (int i = 0; i < this->numberWithoutError; ++i) {
    // cout << "usingTime[" << i << "] = " << usingTime[i] << endl;
    usingRate += ((double)usingTime[i] /
                  ((double)(this->optimalTime) * this->numberWithoutError));
  }
}

int myAPMISRR::isSchedulable() {

  // 打印alpha[i], beta[i]
  double load_sum = 0;
  for (int i = 0; i < this->n; ++i) {
    // cout << "alpha" << i << "=" << alpha[i] << "\t";
    // cout << "beta" << i << "=" << beta[i] << endl;
    load_sum += alpha[i] * (m - 1) * this->V;
    load_sum += beta[i] * this->Vb;
  }
  // load_sum为1
  // cout << "load_sum = " << load_sum << endl;

  int condition1 = 1, condition2 = 1, condition3 = 1;

  // 负载分量必须为正数
  for (int i = 0; i < this->n; ++i) {
    if (!(alpha[i] > 0 && beta[i] > 0)) {
      condition1 = 0;
      cout << "!!!!!!! not schedulable 1 !!!!!!!" << endl;
      break;
    }
  }

  // 引理8
  int num = this->n - 1;
  if (servers[num].getS() + beta[num] * servers[num].getW() * this->Vb <
      beta[num] * servers[num].getG() * this->Vb + +servers[num].getO() +
          alpha[num] * theta * servers[num].getG() * this->V +
          servers[num - 1].getO() +
          beta[num - 1] * theta * servers[num - 1].getG() * this->Vb) {
    condition2 = 0;
    cout << "!!!!!!! not schedulable 2 !!!!!!!" << endl;
  }

  // 处理器P_1的最后一趟负载分量是否能在传输窗口放下
  if (alpha[0] * this->V < beta[0] * this->Vb) {
    condition3 = 0;
    cout << "!!!!!!! not schedulable 3 !!!!!!!" << endl;
  }

  if (condition1 && condition2 && condition3) {
    cout << "```````````` schedulable ````````````" << endl;
    return 1;
  } else {
    return 0;
  }
}

void myAPMISRR::SISinitValue() {
  /**
   * cal eta & gamma (alpha)
   * SIS
   */
  vector<double> mu(this->n, 0);
  vector<double> lambda(this->n, 0);
  // cal mu & lambda
  for (int i = 1; i < this->n; ++i) {
    mu[i] = (servers[i - 1].getS() - servers[i].getO() - servers[i].getS()) /
            servers[i].getW();
    lambda[i] =
        (servers[i - 1].getW() - servers[i - 1].getG()) / servers[i].getW();
  }

  vector<double> eta(this->n, 0);
  vector<double> gamma(this->n, 0);
  // cal eta & gamma
  for (int i = 1; i < this->n; ++i) {
    eta[i] = 1.0;
    for (int j = 1; j <= i; ++j) {
      eta[i] *= lambda[j];
    }
  }
  for (int i = 1; i < this->n; ++i) {
    for (int j = 1; j <= i; j++) {
      double temp = 1.0;
      for (int k = j + 1; k <= i; k++) {
        temp *= lambda[k];
      }
      gamma[i] += (mu[j] * temp);
    }
  }

  // cal alpha
  double sum_2_n_gamma = 0, sum_2_n_eta = 0;
  for (int i = 1; i < this->n; ++i) {
    sum_2_n_eta += eta[i];
    sum_2_n_gamma += gamma[i];
  }
  for (int i = 0; i < this->n; ++i) {
    if (i == 0) {
      alpha[i] = (1.0 - sum_2_n_gamma / this->W) / (1.0 + sum_2_n_eta);
      continue;
    }
    alpha[i] = gamma[i] / this->W + eta[i] * alpha[0];
  }

  double count_alpha = 0;
  for (int i = 0; i < this->n; ++i) {
    count_alpha += alpha[i];
    if (alpha[i] < 0) {
      cout << this->W << " error "
           << "alpha - " << alpha[i] << endl;
    }
  }

  if ((int)((count_alpha + 0.000005) * 100000) != 100000)
    cout << this->W << " count: alpha - " << count_alpha << endl;
}

void myAPMISRR::error(vector<int> &errorPlace, int errorInstallment) {
  // cal left workload
  for (auto i : errorPlace) {
    this->leftW += ((this->m - errorInstallment) * alpha[i - 1] * this->V +
                    beta[i - 1] * this->Vb);
    usingTime[i - 1] = this->P * (errorInstallment - 1);
  }
  cout << "left workload: " << this->leftW << endl;
  this->W = this->leftW;

  // cal left server
  int position = 0;
  vector<Server> oldServers = servers;
  servers.clear();

  for (int i = 0; i < oldServers.size(); i++) {
    auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
    if (iter == errorPlace.end()) {
      servers[position] = oldServers[i];
      position++;
    }
  }
  int serversNumberWithoutError = this->n;
  this->n -= (int)errorPlace.size();

  // 最后加一趟单趟调度
  SISinitValue();

  this->optimalTime += servers[0].getO() + servers[0].getS() +
                       servers[0].getW() * alpha[0] * this->W;
  for (int i = 0; i < this->n; ++i) {
    this->optimalTime += (servers[i].getO() +
                          servers[i].getG() * this->theta * alpha[i] * this->W);
  }

  // cal using time
  for (int i = 0; i < serversNumberWithoutError; ++i) {
    // cout << i << " " << usingTime[i] << endl;
    auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
    if (iter == errorPlace.end()) {
      // 计算正常处理机的使用时间
      usingTime[i] =
          usingTime[i] + (servers[i].getO() + servers[i].getS() +
                          servers[i].getW() * alpha[i] * this->W +
                          servers[i].getG() * alpha[i] * this->W * theta);
      // cout << i << " " << usingTime[i] << endl;
    }
  }

  // cal using rate
  usingRate = 0.0;
  for (int i = 0; i < serversNumberWithoutError; ++i) {
    usingRate += ((double)usingTime[i] /
                  ((double)(this->optimalTime) * serversNumberWithoutError));
    // cout << i << " " << ((double)usingTime[i] /
    // ((double)(this->optimalTime))) << endl;
  }
}

void myAPMISRR::addServer(int add_installment, vector<Server> &new_servers) {
  cout << "**********************start myAPMISRR recover**********************"
       << endl;

  // 重调度开始时间，在最后一个处理机回传结束后才开始新调度
  this->startTime = servers[0].getS() * (add_installment + 1) +
                    servers[0].getW() * alpha[0] * this->V +
                    servers[0].getW() * alpha[0] * add_installment * this->V +
                    servers[0].getG() * beta[0] * this->V * theta;
  cout << "startTime: " << startTime << endl;
  for (int i = 0; i < this->n; ++i) {
    this->startTime += 2 * servers[i].getO() +
                       servers[i].getG() * alpha[i] * this->V * (1 + theta);
  }
  cout << "startTime: " << startTime << endl;

  // 重置usingTime(把新处理机的usingTime加到最后了)
  usingTime.clear();
  usingTime.resize(this->n + new_servers.size(), 0);
  for (int i = 0; i < this->n; ++i) {
    usingTime[i] =
        (servers[i].getS() + alpha[i] * this->V * servers[i].getW()) *
        (add_installment + 1);
  }

  // 给新处理机设置编号
  for (int i = 0; i < new_servers.size(); i++) {
    new_servers[i].setID(this->n + i);
  }

  // 重置处理机
  int num_new_servers = new_servers.size();
  int original_size = this->n;
  servers.resize(original_size + num_new_servers);
  for (int i = num_new_servers; i < this->n + num_new_servers; i++) {
    servers[i].setO(servers[i - num_new_servers].getO());
    servers[i].setS(servers[i - num_new_servers].getS());
    servers[i].setG(servers[i - num_new_servers].getG());
    servers[i].setW(servers[i - num_new_servers].getW());
  }

  // 设置新处理机的参数
  for (int i = 0; i < num_new_servers; i++) {
    servers[i] = new_servers[i];
  }

  this->leftW = 0;
  this->leftW += (this->Vb + (this->m - add_installment - 2) * this->V);
  cout << "this->leftW=" << this->leftW << endl;
  this->W = this->leftW;
  this->m = this->m - add_installment - 1;
  // 若剩余趟数小于3，则设置为3，否则无法调度
  this->m = this->m < 3 ? 3 : this->m;
  this->n += num_new_servers;
  this->numberWithoutError = this->n;

  this->lambda = 1;

  initValue();
  calOptimalTime();
  calUsingRate();
  cout << "this->startTime=" << this->startTime << endl;
  cout << "recover new optimal time: " << this->optimalTime << endl;

  // 重新计算usingRate
  usingRate = 0.0;
  for (int i = 0; i < this->n; ++i) {
    // cout << "usingTime[" << i << "] = " << usingTime[i] << endl;
    usingRate += ((double)usingTime[i] /
                  ((double)(this->optimalTime) * (this->n - num_new_servers) +
                   (this->optimalTime - this->startTime) * num_new_servers));
  }

  cout << "**********************end myAPMISRR recover**********************"
       << endl;
}

double myAPMISRR::getOptimalTime() { return this->optimalTime; }

double myAPMISRR::getUsingRate() { return this->usingRate; }

/**
 * @brief Read o & s & g & w & WTotal from "/data" and set every server.
 *
 */
void myAPMISRR::getDataFromFile(string data_path) {
  FILE *fpo, *fps, *fpg, *fpw, *totalW;
  double valueO[this->n], valueS[this->n], valueG[this->n], valueW[this->n];
  ;

  fpo = fopen((data_path + "o.txt").c_str(), "r");
  fps = fopen((data_path + "s.txt").c_str(), "r");
  fpg = fopen((data_path + "g.txt").c_str(), "r");
  fpw = fopen((data_path + "w.txt").c_str(), "r");
  totalW = fopen((data_path + "WTotal.txt").c_str(), "r");

  if (fpo == nullptr || fps == nullptr || fpg == nullptr || fpw == nullptr ||
      totalW == nullptr) {
    printf("The file can not be opened:\n");
    exit(-1);
  }

  fscanf(totalW, "%lf", &this->W);
  for (int i = 0; fscanf(fpo, "%lf", &valueO[i]) != EOF; ++i)
    ;
  for (int i = 0; fscanf(fps, "%lf", &valueS[i]) != EOF; ++i)
    ;
  for (int i = 0; fscanf(fpg, "%lf", &valueG[i]) != EOF; ++i)
    ;
  for (int i = 0; fscanf(fpw, "%lf", &valueW[i]) != EOF; ++i)
    ;

  for (int i = 0; i < this->n; i++) {
    Server demo(i);

    demo.setO(valueO[i]);
    demo.setS(valueS[i]);
    demo.setG(valueG[i]);
    demo.setW(valueW[i]);

    servers[i] = demo;
  }

  fclose(fpo);
  fclose(fps);
  fclose(fpg);
  fclose(fpw);
  fclose(totalW);
}
