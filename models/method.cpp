/*
 * @FilePath: \InstallmentScheduling\models\method.cpp
 * @Description:
 * @Author: rthete
 * @Date: 2023-05-15 15:51:07
 * @LastEditTime: 2024-01-25 14:25:56
 */

#include "method.h"

namespace ModelRunner {
bool output_using_rate = false; // 控制是否输出ur
/**
 * 运行SIS模型
 */
double run_SIS(int server_num, int workload, double theta, string data_path,
               vector<int> error_server) {
  auto serverN = server_num; // number of servers

  cout << "**********************run SIS**********************" << endl;
  cout << serverN << "servers, theta = " << theta
       << ", total load = " << workload << endl;

  SIS sis(serverN, theta);
  sis.getDataFromFile(data_path);
  sis.setW((double)workload);
  sis.initValue();
  sis.getOptimalModel();
  if (!error_server.empty()) {
    cout << "error server size: " << to_string(error_server.size()) << endl;
    sis.error(error_server);
  }
  sis.calUsingRate();
  cout << "sis.getUsingRate(): " << sis.getUsingRate() << endl;
  cout << "sis.getOptimalTime(): " << sis.getOptimalTime() << endl;
  if (output_using_rate == 1)
    return sis.getUsingRate();
  return sis.getOptimalTime();
}

void run_PMIS() {
  auto workload = 1000; // total workload
  auto serverN = 15;    // number of servers
  auto theta = 0.3;     // Ratio of the output load size to input load size
  auto m = 30;          // installment size

  auto fPMIS = fopen("../output/PMIS.txt", "w");
  fprintf(fPMIS, "----------m: %d----------\n", m);
  fprintf(fPMIS, "workload\tPMIS\n");

  PMIS pmis(serverN, theta);
  pmis.getDataFromFile();
  pmis.setW((double)workload);
  pmis.setM(m);
  pmis.initValue();
  pmis.getOptimalModel();

  fprintf(fPMIS, "%d\t\t%.2lf\n", workload, pmis.getOptimalTime());
  pmis.printResult();

  cout << "done" << endl;
}

/**
 * 运行MISRR模型
 */
double run_MISRR(int server_num, int m, int workload, double theta,
                 string data_path, vector<int> error_server,
                 int error_installment, vector<Server> add_servers,
                 int add_installment) {
  // auto workload = 8000;   // total workload
  // auto theta = 0.3;       // Ratio of the output load size to input load size
  // auto m = 8;            // installment size

  cout << "**********************run MISRR**********************" << endl;
  cout << server_num << " servers, m = " << m << ", theta = " << theta << endl;
  cout << "total load = " << workload << endl;

  MISRR misrr(server_num, theta, m);
  misrr.getDataFromFile(data_path);
  misrr.setW((double)workload);
  misrr.initValue();
  misrr.getOptimalModel();

  cout << "misrr.getUsingRate(): " << misrr.getUsingRate() << endl;
  cout << "misrr.getOptimalTime(): " << misrr.getOptimalTime() << endl;

  // 无故障/恢复，直接返回
  if (error_server.empty() && add_servers.empty()) {
    if (output_using_rate == 1)
      return misrr.getUsingRate();
    return misrr.getOptimalTime();
  }

  // 有故障，重新计算
  if (!error_server.empty()) {
    cout << "error server size: " << to_string(error_server.size())
         << ", error installment: " << error_installment << endl;
    misrr.error_2(error_server, error_installment);
  }

  // 有恢复，重新计算
  if (!add_servers.empty()) {
    cout << "add servers size: " << to_string(add_servers.size())
         << ", add installment: " << add_installment << endl;
    // 打印新增处理机信息(double)
    for (auto &server : add_servers) {
      cout << "server ogsw: " << server.getO() << " " << server.getG() << " "
           << server.getS() << " " << server.getW() << endl;
    }
    if (add_installment >= misrr.getOptimalM() - 2) {
      add_installment = misrr.getOptimalM() - 2;
      cout << "new add installment: " << add_installment << endl;
    }
    misrr.addServer(add_installment, add_servers);
  }

  if (misrr.isSchedulable == 0) {
    cout << "not schedulable" << endl;
    return 0;
  }

  cout << "new misrr.getUsingRate(): " << misrr.getUsingRate() << endl;
  cout << "new misrr.getOptimalTime(): " << misrr.getOptimalTime() << endl;

  if (output_using_rate == 1)
    return misrr.getUsingRate();
  return misrr.getOptimalTime();
}

void run_MISRR_conflict(vector<double> &waiting_time, int server_num, int m,
                        int workload, double theta, string data_path,
                        vector<int> error_server, int error_installment) {
  // auto workload = 8000;   // total workload
  // auto theta = 0.3;       // Ratio of the output load size to input load size
  // auto m = 8;            // installment size

  cout << "**********************run MISRR**********************" << endl;
  cout << server_num << " servers, m = " << m << ", theta = " << theta << endl;
  cout << "total load = " << workload << endl;

  MISRR misrr(server_num, theta, m);
  misrr.getDataFromFile(data_path);
  misrr.setW((double)workload);
  misrr.initValue();
  misrr.getOptimalModel();

  // auto error_installment = 10;
  if (!error_server.empty()) {
    cout << "error server size: " << to_string(error_server.size())
         << ", error installment: " << error_installment << endl;
    misrr.error_2(error_server, error_installment);
  }

  cout << "misrr.getUsingRate(): " << misrr.getUsingRate() << endl;
  cout << "misrr.getOptimalTime(): " << misrr.getOptimalTime() << endl;
  misrr.getWaitingTime(waiting_time);
}

/**
 * 运行不带启动开销的APMISRR算法(original)
 */
double run_APMISRR(double lambda) {
  auto serverN = 14; // number of servers
  auto theta = 0.3;  // Ratio of the output load size to input load size
  auto m = 9;        // installment size

  cout << "**********************run APMISRR**********************" << endl;
  cout << serverN + 1 << " servers, m = " << m << ", theta = " << theta
       << ", lambda = " << lambda << endl;
  cout << "last installment load = " << lambda / m << endl;
  cout << "each internal installment load = " << (m - lambda) / (m * (m - 1))
       << endl;

  APMISRR apmisrr(serverN, theta);
  apmisrr.getDataFromFile();
  apmisrr.setM((int)m);
  apmisrr.setLambda((double)lambda);
  apmisrr.initValue();
  apmisrr.isSchedulable();
  return apmisrr.getOptimalTime();
}

/**
 * 运行带启动开销的APMISRR算法
 */
double run_APMISRR_cost(double lambda, int m) {
  auto serverN = 14; // number of servers
  auto theta = 0.3;  // Ratio of the output load size to input load size
  // auto m = 8;            // installment size

  cout << "**********************run APMISRR**********************" << endl;
  cout << serverN + 1 << " servers, m = " << m << ", theta = " << theta
       << ", lambda = " << lambda << endl;
  cout << "total load = " << 8000 << endl;
  cout << "last installment load = " << lambda / m * 8000 << endl;
  cout << "each internal installment load = "
       << (m - lambda) * 8000 / (m * (m - 1)) << endl;

  APMISRR apmisrr(serverN, theta);
  apmisrr.getDataFromFile();
  apmisrr.setM((int)m);
  apmisrr.setLambda((double)lambda);
  apmisrr.initValue_cost();
  apmisrr.isSchedulable_cost();
  return apmisrr.getOptimalTime_cost(); // test APMISRR total time
  // return apmisrr.getAlpha(); // test APMISRR alpha_1
  // return apmisrr.getBeta(); // test APMISRR beta_1
  return 0;
}

/**
 * 运行带启动开销，非阻塞，去掉P0的APMISRR算法
 */
double run_myAPMISRR(int server_num, double lambda, int m, int workload,
                     double theta, string data_path, vector<int> error_server,
                     int error_installment, vector<Server> add_servers,
                     int add_installment) {
  auto serverN = server_num;

  cout << "**********************run myAPMISRR**********************" << endl;
  cout << serverN << " servers, m = " << m << ", theta = " << theta
       << ", lambda = " << lambda << endl;
  cout << "total load = " << workload << endl;
  cout << "last installment load = " << lambda / m * workload << endl;
  cout << "each internal installment load = "
       << (m - lambda) * workload / (m * (m - 1)) << endl;

  myAPMISRR myapmisrr(serverN, theta);
  myapmisrr.getDataFromFile(data_path);
  myapmisrr.setW(workload);
  myapmisrr.setM((int)m);
  myapmisrr.setLambda((double)lambda);
  myapmisrr.initValue();
  if (myapmisrr.isSchedulable() != 1)
    return 0;
  myapmisrr.calOptimalTime();
  myapmisrr.calUsingRate();

  cout << "myapmisrr.getUsingRate(): " << myapmisrr.getUsingRate() << endl;
  cout << "myapmisrr.getOptimalTime(): " << myapmisrr.getOptimalTime() << endl;

  // 无故障/恢复，直接返回
  if (error_server.empty() && add_servers.empty()) {
    if (output_using_rate == 1)
      return myapmisrr.getUsingRate();
    return myapmisrr.getOptimalTime();
  }

  // 有故障，重新计算
  if (!error_server.empty()) {
    cout << "error server size: " << to_string(error_server.size())
         << ", error installment: " << error_installment << endl;
    myapmisrr.error(error_server, error_installment);
  }

  // 有恢复，重新计算
  if (!add_servers.empty()) {
    cout << "add servers size: " << to_string(add_servers.size())
         << ", add installment: " << add_installment << endl;
    // 限制add_installment的大小
    if (m - add_installment < 2) {
      add_installment = m - 2;
      cout << "new add installment: " << add_installment << endl;
    }
    // 打印新增处理机信息(double)
    for (auto &server : add_servers) {
      cout << "server ogsw: " << server.getO() << " " << server.getG() << " "
           << server.getS() << " " << server.getW() << endl;
    }
    myapmisrr.addServer(add_installment, add_servers);
  }

  cout << "new myapmisrr.getUsingRate(): " << myapmisrr.getUsingRate() << endl;
  cout << "new myapmisrr.getOptimalTime(): " << myapmisrr.getOptimalTime()
       << endl;

  if (output_using_rate == 1)
    return myapmisrr.getUsingRate();
  return myapmisrr.getOptimalTime();
}

/**
 * 运行最后一趟含lambda参数的MISRRL算法
 */
double run_MISRRL(double lambda, int m) {
  auto serverN = 15;
  auto theta = 0.3;
  auto workload = 8000;

  cout << "**********************run MISRRL**********************" << endl;
  cout << serverN << " servers, m = " << m << ", theta = " << theta
       << ", lambda = " << lambda << endl;
  cout << "total load = " << 8000 << endl;
  cout << "last installment load = " << lambda / m * 8000 << endl;
  cout << "each internal installment load = "
       << (m - lambda) * 8000 / (m * (m - 1)) << endl;

  MISRRL misrrl(serverN, theta, m);
  misrrl.getDataFromFile();
  misrrl.setW((double)workload);
  misrrl.setLambda((double)lambda);
  misrrl.initValue();
  misrrl.getOptimalModel();
  // misrrl.setM((int)m);
  // misrrl.isSchedulable();
  misrrl.theLastInstallmentGap(to_string(lambda));
  if (misrrl.isSchedulable == 0) {
    cout << "not schedulable" << endl;
    return 0;
  }
  cout << "misrrl.getUsingRate(): " << misrrl.getUsingRate() << endl;
  cout << "misrrl.getOptimalTime(): " << misrrl.getOptimalTime() << endl;
  // return misrrl.getUsingRate();
  return misrrl.getOptimalTime();
}

/**
 * 运行第一趟、最后一趟含lambda参数的MISRRLL算法
 */
double run_MISRRLL(int server_num, double lambda1, double lambda2, int m,
                   int workload, double theta, string data_path,
                   vector<int> error_server, int error_installment,
                   vector<Server> add_servers, int add_installment) {
  cout << "**********************run MISRRLL**********************" << endl;
  auto serverN = server_num;

  MISRRLL misrrll(serverN, theta, m);
  misrrll.getDataFromFile(data_path);
  misrrll.setW((double)workload);
  misrrll.setLambda((double)lambda1, (double)lambda2);
  misrrll.initValue();
  m = misrrll.getOptimalM();

  cout << serverN << " servers, m = " << m << ", theta = " << theta
       << ", lambda1 = " << lambda1 << ", lambda2 = " << lambda2 << endl;
  cout << "total load = " << workload << endl;
  cout << "last installment load = " << lambda1 / m * workload << endl;
  cout << "first installment load = " << lambda2 / m * workload << endl;
  cout << "each internal installment load = "
       << (m - lambda1 - lambda2) * workload / (m * (m - 2)) << endl;

  misrrll.getOptimalModel();

  misrrll.theLastInstallmentGap(to_string(lambda1));
  if (misrrll.isSchedulable == 0) {
    cout << "not schedulable" << endl;
    return 0;
  }
  cout << "misrrll.getUsingRate(): " << misrrll.getUsingRate() << endl;
  cout << "misrrll.getOptimalTime(): " << misrrll.getOptimalTime() << endl;

  // 无故障/恢复，直接返回
  if (error_server.empty() && add_servers.empty()) {
    if (output_using_rate == 1)
      return misrrll.getUsingRate();
    return misrrll.getOptimalTime();
  }

  // 有故障，重新计算
  if (!error_server.empty()) {
    cout << "error server size: " << to_string(error_server.size())
         << ", error installment: " << error_installment << endl;
    if (error_installment >= misrrll.getOptimalM() - 2) {
      error_installment = misrrll.getOptimalM() - 2;
      cout << "new error installment: " << error_installment << endl;
    }
    misrrll.error(error_server, error_installment);
    // cout << "error optimal m = " << misrrll.getOptimalM() << endl;
  }

  // 有恢复，重新计算
  if (!add_servers.empty()) {
    cout << "add servers size: " << to_string(add_servers.size())
         << ", add installment: " << add_installment << endl;
    // 打印新增处理机信息(double)
    for (auto &server : add_servers) {
      cout << "server ogsw: " << server.getO() << " " << server.getG() << " "
           << server.getS() << " " << server.getW() << endl;
    }
    if (add_installment >= misrrll.getOptimalM() - 2) {
      add_installment = misrrll.getOptimalM() - 2;
      cout << "new add installment: " << add_installment << endl;
    }
    misrrll.addServer(add_installment, add_servers);
  }

  cout << "new misrrll.getUsingRate(): " << misrrll.getUsingRate() << endl;
  cout << "new misrrll.getOptimalTime(): " << misrrll.getOptimalTime() << endl;

  if (misrrll.isSchedulable == 0) {
    cout << "not schedulable" << endl;
    return 0;
  }

  if (output_using_rate == 1)
    return misrrll.getUsingRate();
  return misrrll.getOptimalTime();
}

} // namespace ModelRunner

/**
 * MISRR: 测试theta对可行性的影响
 */
void test_MISRR_theta() {
  // FILE * fpResult;
  // fpResult = fopen("../output/test_MISRR_theta_calM.csv", "w");
  for (int i = 0; i < 20; i += 1) {
    ModelRunner::run_MISRR(i, 8);
    // cout << "theta = " << i;
    // cout<< get<0>(ModelRunner::run_MISRR(i, 8)) << endl;
    // cout << "\ttime = " << get<0>(ModelRunner::run_MISRR(i, 0)) << "\tusing
    // rate = " << get<1>(ModelRunner::run_MISRR(i, 0)) << endl;
    // fprintf(fpResult, "%d,%lf,%lf\n", i, get<0>(ModelRunner::run_MISRR(i,
    // 0)), get<1>(ModelRunner::run_MISRR(i, 0)));
  }
  // fclose(fpResult);
}

/**
 * MISRR:测试趟数m的影响
 */
void test_MISRR_all() {
  FILE *fpResult;
  fpResult = fopen("../output/MISRR/test_MISRR_time.csv", "w");
  for (int m = 3; m <= 40; m++) {
    fprintf(fpResult, "%lf\n", ModelRunner::run_MISRR(m, 8000, 0.3));
  }
  fclose(fpResult);
}

/**
 * 测试total time与\lambda的关系
 */
void test_APMISRR_totalTime() {
  double lambda[10];
  double time[10];
  for (int i = 0; i < 11; i++) {
    lambda[i] = 0.1 * i;
    time[i] = ModelRunner::run_APMISRR(lambda[i]);
  }

  // 实验：总时间与\lambda线性正相关
  cout << "***********"
       << "(time[i] - time[i-1]) / 0.1"
       << "***********" << endl;
  for (int i = 1; i < 10; i++) {
    cout << (time[i] - time[i - 1]) / 0.1 << endl;
  }
  cout << "Total makespan is linearly dependant on \\lambda." << endl;
}

/**
 * 测试alpha[1]与\lambda的关系
 */
void test_APMISRR_alpha() {
  double lambda[10];
  double time[10];
  double alpha1[10];
  for (int i = 0; i < 11; i++) {
    lambda[i] = 0.1 * i;
    alpha1[i] = ModelRunner::run_APMISRR_cost(lambda[i], 8);
  }

  // 实验：alpha[1]与\lambda负线性相关
  cout << "***********"
       << "(alpha1[i] - alpha1[i-1]) / 0.1"
       << "***********" << endl;
  for (int i = 1; i < 10; i++) {
    cout << (alpha1[i] - alpha1[i - 1]) / 0.1 << endl;
  }
  // cout << "alpha_1 is linearly dependant on \\lambda." << endl;
}

/**
 * 测试beta[1]与\lambda的关系
 */
void test_APMISRR_beta() {
  double lambda[10];
  double time[10];
  double beta1[10];
  for (int i = 0; i < 11; i++) {
    lambda[i] = 0.1 * i;
    beta1[i] = ModelRunner::run_APMISRR_cost(lambda[i], 8);
  }

  // 实验：beta[1]与\lambda负线性相关
  cout << "***********"
       << "(beta1[i] - beta1[i-1]) / 0.1"
       << "***********" << endl;
  for (int i = 1; i <= 10; i++) {
    cout << "lambda = " << lambda[i] << ", beta_1 = " << beta1[i] << endl;
    cout << (beta1[i] - beta1[i - 1]) / 0.1 << endl;
  }
  // cout << "beta_1 is linearly dependant on \\lambda." << endl;
}

/**
 * 测试total time与lambda的关系
 */
void test_APMISRR_totalTime_cost() {
  double lambda[10];
  double time[10];
  for (int i = 0; i < 11; i++) {
    lambda[i] = 0.1 * i;
    time[i] = ModelRunner::run_APMISRR_cost(lambda[i], 8);
  }

  // 实验：总时间与\lambda线性正相关
  cout << "***********"
       << "(time[i] - time[i-1]) / 0.1"
       << "***********" << endl;
  for (int i = 1; i < 10; i++) {
    cout << (time[i] - time[i - 1]) / 0.1 << endl;
  }
  // cout << "Total makespan is linearly dependant on \\lambda." << endl;
}

/**
 * APMISRR: 测试总时间与趟数m的关系
 */
void test_APMISRR_installment() {
  FILE *fpResult;
  fpResult = fopen("../output/test.csv", "w");
  for (int m = 3; m <= 20; m++) {
    fprintf(fpResult, "%lf\n", ModelRunner::run_APMISRR_cost(0.2, m));
  }
  fclose(fpResult);
}

/**
 * myAPMISRR: 测试总时间与趟数m的关系
 * lambda < 0.2 || > 0.8不可行
 */
void test_myAPMISRR_installment() {
  FILE *fpResult;
  fpResult = fopen("../output/myAPMISRR/test_myAPMISRR_time.csv", "w");
  for (int m = 3; m <= 40; m++) {
    for (double lambda = 0.1; lambda < 1; lambda += 0.1) {
      fprintf(fpResult, "%lf,",
              ModelRunner::run_myAPMISRR(lambda, m, 8000, 0.3));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

/**
 * MISRRL: 测试总时间与lambda的关系
 * 结论：lambda越小，用时越短
 */
void test_MISRRL_lambda() {
  for (double lambda = 0; lambda < 1; lambda += 0.1) {
    ModelRunner::run_MISRRL(lambda, 8);
    // ModelRunner::run_MISRRL(lambda, 20);
  }
}

/**
 * MISRRL: 测试总时间与m的关系
 * 结论: 有最优趟数m
 */
void test_MISRRL_installment() {
  for (int m = 3; m < 40; m++) {
    ModelRunner::run_MISRRL(0.2, m);
  }
}

void test_MISRRL_all() {
  FILE *fpResult;
  fpResult = fopen("../output/MISRRL/test_MISRRL_time.csv", "w");
  for (int m = 3; m <= 40; m++) {
    for (double lambda = 0.1; lambda < 1; lambda += 0.1) {
      fprintf(fpResult, "%lf,", ModelRunner::run_MISRRL(lambda, m));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void test_MISRRLL_lambda2() {
  for (double lambda = 0; lambda < 5; lambda += 0.5) {
    // ModelRunner::run_MISRRLL(0.6, lambda, 8, 8000);
  }
}

void test_MISRRLL_lambda1() {
  for (double lambda = 0; lambda < 1; lambda += 0.1) {
    // ModelRunner::run_MISRRLL(lambda, 0.6, 8, 8000);
  }
}

/**
 * MISRRL: 测试总时间与m的关系
 */
void test_MISRRLL_installment() {
  vector<double> time_list{};
  for (int m = 3; m < 40; m++) {
    // time_list.push_back(ModelRunner::run_MISRRLL(0.6, 0.6, m, 8000));
  }

  for (auto i = 0; i < time_list.size(); i++) {
    cout << time_list[i] << endl;
  }
}

void test_MISRRLL_all() {
  FILE *fpResult;
  fpResult = fopen("../output/MISRRLL/test_MISRRL_time_lambda2.csv", "w");
  for (int m = 3; m <= 40; m++) {
    for (double lambda = 0.1; lambda < 1; lambda += 0.1) {
      // fprintf(fpResult, "%lf,", ModelRunner::run_MISRRLL(0.2, lambda, m,
      // 8000));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void test_MISRRLL_2_lambda() {
  FILE *fpResult;
  fpResult = fopen("../output/MISRRLL/test_MISRRL_time_2_lambda.csv", "w");
  for (double lambda1 = 0.1; lambda1 < 1; lambda1 += 0.1) {
    for (double lambda2 = 0.1; lambda2 < 1; lambda2 += 0.1) {
      // fprintf(fpResult, "%lf,",
      //         ModelRunner::run_MISRRLL(lambda1, lambda2, 15, 8000));
    }
    fprintf(fpResult, "\n");
  }
  fclose(fpResult);
}

void compare_MISRR_and_MISRRL() {
  for (int m = 3; m < 20; m++) {
    ModelRunner::run_MISRR(0.3, m);
    ModelRunner::run_MISRRL(0.2, m);
  }
}