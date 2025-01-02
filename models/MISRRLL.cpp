/*
 * @FilePath: \InstallmentScheduling\models\MISRRLL.cpp
 * @Description:
 * @Author: rthete
 * @Date: 2023-05-24 13:45:33
 * @LastEditTime: 2024-01-08 17:39:14
 */

#include "MISRRLL.h"
#include "Eigen/src/Core/util/Constants.h"
#include "unsupported/Eigen/src/Polynomials/PolynomialSolver.h"
#include <complex>

MISRRLL::MISRRLL(int valueN, double valueTheta, int installment)
    : n(valueN), m(installment), theta(valueTheta), usingRate(0),
      serversNumberWithoutError(valueN), servers(vector<Server>(MAX_N)),

      // beta
      beta(vector<double>(MAX_N, 0)), mu(vector<double>(MAX_N, 0)),
      eta(vector<double>(MAX_N, 0)),

      // alpha
      alpha(vector<double>(MAX_N, 0)), Delta(vector<double>(MAX_N, 0)),
      Phi(vector<double>(MAX_N, 0)), Epsilon(vector<double>(MAX_N, 0)),

      // gamma
      gamma(vector<double>(MAX_N, 0)), Lambda(vector<double>(MAX_N, 0)),
      Psi(vector<double>(MAX_N, 0)), P(vector<double>(MAX_N, 0)),

      // time gap
      old_beta(vector<double>(MAX_N)),

      // time
      usingTime(vector<double>(MAX_N, 0)), outputName("MISRRLL_print_result") {}

MISRRLL::MISRRLL(int valueN, double valueTheta)
    : n(valueN), theta(valueTheta), usingRate(0),
      serversNumberWithoutError(valueN), servers(vector<Server>(MAX_N)),

      // beta
      beta(vector<double>(MAX_N, 0)), mu(vector<double>(MAX_N, 0)),
      eta(vector<double>(MAX_N, 0)),

      // alpha
      alpha(vector<double>(MAX_N, 0)), Delta(vector<double>(MAX_N, 0)),
      Phi(vector<double>(MAX_N, 0)), Epsilon(vector<double>(MAX_N, 0)),

      // gamma
      gamma(vector<double>(MAX_N, 0)), Lambda(vector<double>(MAX_N, 0)),
      Psi(vector<double>(MAX_N, 0)), P(vector<double>(MAX_N, 0)),

      // time gap
      old_beta(vector<double>(MAX_N)),

      // time
      usingTime(vector<double>(MAX_N, 0)), outputName("MIS-CRR") {}

/**
 * @brief Assign model.
 *
 */
void MISRRLL::initValue() {
  if (m == 0) {
    calOptimalM();
  }

  // cal load for each installment
  // V: 内部调度的每趟任务量
  // cout << "1-this->V" << this->V << endl;
  // cout << "1-this->m=" << this->m << endl;
  // cout << "1-this->W=" << this->W << endl;
  this->V =
      (this->m - this->l_2 - this->l_1) * this->W / (this->m * (this->m - 2));
  // cout << "2-this->V" << this->V << endl;
  // Vb: 最后一趟的任务量
  this->Vb = this->l_2 / this->m * this->W;
  // Va: 第一趟的任务量
  this->Va = this->l_1 / this->m * this->W;

  // mu & eta: 内部调度参数
  // cal mu (beta)
  for (int i = 0; i < this->n; ++i) {
    mu[i] = servers[0].getW() / servers[i].getW();
  }

  // cal eta (beta)
  for (int i = 0; i < this->n; ++i) {
    eta[i] = (servers[0].getS() - servers[i].getS()) / servers[i].getW();
  }

  // 第一趟各项参数
  vector<double> delta(this->n, 0);
  vector<double> phi(this->n, 0);
  vector<double> epsilon(this->n, 0);

  // cal Delta
  for (int i = 1; i < this->n; ++i) {
    delta[i] =
        (servers[i - 1].getW() + servers[i - 1].getG() * (this->theta - 1.0)) /
        servers[i].getW();
  }
  for (int i = 1; i < this->n; ++i) {
    Delta[i] = 1.0;
    for (int j = 1; j <= i; j++) {
      Delta[i] *= delta[j];
    }
  }

  // cal Phi
  for (int i = 1; i < this->n; ++i) {
    // phi[i] = servers[i - 1].getG() / servers[i].getW() * this->V / this->Va;
    phi[i] = servers[i - 1].getG() / servers[i].getW();
  }
  for (int i = 1; i < this->n; ++i) {
    Phi[i] = 0;
    for (int j = 1; j <= i; j++) {
      double temp = 1.0;
      for (int k = j + 1; k <= i; k++) {
        temp *= delta[k];
      }
      Phi[i] += (phi[j] * mu[j - 1] * temp);
    }
  }

  // cal Epsilon
  for (int i = 1; i < this->n; ++i) {
    epsilon[i] =
        (servers[i - 1].getO() + servers[i - 1].getS() - servers[i].getS()) /
        servers[i].getW();
  }
  for (int i = 1; i < this->n; ++i) {
    Epsilon[i] = 0.0;
    for (int j = 1; j <= i; j++) {
      double temp = 1.0;
      for (int k = j + 1; k <= i; k++) {
        temp *= delta[k];
      }
      Epsilon[i] += ((phi[j] * eta[j - 1] + epsilon[j]) * temp);
    }
    // cout << "Epsilon[i]: " << Epsilon[i] << endl;
  }

  // 最后一趟各项参数
  vector<double> lambda(this->n, 0);
  vector<double> psi(this->n, 0);
  vector<double> rho(this->n, 0);

  // cal Lambda
  for (int i = 1; i < this->n; ++i) {
    lambda[i] = (servers[i - 1].getW() - servers[i - 1].getG() +
                 servers[i - 1].getG() * this->theta) /
                servers[i].getW();
  }
  for (int i = 1; i < this->n; ++i) {
    Lambda[i] = 1.0;
    for (int j = 1; j <= i; j++) {
      Lambda[i] *= lambda[j];
    }
  }

  // cal Psi
  for (int i = 1; i < this->n; ++i) {
    psi[i] = (servers[i - 1].getG() * this->theta) / servers[i].getW();
    // cout << "psi[i] = " << psi[i] << endl;
  }
  // cout << "this->V / this->Vb = " << this->V / this->Vb << endl;
  for (int i = 1; i < this->n; ++i) {
    Psi[i] = 0;
    for (int j = 1; j <= i; j++) {
      double value = 1.0;
      for (int k = j + 1; k <= i; k++) {
        value *= lambda[k];
      }
      Psi[i] += ((psi[j] * mu[j - 1]) * value);
      // cout << "Psi[i] = " << Psi[i] << endl;
    }
  }

  // cal Rho
  for (int i = 1; i < this->n; ++i) {
    rho[i] =
        (servers[i - 1].getS() - servers[i - 1].getO() - servers[i].getS()) /
        servers[i].getW();
  }
  for (int i = 1; i < this->n; ++i) {
    P[i] = 0;
    for (int j = 1; j <= i; j++) {
      double value = 1.0;
      for (int k = j + 1; k <= i; k++) {
        value *= lambda[k];
      }
      P[i] += ((rho[j] - psi[j] * eta[j - 1]) * value);
    }
    // cout << "P[i]: " << P[i] << endl;
  }

  // lambda <= 1时，修改最后一趟方程
  if (this->l_2 < 1) {
    // cal Lambda *
    for (int i = 1; i < this->n; ++i) {
      lambda[i] =
          (servers[i - 1].getW() + servers[i - 1].getG() * this->theta) /
          servers[i].getW();
    }
    for (int i = 1; i < this->n; ++i) {
      Lambda[i] = 1.0;
      for (int j = 1; j <= i; j++) {
        Lambda[i] *= lambda[j];
      }
    }

    // cal Psi *
    for (int i = 1; i < this->n; ++i) {
      // psi[i] = (servers[i - 1].getG() * (this->theta + 1)) /
      // servers[i].getW() * (this->V / this->Vb);
      psi[i] = (servers[i - 1].getG() * (this->theta + 1)) / servers[i].getW();
      // cout << "psi[i] = " << psi[i] << endl;
    }
    // cout << "this->V / this->Vb = " << this->V / this->Vb << endl;
    for (int i = 1; i < this->n; ++i) {
      Psi[i] = 0;
      for (int j = 1; j <= i; j++) {
        double value = 1.0;
        for (int k = j + 1; k <= i; k++) {
          value *= lambda[k];
        }
        Psi[i] += ((psi[j] * mu[j - 1]) * value);
        // cout << "Psi[i] = " << Psi[i] << endl;
      }
    }

    // cal Rho *
    for (int i = 1; i < this->n; ++i) {
      rho[i] =
          (servers[i - 1].getS() - servers[i - 1].getO() - servers[i].getS()) /
          servers[i].getW();
    }
    for (int i = 1; i < this->n; ++i) {
      P[i] = 0;
      for (int j = 1; j <= i; j++) {
        double value = 1.0;
        for (int k = j + 1; k <= i; k++) {
          value *= lambda[k];
        }
        P[i] += ((rho[j] - psi[j] * eta[j - 1]) * value);
      }
    }
  }

  // cal optimal installment size m
  if (m == 0) {
    calOptimalM();
  }

  // cal beta & alpha & gamma
  calBeta();
  calAlpha();
  calGamma();

  // ensure each load fraction > 0 && all load fraction sum up to 1
  double count_alpha = 0, count_beta = 0, count_gamma = 0;
  for (int i = 0; i < this->n; ++i) {
    count_alpha += alpha[i];
    count_beta += beta[i];
    count_gamma += gamma[i];
    // 此处改为分量 * 总任务量 > 0
    if (alpha[i] * this->Va < 0 || beta[i] * this->V < 0 ||
        gamma[i] * this->Vb < 0) {
      isSchedulable = 0;
      cout << this->W << " error "
           << "alpha * Va - " << alpha[i] * this->Va << " beta * V - "
           << beta[i] * this->V << " gamma * Vb - " << gamma[i] * this->Vb
           << endl;
    }
  }

  // cout << "count_alpha: " << count_alpha << endl;
  // cout << "count_beta: " << count_beta << endl;
  // cout << "count_gamma: " << count_gamma << endl;

  if ((int)((count_alpha + 0.000005) * 100000) != 100000) {
    isSchedulable = 0;
    cout << this->W << " count: alpha - " << count_alpha << endl;
  }

  if ((int)((count_beta + 0.000005) * 100000) != 100000) {
    isSchedulable = 0;
    cout << this->W << " count: beta - " << count_beta << endl;
  }

  if ((int)((count_gamma + 0.000005) * 100000) != 100000) {
    isSchedulable = 0;
    cout << this->W << " count: gamma - " << count_gamma << endl;
  }

  //   if (beta[0] * this->V < gamma[0] * this->Vb) {
  //     isSchedulable = 0;
  //   }
}

/**
 * @brief Cal load fraction beta for every internal installment.
 *
 */
void MISRRLL::calBeta() {
  // cal beta
  double sum_2_n_eta = 0, sum_2_n_mu = 0;
  for (int i = 1; i < this->n; ++i) {
    sum_2_n_mu += mu[i];
    sum_2_n_eta += eta[i];
  }
  double sum = 0;
  for (int i = 0; i < this->n; ++i) {
    if (i == 0) {
      beta[i] = (1.0 - sum_2_n_eta / this->V) / (1.0 + sum_2_n_mu);
      sum += beta[i];
      continue;
    }
    beta[i] = mu[i] * beta[0] + (eta[i] / this->V);
    // 计算每趟内部调度所用时间
    // cout << "i=" << i << " beta[i]: " << beta[i] << endl;
    // cout << "\tserver[i].getS(): " << servers[i].getS() << "\ttime: "
    //      << beta[i] * servers[i].getW() * this->V + servers[i].getS() <<
    //      endl;
    // cout << "this->V=" << this->V << endl;
    sum += beta[i];
  }
  // beta[i]之和为1
  // cout << "this->V: " << this->V << endl;
  // cout << "sum: " << sum << endl;
  // cout << "sum * this->V: " << sum * this->V << endl;
}

/**
 * @brief Cal load fraction alpha for the first installment.
 *
 */
void MISRRLL::calAlpha() {
  // cal alpha
  double sum_2_n_phi = 0, sum_2_n_epsilon = 0, sum_2_n_delta = 0;
  for (int i = 1; i < this->n; ++i) {
    sum_2_n_phi += Phi[i];
    sum_2_n_epsilon += Epsilon[i];
    sum_2_n_delta += Delta[i];
  }
  for (int i = 0; i < this->n; ++i) {
    if (i == 0) {
      // alpha[i] = (1.0 - beta[0] * sum_2_n_phi - (sum_2_n_epsilon / this->Va))
      // / (1.0 + sum_2_n_delta);
      alpha[i] = (1.0 - beta[0] * sum_2_n_phi * this->V / this->Va -
                  (sum_2_n_epsilon / this->Va)) /
                 (1.0 + sum_2_n_delta);
      continue;
    }
    // alpha[i] = Delta[i] * alpha[0] + Phi[i] * beta[0] + (Epsilon[i] /
    // this->Va);
    alpha[i] = Delta[i] * alpha[0] + Phi[i] * beta[0] * this->V / this->Va +
               (Epsilon[i] / this->Va);
    // cout << "i=" << i << " "
    //      << "alpha[i] : " << alpha[i] << endl;
  }
}

/**
 * @brief Cal load fraction gama for the last installment.
 *
 */
void MISRRLL::calGamma() {
  // cal gamma
  double sum_2_n_psi = 0, sum_2_n_p = 0, sum_2_n_lambda = 0;
  for (int i = 1; i < this->n; ++i) {
    sum_2_n_psi += Psi[i];
    sum_2_n_p += P[i];
    sum_2_n_lambda += Lambda[i];
  }
  for (int i = 0; i < this->n; ++i) {
    if (i == 0) {
      // gamma[i] = (1.0 + beta[0] * sum_2_n_psi - (sum_2_n_p / this->Vb)) /
      // (1.0 + sum_2_n_lambda);
      gamma[i] = (1.0 + beta[0] * sum_2_n_psi * (this->V / this->Vb) -
                  (sum_2_n_p / this->Vb)) /
                 (1.0 + sum_2_n_lambda);
      continue;
    }
    gamma[i] = Lambda[i] * gamma[0] - Psi[i] * beta[0] * (this->V / this->Vb) +
               P[i] / this->Vb;
    // cout << "gamma[i]: " << gamma[i] << endl;
  }
}

/**
 * @brief Cal optimal installment size m using formula (4.35).
 *
 */
void MISRRLL::calOptimalM() {
  double A = 0, B = 0, C = 0, D = 0;
  double A_prime = 0, C_prime = 0;
  double sum_2_n_mu = 0, sum_2_n_eta = 0;
  double sum_2_n_delta = 0, sum_2_n_phi = 0, sum_2_n_epsilon = 0;
  double sum_2_n_lambda = 0, sum_2_n_psi = 0, sum_2_n_p = 0;
  for (int i = 1; i < this->n; i++) {
    sum_2_n_mu += mu[i];
    sum_2_n_eta += eta[i];

    sum_2_n_delta += Delta[i];
    sum_2_n_phi += Phi[i];
    sum_2_n_epsilon += Epsilon[i];

    sum_2_n_lambda += Lambda[i];
    sum_2_n_psi += Psi[i];
    sum_2_n_p += P[i];
  }

  double GAMMA = 1 / (1.0 + sum_2_n_delta);
  double ZETA =
      (sum_2_n_epsilon * (1 + sum_2_n_mu) - sum_2_n_phi * sum_2_n_eta) /
      ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_delta));
  double UPSILON = -sum_2_n_phi / ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_delta));

  double IOT = 1 / (1.0 + sum_2_n_mu);
  double KAPPA = sum_2_n_eta / (1.0 + sum_2_n_mu);

  double CHI = 1 / (1.0 + sum_2_n_lambda);
  double OMEGA = (sum_2_n_p * (1 + sum_2_n_mu) + sum_2_n_psi * sum_2_n_eta) /
                 ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_lambda));
  double TAU = sum_2_n_psi / ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_lambda));

  // cout << "TAU: " << TAU << endl;
  // cout << "UPSILON: " << UPSILON << endl;

  // cal A
  A += (servers[0].getG() * theta + servers[0].getW()) * CHI;
  for (int i = 1; i < this->n; ++i) {
    A += (servers[i].getG() * theta * Lambda[i] * CHI);
  }
  A *= l_1;
  A += (servers[0].getW() * GAMMA * l_2);
  A -= (servers[0].getW() * IOT * (l_1 + l_2));
  A *= this->W;

  // cal B
  B = servers[0].getS() - servers[0].getW() * KAPPA;

  // cal C
  C += (servers[0].getG() * theta * TAU);
  // cout << C << endl;
  for (int i = 1; i < this->n; ++i) {
    C += (servers[i].getG() * theta * (Lambda[i] * TAU - Psi[i] * IOT));
    // cout << "Lambda[i] * TAU - Psi[i] * IOT" << Lambda[i] * TAU - Psi[i] *
    // IOT << endl;
  }
  // cout << C << endl;
  C += servers[0].getW() * (UPSILON + TAU);
  // cout << C << endl;
  C *= this->W;

  // cal D
  D += servers[0].getW() * (-ZETA - OMEGA + 2 * KAPPA + this->W * IOT);
  D -= servers[0].getG() * OMEGA * theta;
  for (int i = 0; i < this->n; ++i) {
    D += servers[i].getO();
  }
  for (int i = 1; i < this->n; ++i) {
    D +=
        servers[i].getG() * theta * (Psi[i] * KAPPA - Lambda[i] * OMEGA + P[i]);
  }

  // cal A_prime and C_prime
  A_prime = A + (l_1 + l_2) / 2 * C;
  C_prime = (2 - l_1 - l_2) / 2 * C;

  //// cal T1
  // double T1 = A_prime / m + C_prime / (m - 2) + B * m;
  // cout << "```````````` cal optimal M ````````````" << endl;
  // cout << "A_prime: " << A_prime << ", C_prime: " << C_prime << ", B: " << B
  //      << ", D: " << D << endl;
  // cout << "T1: " << T1 << endl;
  // cout << "T(m): " << T1 + D << endl;

  this->m = findFirstPositiveRealRoot(B, 2 * B, -A - B - C, 3 * A + 2 * B - C,
                                      -2 * A);
  if (this->m < 3)
    this->m = 3;
  cout << "cal optimal m: " << m << endl;
}

/**
 * @brief Cal optimal makespan & every server's using time.
 *
 */
void MISRRLL::getOptimalModel() {
  // cal optimal time(4.15)
  this->optimalTime = 0;
  this->optimalTime += this->startTime;

  this->optimalTime += servers[0].getO();
  this->optimalTime +=
      (servers[0].getS() + this->Va * alpha[0] * servers[0].getW());
  this->optimalTime += (this->m - 2) * (servers[0].getS() +
                                        this->V * beta[0] * servers[0].getW());
  this->optimalTime +=
      (servers[0].getS() + this->Vb * gamma[0] * servers[0].getW());

  this->optimalTime += servers[0].getG() * gamma[0] * this->Vb * this->theta;
  for (int i = 1; i < this->n; i++) {
    this->optimalTime += servers[i].getO();
    this->optimalTime += gamma[i] * this->Vb * servers[i].getG() * this->theta;
  }

  // update usingTime
  for (int i = 0; i < this->n; i++) {
    int id = servers[i].getId();
    usingTime[id] +=
        (alpha[i] * Va * servers[i].getW() + servers[i].getS() +
         (m - 2) * (beta[i] * this->V * servers[i].getW() + servers[i].getS()) +
         gamma[i] * Vb * servers[i].getW() + servers[i].getS() +
         gamma[i] * Vb * servers[i].getG() * theta);
  }

  // cal using rate
  usingRate = 0.0;
  for (int i = 0; i < this->serversNumberWithoutError; ++i) {
    // cout << "usingTime[" << i << "] = " << usingTime[i] << endl;
    usingRate += ((double)usingTime[i] / ((double)(this->optimalTime) *
                                          this->serversNumberWithoutError));
  }

  // for (int i = 0; i < this->n; i++) {
  //     cout << "alpha[" << i << "] = " << alpha[i];
  //     cout << ", beta[" << i << "] = " << beta[i];
  //     cout << ", gamma[" << i << "] = " << gamma[i] << endl;
  // }
}

/**
 * @brief Read o & s & g & w & WTotal from "/data" and set every server.
 *
 */
void MISRRLL::getDataFromFile(string data_path) {
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

/**
 * @brief Output scheduling's using time & using rate to /output.
 *
 */
void MISRRLL::printResult() {
  FILE *fpResult;
  fpResult = fopen(("../output/" + outputName + ".txt").c_str(), "a+");

  if (fpResult == nullptr) {
    printf("The file results.txt can not be opened:\n");
    exit(-1);
  }

  fprintf(fpResult, "time: %lf : %lf + %lf + %lf\n", this->optimalTime,
          this->startTime, this->timeGap,
          this->optimalTime - this->startTime - this->timeGap);
  fprintf(fpResult, "installment: %d : %d + %d\n", this->m,
          this->beforeInstallment, this->m - this->beforeInstallment);
  fprintf(fpResult, "workload: %lf : %lf + %lf\n", this->WWithoutError,
          this->WWithoutError - this->W, this->W);

  // print usingTime
  fprintf(fpResult, "using time: (server id : using time)\n");
  usingRate = 0.0;
  for (int i = 0; i < this->serversNumberWithoutError; ++i) {
    fprintf(fpResult, "%d : %lf\n", i, usingTime[i]);
    usingRate += ((double)usingTime[i] / ((double)(this->optimalTime) *
                                          this->serversNumberWithoutError));
  }

  // print usingRate
  fprintf(fpResult, "using rate: %lf\n", usingRate);
  fprintf(fpResult, "\n\n");

  fclose(fpResult);
}

void MISRRLL::setW(double value) {
  this->W = value;
  this->WWithoutError = value;
}

/**
 * @brief: Set 2 lambdas.
 * @param {double} value1 first installment lambda
 * @param {double} value2 last installment lambda
 * @return {*}
 */
void MISRRLL::setLambda(double value1, double value2) {
  this->l_1 = value1;
  this->l_2 = value2;
}

double MISRRLL::getOptimalTime() const { return this->optimalTime; }

double MISRRLL::getUsingRate() const { return this->usingRate; }

int MISRRLL::getOptimalM() { return this->m; }

/**
 * @brief Calculate the idle time between servers in last installment.
 * Output result to "/output/last_installment_gap.txt".
 *
 */
void MISRRLL::theLastInstallmentGap(string title) {
  vector<double> freeTime(this->n, 0);  // 倒数第二趟计算结束时间
  vector<double> startTime(this->n, 0); // 最后一趟计算开始时间
  vector<double> timeGap(this->n, 0);

  double preTime = 0;
  for (int i = 0; i < this->n; ++i) {
    preTime += servers[i].getO();
    freeTime[i] = preTime;
    preTime += (servers[i].getO() +
                (1 + theta) * servers[i].getG() * this->V * beta[i]);
  }

  preTime = 0;
  for (int i = 0; i < this->n; ++i) {
    preTime += servers[i].getO();
    startTime[i] = preTime;
    preTime +=
        (servers[i].getO() + theta * servers[i].getG() * this->V * beta[i] +
         servers[i].getG() * Vb * gamma[i]);
  }

  if (this->l_2 < 1) {
    double preTime = 0;
    for (int i = 0; i < this->n; ++i) {
      preTime += servers[i].getO();
      freeTime[i] = preTime;
      preTime += (servers[i].getO() + servers[i].getG() * Vb * gamma[i] +
                  theta * servers[i].getG() * this->V * beta[i]);
    }

    preTime = 0;
    for (int i = 0; i < this->n; ++i) {
      preTime += servers[i].getO();
      startTime[i] = preTime;
      preTime += (servers[i].getO() +
                  (1 + theta) * servers[i].getG() * this->V * beta[i]);
    }
  }

  FILE *fresult;
  string fname =
      "../output/MISRRLL/last_installment_gap_lambda" + title + ".txt";
  fresult = fopen(fname.c_str(), "w");

  for (int i = 0; i < this->n; ++i) {
    timeGap[i] = startTime[i] - freeTime[i];
    if (timeGap[i] < 0) {
      isSchedulable = 0;
    }
    fprintf(fresult, "gap[%d]: \t%.2f\n", i, timeGap[i]);
  }
}

void MISRRLL::addServer(int installment, vector<Server> &new_servers) {
  // 重调度开始时间
  this->startTime = servers[0].getO() + servers[0].getS() * (installment + 1) +
                    servers[0].getW() * alpha[0] * this->Va +
                    servers[0].getW() * beta[0] * installment * this->V +
                    servers[0].getG() * beta[0] * this->V * theta;

  // 第x+1趟各处理机释放时间(相对于startTime)
  vector<double> release_time(this->n + 1, 0);
  release_time[1] = -servers[1].getG() * beta[1] * this->V * theta;
  for (int i = 2; i < this->n + 1; i++) {
    release_time[i] =
        release_time[i - 1] + servers[i - 1].getO() +
        servers[i - 1].getG() * alpha[i - 1] * this->Va * (1 + theta) +
        servers[i - 1].getG() * beta[i - 1] * this->V * (1 + theta) +
        servers[i].getO();
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
  this->leftW += (this->Vb + (this->m - installment - 2) * this->V);
  this->W = this->leftW;
  this->m = this->m - installment - 1;
  // this->m = this->m - installment - 17;
  //   this->m = 12;
  this->n += num_new_servers;

  initValue();
  getOptimalModel();
  this->optimalTime += this->startTime;
  cout << "new optimal time: " << this->optimalTime << endl;

  // 重调度第一趟各处理机开始计算的时间
  vector<double> restart_time(this->n, 0);
  restart_time[0] = servers[0].getO();
  for (int i = 1; i < this->n; i++) {
    restart_time[i] = restart_time[i - 1] +
                      servers[i - 1].getG() * alpha[i - 1] * this->Va +
                      servers[i].getO();
  }

  vector<double> gap(this->n, 0);
  for (int i = 0; i < this->n; i++) {
    gap[i] = restart_time[i] - release_time[i];
    std::cout << release_time[i] << ", " << restart_time[i] << std::endl;
  }

  double waiting_time = *(std::min_element(gap.begin(), gap.end()));
  std::cout << "wating_time: " << waiting_time << std::endl;
  waiting_time = waiting_time < 0 ? -waiting_time : 0;
  std::cout << "wating_time: " << waiting_time << std::endl;
  std::cout << "final optimal time: " << this->optimalTime + waiting_time
            << std::endl;

  for (int i = 0; i < this->n; i++) {
    std::cout << alpha[i] << ", ";
    std::cout << beta[i] << ", ";
    std::cout << gamma[i] << std::endl;
  }
}

void MISRRLL::error(vector<int> &errorPlace, int errorInstallment) {
  // 最优参数
  this->l_1 = 1;
  this->l_2 = 1;
  this->errorPlace = errorPlace;
  this->errorInstallment = errorInstallment;

  outputName = outputName + '_' + to_string((int)errorPlace.size()) +
               "_installment_" + to_string(errorInstallment) + "_ server";
  for (auto i : errorPlace) {
    outputName = outputName + '_' + to_string(i);
  }

  // cal re-schedule start time
  if (errorInstallment == 1) {
    this->startTime = servers[0].getS() * 2 +
                      servers[0].getW() * (alpha[0] + beta[0]) * this->V;
  } else if (errorInstallment == this->m || errorInstallment == this->m - 1) {
    this->startTime = optimalTime;
  } else {
    // 。。。很难搞
    this->startTime =
        servers[0].getS() * (errorInstallment + 1) +
        servers[0].getW() *
            (alpha[0] * this->Va + beta[0] * errorInstallment * this->V);
  }

  // cal each server's release time
  vector<double> busyTime(this->n, 0);
  // double preTime = servers[0].getO();
  double preTime = 0;
  for (int i = 0; i < this->n; ++i) {
    if (i != 0) {
      preTime += (servers[i - 1].getO() +
                  servers[i - 1].getG() * beta[i - 1] * this->V * (1 + theta));
      // servers[i].getO());
    }
    busyTime[i] = preTime;
  }

  // cal each server's using time
  for (int i = 0; i < n; ++i) {
    auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
    if (iter == errorPlace.end()) {
      // 正常处理机
      usingTime[i] =
          (servers[i].getS() + alpha[i] * this->Va * servers[i].getW() +
           errorInstallment *
               (servers[i].getS() + beta[i] * this->V * servers[i].getW()));
    } else {
      // 出错的机器
      usingTime[i] =
          (servers[i].getS() + alpha[i] * this->Va * servers[i].getW() +
           (errorInstallment - 2) *
               (servers[i].getS() + beta[i] * this->V * servers[i].getW()));
    }
  }

  // cal left operating server
  int position = 0;
  vector<Server> oldServers = servers;
  servers.clear();
  for (int i = 0; i < oldServers.size(); ++i) {
    auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
    if (iter == errorPlace.end()) {
      // 这里储存了处理机的原编号
      servers[position] = oldServers[i];
      position++;
    }
  }
  this->n -= (int)errorPlace.size();

  // cal left workload
  this->leftW = 0.0;
  if (errorInstallment == 1) {
    this->leftW += ((this->m - errorInstallment - 1) * this->V);
    for (auto i : errorPlace) {
      this->leftW += (alpha[i - 1] * this->Va + beta[i - 1] * this->V);
    }
  } else if (errorInstallment == this->m - 1) {
    for (auto i : errorPlace) {
      this->leftW += (beta[i - 1] * this->V + gamma[i - 1] * this->Vb);
    }
  } else if (errorInstallment == this->m) {
    for (auto i : errorPlace) {
      this->leftW += gamma[i - 1] * this->Vb;
    }
  } else {
    // this->leftW += ((this->m - errorInstallment - 1) * this->V);
    // 适配新模型
    this->leftW += ((this->m - errorInstallment - 2) * this->V + this->Vb);
    for (auto i : errorPlace) {
      this->leftW += (2 * beta[i - 1] * this->V);
    }
  }

  cout << "left workload: " << this->leftW << endl;

  // re-schedule
  this->W = this->leftW;
  beforeInstallment = errorInstallment + 1;

  old_V = V;
  old_beta = beta;

  m = 0;
  initValue();

  cout << "error optimal m=" << m << endl;
  cout << "error this->optimalTime=" << this->optimalTime << endl;
  cout << "this->startTime=" << this->startTime << endl;

  getOptimalModel();

  // this->optimalTime += this->startTime;
  this->m += beforeInstallment;

  // cal each server's restart time
  position = 0;
  // preTime = servers[0].getO();
  preTime = 0;
  vector<double> reTime(serversNumberWithoutError, 0);
  for (int i = 0; i < this->n; ++i) {
    if (i != 0) {
      preTime += (servers[i - 1].getO() +
                  servers[i - 1].getG() * alpha[i - 1] * this->Va);
    }
    // cout << "i=" << i << " servers[i].getId()=" << servers[i].getId() <<
    // endl;
    reTime[servers[i].getId()] = preTime;
  }

  // cal each server's conflict time
  vector<double> codeTimeGap(serversNumberWithoutError, 0);
  for (int i = 0; i < serversNumberWithoutError; ++i) {
    if (reTime[i] != 0)
      codeTimeGap[i] = busyTime[i] - reTime[i];
    // 验证没什么问题
    // cout << "Ts time[" << i << "]: " << reTime[i] << endl;
    // cout << "Tf time[" << i << "]: " << busyTime[i] << endl;
    // cout << "conflict time[" << i << "]: " << codeTimeGap[i] << endl;
  }

  // // mathTimeGap
  // vector<double> mathTimeGap(serversNumberWithoutError, 0);
  // preTime = 0, position = 1;
  // for (int i = 0; i < serversNumberWithoutError; ++i) {
  //   if (i == 0)
  //     continue;
  //   if (find(errorPlace.begin(), errorPlace.end(), i + 1) !=
  //   errorPlace.end()) {
  //     preTime +=
  //         (2 * oldServers[i - 1].getO() +
  //          (theta + 1.0) * oldServers[i - 1].getG() * old_beta[i - 1] *
  //          old_V);
  //     continue;
  //   }
  //   int demo = servers[position].getId();
  //   preTime +=
  //       (2 * oldServers[i - 1].getO() +
  //        (theta + 1.0) * oldServers[i - 1].getG() * old_beta[i - 1] * old_V -
  //        servers[position - 1].getO() -
  //        servers[position - 1].getG() * alpha[position - 1] * this->Va);

  //   mathTimeGap[demo] = preTime;
  //   ++position;
  // }

  // compare codeTimeGap & mathTimeGap
  // cout << "code: " << codeTimeGap[serversNumberWithoutError - 1] <<
  //     "math: " << mathTimeGap[serversNumberWithoutError - 1] << endl;

  // cal server[0]'s optimal waiting time
  timeGap = *max_element(codeTimeGap.begin(), codeTimeGap.end());

  // get whole schedule's optimal makespan
  this->optimalTime += timeGap;
}

int MISRRLL::findFirstPositiveRealRoot(double a, double b, double c, double d,
                                       double e) {
  // 定义四次方程的系数
  Eigen::VectorXd coefficients(5);
  coefficients << e, d, c, b, a;
  // cout << e << " " << d << " " << c << " " << b << " " << a << endl;

  // 使用 Eigen 的多项式求根功能
  Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
  solver.compute(coefficients);

  // 获取方程的根（直接使用 Eigen 的 RootsType 遍历）
  const auto &roots = solver.roots();

  // 遍历根，筛选出第一个大于0的实数根
  for (int i = 0; i < roots.size(); ++i) {
    const auto &root = roots[i];
    // cout << root.real() << " + " << root.imag() << "i" << endl;
    if (std::abs(root.imag()) < 1e-8 &&
        root.real() > 3) { // 判断是否为大于3的实数根
      return std::lround(root.real());
    }
  }

  // 如果没有找到符合条件的根，返回 std::nullopt
  return 0;
}

// /**
//  * @brief Calculate the idle time between servers in first installment.
//  * Output result to "/output/first_installment_gap.txt".
//  *
//  */
// void MISRRLL::theFirstInstallmentGap(string title) {
//     vector<double> freeTime(this->n, 0);    // 第一趟计算结束时间
//     vector<double> startTime(this->n, 0);   // 第二趟计算开始时间
//     vector<double> timeGap(this->n, 0);

//     double preTime = 0;
//     for (int i = 0; i < this->n; ++i) {
//         preTime += servers[i].getO();
//         freeTime[i] = preTime;
//         preTime += (servers[i].getO() + servers[i].getG() * Va * alpha[i]);
//     }

//     preTime = 0;
//     for (int i = 0; i < this->n; ++i) {
//         preTime += servers[i].getO();
//         startTime[i] = preTime;
//         preTime += (servers[i].getO() + theta * servers[i].getG() * V *
//         beta[i] +
//                 servers[i].getG() * Vb * gamma[i]);
//     }

//     if(this->l_2 < 1) {
//         double preTime = 0;
//         for (int i = 0; i < this->n; ++i) {
//             preTime += servers[i].getO();
//             freeTime[i] = preTime;
//             preTime += (servers[i].getO() + servers[i].getG() * Vb * gamma[i]
//             + theta * servers[i].getG() * V * beta[i]);
//         }

//         preTime = 0;
//         for (int i = 0; i < this->n; ++i) {
//             preTime += servers[i].getO();
//             startTime[i] = preTime;
//             preTime += (servers[i].getO() + (1 + theta) * servers[i].getG() *
//             V * beta[i]);
//         }
//     }

//     FILE *fresult;
//     string fname = "../output/MISRRLL/first_installment_gap_lambda" + title +
//     ".txt"; fresult = fopen(fname.c_str(), "w");

//     for (int i = 0; i < this->n; ++i) {
//         timeGap[i] = startTime[i] - freeTime[i];
//         if (timeGap[i] < 0) {
//             isSchedulable = 0;
//         }
//         fprintf(fresult, "gap[%d]: \t%.2f\n", i, timeGap[i]);
//     }
// }