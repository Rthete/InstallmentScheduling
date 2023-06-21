#include "MISRRL.h"

MISRRL::MISRRL(int valueN, double valueTheta, int installment) :
        n(valueN),
        m(installment),
        theta(valueTheta),
        usingRate(0),
        serversNumberWithoutError(valueN),
        servers(vector<Server>(MAX_N)),

        //beta
        beta(vector<double>(MAX_N, 0)),
        mu(vector<double>(MAX_N, 0)),
        eta(vector<double>(MAX_N, 0)),

        // alpha
        alpha(vector<double>(MAX_N, 0)),
        Delta(vector<double>(MAX_N, 0)),
        Phi(vector<double>(MAX_N, 0)),
        Epsilon(vector<double>(MAX_N, 0)),

        // gamma
        gamma(vector<double>(MAX_N, 0)),
        Lambda(vector<double>(MAX_N, 0)),
        Psi(vector<double>(MAX_N, 0)),
        P(vector<double>(MAX_N, 0)),

        // time gap
        old_beta(vector<double>(MAX_N)),

        // time
        usingTime(vector<double>(MAX_N, 0)),
        outputName("MISRRL_print_result"){}

MISRRL::MISRRL(int valueN, double valueTheta) :
        n(valueN),
        theta(valueTheta),
        usingRate(0),
        serversNumberWithoutError(valueN),
        servers(vector<Server>(MAX_N)),

        //beta
        beta(vector<double>(MAX_N, 0)),
        mu(vector<double>(MAX_N, 0)),
        eta(vector<double>(MAX_N, 0)),

        // alpha
        alpha(vector<double>(MAX_N, 0)),
        Delta(vector<double>(MAX_N, 0)),
        Phi(vector<double>(MAX_N, 0)),
        Epsilon(vector<double>(MAX_N, 0)),

        // gamma
        gamma(vector<double>(MAX_N, 0)),
        Lambda(vector<double>(MAX_N, 0)),
        Psi(vector<double>(MAX_N, 0)),
        P(vector<double>(MAX_N, 0)),

        // time gap
        old_beta(vector<double>(MAX_N)),

        // time
        usingTime(vector<double>(MAX_N, 0)),
        outputName("MIS-CRR") {}

/**
 * @brief Assign model.
 * 
 */
void MISRRL::initValue() {

    // cal load for each installment
    // V: 第一趟及内部调度的每趟任务量
    this->V = (this->m - this->l) * this->W / (this->m * (this->m - 1));
    // Vb: 最后一趟的任务量
    this->Vb = this->l / this->m * this->W;

    // cal mu (beta)
    for (int i = 0; i < this->n; ++i) {
        mu[i] = servers[0].getW() / servers[i].getW();
    }

    // cal eta (beta)
    for (int i = 0; i < this->n; ++i) {
        eta[i] = (servers[0].getS() - servers[i].getS()) / servers[i].getW();
    }

    vector<double> delta(this->n, 0);
    vector<double> phi(this->n, 0);
    vector<double> epsilon(this->n, 0);

    // cal Delta
    for (int i = 1; i < this->n; ++i) {
        delta[i] = (servers[i - 1].getW() + servers[i - 1].getG() * (this->theta - 1.0)) / servers[i].getW();
    }
    for (int i = 1; i < this->n; ++i) {
        Delta[i] = 1.0;
        for (int j = 1; j <= i; j++) {
            Delta[i] *= delta[j];
        }
    }

    // cal Phi
    for (int i = 1; i < this->n; ++i) {
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
        epsilon[i] = (servers[i - 1].getO() + servers[i - 1].getS() - servers[i].getS()) / servers[i].getW();
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
    }

    // 最后一趟各项参数
    vector<double> lambda(this->n, 0);
    vector<double> psi(this->n, 0);
    vector<double> rho(this->n, 0);

    // cal Lambda
    for (int i = 1; i < this->n; ++i) {
        lambda[i] = (servers[i - 1].getW() - servers[i - 1].getG() + servers[i - 1].getG() * this->theta) / servers[i].getW();
    }
    for (int i = 1; i < this->n; ++i) {
        Lambda[i] = 1.0;
        for (int j = 1; j <= i; j++) {
            Lambda[i] *= lambda[j];
        }
    }

    // cal Psi
    for (int i = 1; i < this->n; ++i) {
        psi[i] = (servers[i - 1].getG() * this->theta) / servers[i].getW() * (this->V / this->Vb);
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
        rho[i] = (servers[i - 1].getS() - servers[i - 1].getO() - servers[i].getS()) / servers[i].getW();
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

    // lambda < 1时，修改最后一趟方程
    if(this->l < 1) {
        // cal Lambda *
        for (int i = 1; i < this->n; ++i) {
            lambda[i] = (servers[i - 1].getW() + servers[i - 1].getG() * this->theta) / servers[i].getW();
        }
        for (int i = 1; i < this->n; ++i) {
            Lambda[i] = 1.0;
            for (int j = 1; j <= i; j++) {
                Lambda[i] *= lambda[j];
            }
        }

        // cal Psi *
        for (int i = 1; i < this->n; ++i) {
            psi[i] = (servers[i - 1].getG() * (this->theta + 1)) / servers[i].getW() * (this->V / this->Vb);
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
            rho[i] = (servers[i - 1].getS() - servers[i - 1].getO() - servers[i].getS()) / servers[i].getW();
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
        if (alpha[i] * this->V < 0 || beta[i] * this->V < 0 || gamma[i] * this->Vb < 0) {
            isSchedulable = 0;
            cout << this->W << " error " << "alpha * V - " << alpha[i] * this->V << " beta * V - " << beta[i] * this->V
                << " gamma * Vb - " << gamma[i] * this->Vb << endl;
        }
    }
    if ((int)((count_alpha + 0.000005) * 100000) != 100000) {
        isSchedulable = 0;
        cout << this->W << " count: alpha - " << count_alpha << endl;
    }
        
    if ((int)((count_beta + 0.000005) * 100000)  != 100000) {
        isSchedulable = 0;
        cout << this->W << " count: beta - " << count_beta << endl;
    }

    if ((int)((count_gamma + 0.000005) * 100000)  != 100000) {
        isSchedulable = 0;
        cout << this->W << " count: gamma - " << count_gamma << endl;
    }
}

/**
 * @brief Cal load fraction beta for every internal installment.
 * 
 */
void MISRRL::calBeta() {

    //cal beta
    double sum_2_n_eta = 0, sum_2_n_mu = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_mu += mu[i];
        sum_2_n_eta += eta[i];
    }
    double sum = 0;
    for (int i = 0 ; i < this->n; ++i) {
        if (i == 0) {
            beta[i] = (1.0 - sum_2_n_eta / this->V) / (1.0 + sum_2_n_mu);
            sum += beta[i];
            continue;
        }
        beta[i] = mu[i] * beta[0] + (eta[i] / this->V);
        // 计算每趟内部调度所用时间
        // cout << "beta[i]: " << beta[i] << endl;
        // cout << "\tserver[i].getS(): " << servers[i].getS() << "\ttime: " 
        //     << beta[i] * servers[i].getW() * this->V + servers[i].getS() << endl;

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
void MISRRL::calAlpha() {

    // cal alpha
    double sum_2_n_phi = 0, sum_2_n_epsilon = 0, sum_2_n_delta = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_phi += Phi[i];
        sum_2_n_epsilon += Epsilon[i];
        sum_2_n_delta += Delta[i];
    }
    for (int i = 0; i < this->n; ++i) {
        if (i == 0) {
            alpha[i] = (1.0 - beta[0] * sum_2_n_phi - (sum_2_n_epsilon / this->V)) / (1.0 + sum_2_n_delta);
            continue;
        }
        alpha[i] = Delta[i] * alpha[0] + Phi[i] * beta[0] + (Epsilon[i] / this->V);
        // cout << "alpha[i]: " << alpha[i] << endl;
    }
}

/**
 * @brief Cal load fraction gama for the last installment.
 * 
 */
void MISRRL::calGamma() {
    // cal gamma
    double sum_2_n_psi = 0, sum_2_n_p = 0, sum_2_n_lambda = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_psi += Psi[i];
        sum_2_n_p += P[i];
        sum_2_n_lambda += Lambda[i];
    }
    for (int i = 0; i < this->n; ++i) {
        if (i == 0) {
            gamma[i] = (1.0 + beta[0] * sum_2_n_psi - (sum_2_n_p / this->Vb)) / (1.0 + sum_2_n_lambda);
            continue;
        }
        gamma[i] = Lambda[i] * gamma[0] - Psi[i] * beta[0] + P[i] / this->Vb;
        // cout << "gamma[i]: " << gamma[i] << endl;
    }
}

/**
 * @brief Cal optimal installment size m using formula (4.35).
 * 
 */
void MISRRL::calOptimalM() {
    double A = 0, B = 0;
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

    double GAMMA = ((1.0 + sum_2_n_mu - sum_2_n_phi) * (1.0 + sum_2_n_lambda)) /
            ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_delta) * (1.0 + sum_2_n_lambda));
    double IOT = ((1.0 + sum_2_n_lambda) * (1.0 + sum_2_n_delta)) /
            ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_delta) * (1.0 + sum_2_n_lambda));
    double CHI = ((1.0 + sum_2_n_mu + sum_2_n_psi) * (1.0 + sum_2_n_delta)) /
            ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_delta) * (1.0 + sum_2_n_lambda));
    double KAPPA = (sum_2_n_eta * (1.0 + sum_2_n_lambda) * (1.0 + sum_2_n_delta)) /
            ((1.0 + sum_2_n_mu) * (1.0 + sum_2_n_delta) * (1.0 + sum_2_n_lambda));

    // cal A
    A += servers[0].getG() * theta * CHI;
    for (int i = 1; i < this->n; ++i) {
        A += (servers[i].getG() * theta * (Lambda[i] * CHI - Psi[i] * IOT));
    }
    A += ((GAMMA + CHI - 2 * IOT) * servers[0].getW());
    A *= this->W;

    // cal B
    B = servers[0].getS() - servers[0].getW() * KAPPA;

    // cal m
    this->m = (int)(sqrt(1.0 / 4.0 + A / B) + 1.0 / 2.0);
    // cout << "optimal m = " << this->m << endl;
}

/**
 * @brief Cal optimal makespan & every server's using time.
 * 
 */
void MISRRL::getOptimalModel() {
    
    // cal optimal time(4.15)
    this->optimalTime = 0;
    this->optimalTime += (servers[0].getS() + this->V * alpha[0] * servers[0].getW());
    this->optimalTime += (this->m - 2) * (servers[0].getS() + this->V * beta[0] * servers[0].getW());
    this->optimalTime += (servers[0].getS() + this->Vb * gamma[0] * servers[0].getW());

    this->optimalTime += servers[0].getG() * gamma[0] * this->Vb * this->theta;
    for (int i = 1; i < this->n; i++) {
        this->optimalTime += servers[i].getO();
        this->optimalTime += gamma[i] * this->Vb * servers[i].getG() * this->theta;
    }

    
    // update usingTime
    for (int i = 0; i < this->n; i++) {
        int id = servers[i].getId();
        usingTime[id] += (alpha[i] * V * servers[i].getW() + servers[i].getS() +
                (m - 2) * (beta[i] * V * servers[i].getW() + servers[i].getS()) +
                gamma[i] * Vb * servers[i].getW() + servers[i].getS() + gamma[i] * Vb * servers[i].getG() * theta);
    }

    // cal using rate
    usingRate = 0.0;
    for (int i = 0; i < this->serversNumberWithoutError; ++i) {
        usingRate += ((double)usingTime[i] / ((double)(this->optimalTime) * this->serversNumberWithoutError));
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
void MISRRL::getDataFromFile() {
    FILE *fpo, *fps, *fpg, *fpw, *totalW;
    double valueO[this->n], valueS[this->n], valueG[this->n], valueW[this->n];

    fpo = fopen("../data/15-servers-w-20/o.txt", "r");
    fps = fopen("../data/15-servers-w-20/s.txt", "r");
    fpg = fopen("../data/15-servers-w-20/g.txt", "r");
    fpw = fopen("../data/15-servers-w-20/w.txt", "r");
    totalW  = fopen("../data/15-servers-w-20/WTotal.txt", "r");

    if (fpo == nullptr || fps == nullptr || fpg == nullptr || fpw == nullptr || totalW == nullptr) {
        printf("The file can not be opened:\n");
        exit(-1);
    }

    fscanf(totalW, "%lf", &this->W);
    for (int i = 0; fscanf(fpo, "%lf", &valueO[i]) != EOF; ++i);
    for (int i = 0; fscanf(fps, "%lf", &valueS[i]) != EOF; ++i);
    for (int i = 0; fscanf(fpg, "%lf", &valueG[i]) != EOF; ++i);
    for (int i = 0; fscanf(fpw, "%lf", &valueW[i]) != EOF; ++i);

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
void MISRRL::printResult() {
    FILE * fpResult;
    fpResult = fopen(("../output/" + outputName + ".txt").c_str(), "a+");

    if (fpResult == nullptr) {
        printf("The file results.txt can not be opened:\n");
        exit(-1);
    }

    fprintf(fpResult, "time: %lf : %lf + %lf + %lf\n", this->optimalTime, this->startTime, this->timeGap,
            this->optimalTime - this->startTime - this->timeGap);
    fprintf(fpResult, "installment: %d : %d + %d\n", this->m, this->beforeInstallment, this->m - this->beforeInstallment);
    fprintf(fpResult, "workload: %lf : %lf + %lf\n", this->WWithoutError, this->WWithoutError - this->W, this->W);

    // print usingTime
    fprintf(fpResult, "using time: (server id : using time)\n");
    usingRate = 0.0;
    for (int i = 0; i < this->serversNumberWithoutError; ++i) {
        fprintf(fpResult, "%d : %lf\n", i, usingTime[i]);
        usingRate += ((double)usingTime[i] / ((double)(this->optimalTime) * this->serversNumberWithoutError));
    }

    // print usingRate
    fprintf(fpResult, "using rate: %lf\n", usingRate);
    fprintf(fpResult, "\n\n");

    fclose(fpResult);
}

/**
 * @brief Dynamic fault-tolerant algorithm.
 * Calculate the optimal waiting time and output reschedule result.
 * 
 * @param errorPlace 
 * @param errorInstallment 
 */
void MISRRL::error(vector<int> &errorPlace, int errorInstallment) {
    outputName = outputName + '_' + to_string((int)errorPlace.size()) +
            "_installment_" + to_string(errorInstallment) + "_ server";
    for (auto i : errorPlace) {
        outputName = outputName + '_' + to_string(i);
    }

    // cal re-schedule start time
    if (errorInstallment == 1) {
        this->startTime = servers[0].getS() * 2 + servers[0].getW() * (alpha[0] + beta[0]) * this->V;
    } else if (errorInstallment == this->m || errorInstallment == this->m - 1) {
        this->startTime = optimalTime;
    } else {
        this->startTime = servers[0].getS() * (errorInstallment + 1) +
                servers[0].getW() * (alpha[0] + beta[0] * errorInstallment) * this->V;
    }

    // cal each server's release time
    vector<double> busyTime(this->n, 0);
    double preTime = servers[0].getO();
    for (int i = 0; i < this->n; ++i) {
        if (i != 0) {
            preTime += (servers[i - 1].getO() + servers[i - 1].getG() * beta[i - 1] * this->V * (1 + theta) + servers[i].getO());
        }
        busyTime[i] = preTime;
    }

    // cal each server's using time
    for (int i = 0; i < n; ++i) {
        auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
        if (iter == errorPlace.end()) {
            usingTime[i] = (servers[i].getS() + alpha[i] * V * servers[i].getW() +
                    errorInstallment * (servers[i - 1].getS() + beta[i] * V * servers[i].getW()));
        } else {
            usingTime[i] = (servers[i].getS() + alpha[i] * V * servers[i - 1].getW() +
                    errorInstallment * (servers[i].getS() + beta[i] * V * servers[i].getW()));
        }
    }

    // cal left operating server
    int position = 0;
    vector<Server> oldServers = servers;
    servers.clear();
    for (int i = 0; i < oldServers.size(); ++i) {
        auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
        if (iter == errorPlace.end()) {
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
            this->leftW += (alpha[i - 1] * this->V + beta[i - 1] * this->V);
        }
    } else if (errorInstallment == this->m - 1) {
        for (auto i : errorPlace) {
            this->leftW += (beta[i - 1] * this->V + gamma[i - 1] * this->V);
        }
    } else if (errorInstallment == this->m) {
        for (auto i : errorPlace) {
            this->leftW += gamma[i - 1] * this->V;
        }
    } else {
        this->leftW += ((this->m - errorInstallment - 1) * this->V);
        for (auto i : errorPlace) {
            this->leftW += (2 * beta[i - 1] * this->V);
        }
    }

    // re-schedule
    this->W = this->leftW;
    beforeInstallment = errorInstallment + 1;

    old_V = V;
    old_beta = beta;

    m = 0;
    initValue();
    getOptimalModel();

    this->optimalTime += this->startTime;
    this->m += beforeInstallment;

    // cal each server's restart time
    position = 0;
    preTime = servers[0].getO();
    vector<double> reTime(serversNumberWithoutError, 0);
    for (int i = 0; i < this->n; ++i) {
        if (i != 0) {
            preTime += (servers[i].getO() + servers[i - 1].getG() * alpha[i - 1] * V);
        }
        reTime[servers[i].getId()] = preTime;
    }

    // cal each server's conflict time
    vector<double> codeTimeGap(serversNumberWithoutError, 0);
    for (int i = 0; i < serversNumberWithoutError; ++i) {
        if (reTime[i] != 0)
            codeTimeGap[i] = busyTime[i] - reTime[i];
    }

    // mathTimeGap
    vector<double> mathTimeGap(serversNumberWithoutError, 0);
    preTime = 0, position = 1;
    for (int i = 0; i < serversNumberWithoutError; ++i) {
        if (i == 0) continue;
        if (find(errorPlace.begin(), errorPlace.end(), i + 1) != errorPlace.end()) {
            preTime += (2 * oldServers[i - 1].getO() + (theta + 1.0) * oldServers[i - 1].getG() * old_beta[i - 1] * old_V);
            continue;
        }
        int demo = servers[position].getId();
        preTime += (2 * oldServers[i - 1].getO() + (theta + 1.0) * oldServers[i - 1].getG() * old_beta[i - 1] * old_V
                    - servers[position - 1].getO() - servers[position - 1].getG() * alpha[position - 1] * V);

        mathTimeGap[demo] = preTime;
        ++position;
    }

    // compare codeTimeGap & mathTimeGap
    // cout << "code: " << codeTimeGap[serversNumberWithoutError - 1] <<
    //     "math: " << mathTimeGap[serversNumberWithoutError - 1] << endl;

    // cal server[0]'s optimal waiting time
    timeGap = *max_element(codeTimeGap.begin(), codeTimeGap.end());

    // get whole schedule's optimal makespan
    this->optimalTime += timeGap;
}

void MISRRL::setW(double value) {
    this->W = value;
    this->WWithoutError = value;
}

void MISRRL::setLambda(double value) {
    this->l = value;
}

double MISRRL::getOptimalTime() const {
    return  this->optimalTime;
}

double MISRRL::getUsingRate() const {
    return this->usingRate;
}

int MISRRL::getOptimalM() const {
    return this->m;
}

/**
 * @brief Calculate the idle time between servers in last installment.
 * Output result to "/output/last_installment_gap.txt".
 * 
 */
void MISRRL::theLastInstallmentGap(string title) {
    vector<double> freeTime(this->n, 0);    // 倒数第二趟计算结束时间
    vector<double> startTime(this->n, 0);   // 最后一趟计算开始时间
    vector<double> timeGap(this->n, 0);

    double preTime = 0;
    for (int i = 0; i < this->n; ++i) {
        preTime += servers[i].getO();
        freeTime[i] = preTime;
        preTime += (servers[i].getO() + (1 + theta) * servers[i].getG() * V * beta[i]);
    }

    preTime = 0;
    for (int i = 0; i < this->n; ++i) {
        preTime += servers[i].getO();
        startTime[i] = preTime;
        preTime += (servers[i].getO() + theta * servers[i].getG() * V * beta[i] +
                servers[i].getG() * Vb * gamma[i]);
    }

    if(this->l < 1) {
        double preTime = 0;
        for (int i = 0; i < this->n; ++i) {
            preTime += servers[i].getO();
            freeTime[i] = preTime;
            preTime += (servers[i].getO() + servers[i].getG() * Vb * gamma[i] + theta * servers[i].getG() * V * beta[i]);
        }

        preTime = 0;
        for (int i = 0; i < this->n; ++i) {
            preTime += servers[i].getO();
            startTime[i] = preTime;
            preTime += (servers[i].getO() + (1 + theta) * servers[i].getG() * V * beta[i]);
        }
    }

    FILE *fresult;
    string fname = "../output/MISRRL/last_installment_gap_lambda" + title + ".txt";
    fresult = fopen(fname.c_str(), "w");

    for (int i = 0; i < this->n; ++i) {
        timeGap[i] = startTime[i] - freeTime[i];
        if (timeGap[i] < 0) {
            isSchedulable = 0;
        }
        fprintf(fresult, "gap[%d]: \t%.2f\n", i, timeGap[i]);
    }
}
