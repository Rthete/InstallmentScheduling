//
// Created by xd_wxl on 2023/1/17.
//

//
// Created by kqzhang on 2022/12/8.
//

#include "RRMIS.h"

bool cmp(Server & A, Server & B) {
    return A.getW() > B.getW();
}

RRMIS::RRMIS(int valueN, double valueTheta) :
        n(valueN),
        theta(valueTheta),
        serversNumberWithoutError(valueN),
        servers(vector<Server>(MAX_N)),

        //beta
        beta(vector<double>(MAX_N)),
        Mu(vector<double>(MAX_N)),
        Eta(vector<double>(MAX_N)),

        // alpha
        alpha(vector<double>(MAX_N)),
        Delta(vector<double>(MAX_N)),
        Phi(vector<double>(MAX_N)),
        Epsilon(vector<double>(MAX_N)),

        // gamma
        gamma(vector<double>(MAX_N, 0)),
        Lambda(vector<double>(MAX_N, 0)),
        Psi(vector<double>(MAX_N, 0)),
        P(vector<double>(MAX_N, 0)),

        // time
        usingTime(vector<double>(MAX_N, 0)),

        //
        old_beta(vector<double>(MAX_N)){ }

/**
 * assign model
 */
void RRMIS::initValue() {
    /**
      * cal best M & get mean value
      */
//    calOptimalM();
//    this->m = getBestM();
//    this->m = 15;
    this->V = this->W / this->m;

    /**
      * cal alpha & beta & gamma
      */
    calBeta();
    calAlpha();
    calGamma();

    double count_alpha = 0, count_beta = 0, count_gamma = 0;
    for (int i = 0; i < this->n; ++i) {
        count_alpha += alpha[i];
        count_beta += beta[i];
        count_gamma += gamma[i];
        if (alpha[i] < 0 || beta[i] < 0 || gamma[i] < 0) {
            cout << this->W << " error " << "alpha - " << alpha[i] << " beta - " << beta[i]
                 << " gamma - " << gamma[i] << endl;
        }
    }
    if ((int)((count_alpha + 0.000005) * 100000) != 100000)
        cout << this->W << " count: alpha - " << count_alpha << endl;
    if ((int)((count_beta + 0.000005) * 100000)  != 100000)
        cout << this->W << " count: beta - " << count_beta << endl;
    if ((int)((count_gamma + 0.000005) * 100000)  != 100000)
        cout << this->W << " count: gamma - " << count_gamma << endl;
}

void RRMIS::calAlpha() {
    /**
     * cal delta & phi & epsilon (alpha)
     */
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
            Phi[i] += (phi[j] * Mu[j - 1] * temp);
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
            Epsilon[i] += ((phi[j] * Eta[j - 1] + epsilon[j]) * temp);
        }
    }

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
    }

}

/*void RRMIS::calBeta() {
    vector<double> mu(this->n, 0);
    vector<double> eta(this->n, 0);
    for (int i = 1; i < this->n; ++i) {
        mu[i] = (servers[i - 1].getW() + servers[i - 1].getG() * theta) /
                (servers[i].getW() + servers[i].getG() * theta);
    }
    for (int i = 1; i < this->n; ++i) {
        eta[i] = (servers[i - 1].getS() - servers[i].getS() + servers[i - 1].getO() - servers[i].getO()) /
                (servers[i].getW() + servers[i].getG() * theta);
    }

    // cal Mu (beta)
    for (int i = 1; i < this->n; i++) {
        Mu[i] = 1.0;
        for (int j = 1; j <= i; j++) {
            Mu[i] *= mu[j];
        }
    }
    // cal Eta (beta)
    for (int i = 1; i < this->n; ++i) {
        Eta[i] = 0;
        for (int j = 1; j <= i; ++j) {
            double temp = 1.0;
            for (int k = j + 1; k <= i; ++k) {
                temp *= mu[k];
            }
            Eta[i] += (eta[j] * temp);
        }
    }

    double sum_2_n_eta = 0, sum_2_n_mu = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_mu += Mu[i];
        sum_2_n_eta += Eta[i];
    }

    for (int i = 0 ; i < this->n; ++i) {
        if (i == 0) {
            beta[i] = (1.0 - sum_2_n_eta / this->V) / (1.0 + sum_2_n_mu);
            continue;
        }
        beta[i] = Mu[i] * beta[0] + (Eta[i] / this->V);
    }
}*/

void RRMIS::calBeta() {
    Mu[0] = 1.0;
    for (int i = 1; i < this->n; ++i) {
        Mu[i] = (servers[0].getW() + servers[0].getG() * theta) / (servers[i].getW() + servers[i].getG() * theta);
    }
    for (int i = 1; i < this->n; ++i) {
        Eta[i] = (servers[0].getS() - servers[i].getS() + servers[0].getO() - servers[i].getO()) /
                 (servers[i].getW() + servers[i].getG() * theta);
    }

    double sum_2_n_eta = 0, sum_2_n_mu = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_mu += Mu[i];
        sum_2_n_eta += Eta[i];
    }

    for (int i = 0 ; i < this->n; ++i) {
        if (i == 0) {
            beta[i] = (1.0 - sum_2_n_eta / this->V) / (1.0 + sum_2_n_mu);
            continue;
        }
        beta[i] = Mu[i] * beta[0] + (Eta[i] / this->V);
    }

}

void RRMIS::calGamma() {
/**
     * cal lambda & psi & rho (gamma)
     */
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
        psi[i] = (servers[i].getG() * this->theta) / servers[i].getW();
    }
    for (int i = 1; i < this->n; ++i) {
        Psi[i] = 0;
        for (int j = 1; j <= i; j++) {
            double value = 1.0;
            for (int k = j + 1; k <= i; k++) {
                value *= lambda[k];
            }
            Psi[i] += ((psi[j] * Mu[j]) * value);
        }
    }
    // cal Rho
    for (int i = 1; i < this->n; ++i) {
        rho[i] = (servers[i - 1].getS() - servers[i].getO() - servers[i].getS()) / servers[i].getW();
    }
    for (int i = 1; i < this->n; ++i) {
        P[i] = 0;
        for (int j = 1; j <= i; j++) {
            double value = 1.0;
            for (int k = j + 1; k <= i; k++) {
                value *= lambda[k];
            }
            P[i] += ((rho[j] - psi[j] * Eta[j]) * value);
        }
    }

    // cal gamma
    double sum_2_n_psi = 0, sum_2_n_p = 0, sum_2_n_lambda = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_psi += Psi[i];
        sum_2_n_p += P[i];
        sum_2_n_lambda += Lambda[i];
    }
    for (int i = 0; i < this->n; ++i) {
        if (i == 0) {
            gamma[i] = (1.0 + beta[0] * sum_2_n_psi - (sum_2_n_p / this->V)) / (1.0 + sum_2_n_lambda);
            continue;
        }
        gamma[i] = Lambda[i] * gamma[0] - Psi[i] * beta[0] + P[i] / this->V;
    }

}

void RRMIS::calOptimalM() {
    this->m = 10;
}

void RRMIS::getOptimalModel() {
    this->optimalTime = 0;
    this->optimalTime += (servers[0].getO() + servers[0].getS() + servers[0].getW() * alpha[0] * this->V +
            servers[0].getG() * alpha[0] * this->V * theta);
    this->optimalTime += ((servers[0].getO() + servers[0].getS() + servers[0].getW() * beta[0] * this->V +
            servers[0].getG() * beta[0] * this->V * theta) * (this->m - 2));
    this->optimalTime += (servers[0].getO() + servers[0].getS() + servers[0].getW() * gamma[0] * this->V);
    for (int i = 0; i < this->n; ++i) {
        this->optimalTime += (servers[i].getO() + servers[i].getG() * gamma[i] * this->V * this->theta);
    }

    /**
     * update usingTime
     */
    for (int i = 0; i < this->n; i++) {
        int id = servers[i].getId();
        usingTime[id] = (alpha[i] * V * servers[i].getW() + servers[i].getS() +
                (m - 2) * (beta[i] * V * servers[i].getW() + servers[i].getS()) +
                    gamma[i] * V * servers[i].getW() + servers[i].getS());
        usingTime[id] += (alpha[i] + (m - 2) * beta[i] + gamma[i]) * servers[i].getG() * theta * V;
    }
}

void RRMIS::getDataFromFile() {
    FILE *fpo, *fps, *fpg, *fpw, *totalW;
    double valueO[this->n], valueS[this->n], valueG[this->n], valueW[this->n];

    fpo = fopen("../o.txt", "r");
    fps = fopen("../s.txt", "r");
    fpg = fopen("../g.txt", "r");
    fpw = fopen("../w.txt", "r");
    totalW  = fopen("../WTotal.txt", "r");

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

void RRMIS::printResult() {
    FILE * fpResult;
    fpResult = fopen((outputName + ".txt").c_str(), "a+");

    if (fpResult == nullptr) {
        printf("The file results.txt can not be opened:\n");
        exit(-1);
    }

    fprintf(fpResult, "time: %lf\n", this->optimalTime);
    fprintf(fpResult, "installment: %d\n", this->m);
    fprintf(fpResult, "workload: %lf\n", this->WWithoutError);

    fprintf(fpResult, "using time: (server id : using time)\n");
    usingRate = 0.0;
    for (int i = 0; i < this->n; ++i) {
        fprintf(fpResult, "%d : %lf\n", i, usingTime[i]);
        usingRate += ((double)usingTime[i] / ((double)(this->optimalTime) * this->n));
    }
    fprintf(fpResult, "using rate: %lf\n", usingRate);

    fprintf(fpResult, "\n\n");

    fclose(fpResult);
}

void RRMIS::error(vector<int> &errorPlace, int errorInstallment) {
    outputName = outputName + '_' + to_string((int)errorPlace.size()) +
                 "_installment_" + to_string(errorInstallment) + "_ server";
    for (auto i : errorPlace) {
        outputName = outputName + '_' + to_string(i);
    }

    // cal startTime
    if (errorInstallment == 1) {
        this->startTime = (servers[0].getS() + servers[0].getO()) * 2 +
                servers[0].getW() * (alpha[0] + beta[0]) * this->V * (1 + theta);
    } else if (errorInstallment == this->m || errorInstallment == this->m - 1) {
        // more
    } else {
        this->startTime = servers[0].getS() + servers[0].getO() + servers[0].getW() * alpha[0] * this->V * (1 + theta) +
                (servers[0].getS() + servers[0].getO() + servers[0].getW() * beta[0] * this->V * (1 + theta)) * errorInstallment;
    }

    // release time
    vector<double> busyTime(this->n, 0);
    double preTime = -(servers[0].getO() + servers[0].getG() * beta[0] * V * theta);
    for (int i = 0; i < this->n; ++i) {
        preTime += (servers[i].getO());
        busyTime[i] = preTime;
        preTime += (servers[i].getO() + servers[i].getG() * V * beta[i] + servers[i].getG() * beta[i - 1] * V * theta);
    }

    /*for (int i = 0; i < n; ++i) {
        auto iter = find(errorPlace.begin(), errorPlace.end(), i + 1);
        if (iter == errorPlace.end()) {
            usingTime[i] = (servers[i].getS() + alpha[i] * V * servers[i].getW() +
                            errorInstallment * (servers[i - 1].getS() + beta[i] * V * servers[i].getW()));
        } else {
            usingTime[i] = (servers[i].getS() + alpha[i] * V * servers[i - 1].getW() +
                            errorInstallment * (servers[i].getS() + beta[i] * V * servers[i].getW()));
        }
    }*/

    // cal leftServer
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

    // cal leftW
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

    this->W = this->leftW;
//    beforeInstallment = errorInstallment + 1;
    if (errorInstallment == this->m || errorInstallment == this->m - 1)     // more
        this->startTime = optimalTime;

    old_V = V;
    old_beta = beta;

    initValue();
    getOptimalModel();

    this->optimalTime += this->startTime;
//    this->m += beforeInstallment;

    // cal startTime
    position = 0;
    preTime = servers[0].getO();
    vector<double> reTime(serversNumberWithoutError, 0);
    for (int i = 0; i < this->n; ++i) {
        if (i != 0) {
            preTime += (servers[i].getO() + servers[i - 1].getG() * alpha[i - 1] * V);
        }
        reTime[servers[i].getId()] = preTime;
    }

    vector<double> codeTimeGap(serversNumberWithoutError, 0);
    for (int i = 0; i < serversNumberWithoutError; ++i) {
        if (reTime[i] != 0)
            codeTimeGap[i] = busyTime[i] - reTime[i];
    }

    // math
   /* vector<double> mathTimeGap(serversNumberWithoutError, 0);
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
    }*/
//    cout << "code: " << codeTimeGap[serversNumberWithoutError - 1] <<
//        "math: " << mathTimeGap[serversNumberWithoutError - 1] << endl;

//    timeGap = *max_element(codeTimeGap.begin(), codeTimeGap.end());
//    this->optimalTime += timeGap;
}


void RRMIS::setW(double value) {
    this->W = value;
    this->WWithoutError = value;
}

int RRMIS::getM() {
    return this->m;
}

void RRMIS::theLastInstallmentGap() {
    vector<double> freeTime(this->n, 0);
    vector<double> startTime(this->n, 0);
    vector<double> timeGap(this->n, 0);

    double preTime = -(servers[0].getO() + servers[0].getG() * beta[0] * V * theta);
    for (int i = 0; i < this->n; ++i) {
        preTime += servers[i].getO();
        freeTime[i] = preTime;
        preTime += (servers[i].getO() + (1 + theta) * servers[i].getG() * V * beta[i]);
    }

    preTime = -(servers[0].getO() + servers[0].getG() * beta[0] * V * theta);
    for (int i = 0; i < this->n; ++i) {
        preTime += (2 * servers[i].getO() + servers[i].getG() * beta[i] * V * theta);
        startTime[i] = preTime;
        preTime += servers[i].getG() * V * gamma[i];
    }

    for (int i = 0; i < this->n; ++i) {
        timeGap[i] = startTime[i] - freeTime[i] - servers[i].getO() -
                servers[i].getG() * beta[i] * V * theta;
        cout << i << " : " << timeGap[i] << endl;
    }
}

int RRMIS::getBestM() {
    double _time = INT_MAX, ans = 0;
    for (int i = 0; i < 100; ++i) {
        this->m = i;
        this->V = this->W / this->m;

        /**
          * cal alpha & beta & gamma
          */
        calBeta();
        calAlpha();
        calGamma();

        getOptimalModel();

        if (this->optimalTime < _time) {
            ans = this->m;
            _time = this->optimalTime;
        }
    }

    return ans;
}

double RRMIS::getOptimalTime() const {
    return this->optimalTime;
}

double RRMIS::getUsingRate() const {
    return this->usingRate;
}

void RRMIS::setM(int value) {
    this->m = value;
}