//
// Created by kqzhang on 2022/5/13.
//

#include "SIS.h"
SIS::SIS(int valueN, double valueTheta) :
    n(valueN),
    theta(valueTheta),
    usingRate(0),
    servers(vector<Server>(MAX_N)),

    // error
    numberWithoutError(valueN),

    // alpha
    alpha(vector<double>(MAX_N, 0)),

    // time
    usingTimes(vector<double>(MAX_N, 0)),

    // print
    outputName("SIS") {
}

// init value
void SIS::setW(double value) {
    this->W = value;
    this->WWithoutError = value;
}

/**
 * assign model
 */
void SIS::initValue() {
    /**
     * cal eta & gamma (alpha)
     * SIS
     */
    vector<double> mu(this->n, 0);
    vector<double> lambda(this->n, 0);
    // cal mu & lambda
    for (int i = 1; i < this->n; ++i) {
        mu[i] = (servers[i - 1].getS() - servers[i].getO() - servers[i].getS()) / servers[i].getW();
        lambda[i] = (servers[i - 1].getW() - servers[i - 1].getG()) / servers[i].getW();
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
            cout << this->W << " error " << "alpha - " << alpha[i] << endl;
        }
    }

    if ((int)((count_alpha + 0.000005) * 100000) != 100000)
        cout << this->W << " count: alpha - " << count_alpha << endl;
}

void SIS::getDataFromFile(string data_path) {
    FILE *fpo, *fps, *fpg, *fpw, *totalW;
    double valueO[this->n], valueS[this->n], valueG[this->n], valueW[this->n];;

    fpo = fopen((data_path + "o.txt").c_str(), "r");
    fps = fopen((data_path + "s.txt").c_str(), "r");
    fpg = fopen((data_path + "g.txt").c_str(), "r");
    fpw = fopen((data_path + "w.txt").c_str(), "r");
    totalW  = fopen((data_path + "WTotal.txt").c_str(), "r");

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

void SIS::getOptimalModel() {
    this->m += 1;

    this->optimalTime = servers[0].getO() + servers[0].getS() + servers[0].getW() * alpha[0] * this->W;

    // result return
    double result_return = 0.0;
    for (int i = 0; i < this->n; ++i) {
        result_return += (servers[i].getO() + servers[i].getG() * this->theta * alpha[i] * this->W);
    }

    this->optimalTime += result_return;

    for (int i = 0; i < this->n; ++i) {
        usingTimes[servers[i].getId()] += (servers[i].getS() + servers[i].getW() * alpha[i] * this->W +
                servers[i].getG() * this->theta * alpha[i] * this->W);
    }
}

void SIS::error(vector<int> &errorPlace) {
    this->leftW = 0.0;
    outputName += "_errorNum_" + to_string((int)errorPlace.size()) + "_server";
    for (auto & i : errorPlace) {
        leftW += this->W * alpha[i - 1];
        outputName += '_' + to_string(i);
    }
    this->W = this->leftW;

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
    this->n -= (int)errorPlace.size();

    this->startTime = this->optimalTime;
    initValue();
    getOptimalModel();
    this->optimalTime += this->startTime;

    // restore time
    for (auto & i : errorPlace) {
        usingTimes[i - 1] = 0;
    }
}

void SIS::printResult() {
    FILE * fpResult;
    fpResult = fopen((outputName + ".txt").c_str(), "a+");

    if (fpResult == nullptr) {
        printf("The file result.txt can not be opened:\n");
        exit(-1);
    }

    fprintf(fpResult, "time: %lf : %lf + %lf\n", this->optimalTime, this->startTime, this->optimalTime - this->startTime);
    fprintf(fpResult, "installment: %d\n", this->m);
    fprintf(fpResult, "workload: %lf : %lf + %lf\n", this->WWithoutError, this->WWithoutError - this->W, this->W);

    // print usingTime
    fprintf(fpResult, "using time: (server id : using time)\n");
    this->usingRate = 0.0;
    for (int i = 0; i < this->numberWithoutError; ++i) {
        fprintf(fpResult, "%d : %lf\n", i, usingTimes[i]);
        usingRate += (usingTimes[i] / (this->optimalTime * this->numberWithoutError));
    }
    fprintf(fpResult, "using rate: %lf\n", usingRate);

    fprintf(fpResult, "\n\n");

    fclose(fpResult);
}

void SIS::calUsingRate() {
    this->usingRate = 0.0;
    for (int i = 0; i < this->numberWithoutError; ++i) {
        usingRate += (usingTimes[i] / (this->optimalTime * this->numberWithoutError));
    }
}

double SIS::getOptimalTime() {
    return this->optimalTime;
}

double SIS::getUsingRate() {
    return this->usingRate;
}
