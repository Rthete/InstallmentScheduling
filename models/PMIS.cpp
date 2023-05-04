#include "PMIS.h"

PMIS::PMIS(int valueN, double valueTheta) :
        n(valueN),
        theta(valueTheta),
        usingRate(0),
        servers(vector<Server>(MAX_N)),

        // error
        numberWithoutError(valueN),

        // alpha & beta
        alpha(vector<double>(MAX_N, 0)),
        beta(vector<double>(MAX_N, 0)),

        // time
        usingTimes(vector<double>(MAX_N, 0)),

        // print
        outputName("PMIS_printResult") {}

// init value
void PMIS::setW(double value) {
    this->W = value;
    this->WWithoutError = value;
}

void PMIS::setM(int value) {
    this->m = value;
}

/**
 * assign model
 */
void PMIS::initValue() {
    this->V = this->W / this->m;

    initAlpha();
    initBeta();

    double count_alpha = 0, count_beta = 0;
    for (int i = 0; i < this->n; ++i) {
        count_alpha += alpha[i];
        count_beta += beta[i];
        if (alpha[i] < 0 || beta[i] < 0) {
            cout << this->W << " error " << "alpha - " << alpha[i] << " beta - " << beta[i] << endl;
        }
    }
    if ((int)((count_alpha + 0.000005) * 100000) != 100000)
        cout << this->W << " count: alpha - " << count_alpha << endl;
    if ((int)((count_beta + 0.000005) * 100000)  != 100000)
        cout << this->W << " count: beta - " << count_beta << endl;
}

void PMIS::initAlpha() {
    /**
     * cal alpha
     * PMIS
     */
    vector<double> Delta(this->n, 0);
    vector<double> Phi(this->n, 0);
    // cal Delta & Phi
    Delta[0] = 1;
    for (int i = 1; i < this->n; ++i) {
        Delta[i] = servers[0].getW() / servers[i].getW();
        Phi[i] = (servers[0].getS() - servers[i].getS()) / servers[i].getW();
    }

    // cal alpha
    double sum_2_n_Delta = 0, sum_2_n_Phi = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_Delta += Delta[i];
        sum_2_n_Phi += Phi[i];
    }
    for (int i = 0; i < this->n; ++i) {
        if (i == 0) {
            alpha[i] = (1.0 - sum_2_n_Phi / this->V) / (1.0 + sum_2_n_Delta);
            continue;
        }
        alpha[i] = Phi[i] / this->V + Delta[i] * alpha[0];
    }
}

void PMIS::initBeta() {
    /**
     * cal beta
     * PMIS
     */
    vector<double> delta(this->n, 0);
    vector<double> epsilon(this->n, 0);
    // cal delta & epsilon
    for (int i = 1; i < this->n; ++i) {
        delta[i] = (servers[i - 1].getS() - servers[i].getO() - servers[i].getS()) / servers[i].getW();
        epsilon[i] = (servers[i - 1].getW() - servers[i - 1].getG()) / servers[i].getW();
    }

    vector<double> Epsilon(this->n, 0);
    vector<double> Gamma(this->n, 0);
    // cal eta & gamma
    for (int i = 1; i < this->n; ++i) {
        Epsilon[i] = 1.0;
        for (int j = 1; j <= i; ++j) {
            Epsilon[i] *= epsilon[j];
        }
    }
    for (int i = 1; i < this->n; ++i) {
        Gamma[i] = 0;
        for (int j = 1; j <= i; j++) {
            double temp = 1.0;
            for (int k = j + 1; k <= i; k++) {
                temp *= epsilon[k];
            }
            Gamma[i] += (delta[j] * temp);
        }
    }

    // cal alpha
    double sum_2_n_Gamma = 0, sum_2_n_Epsilon = 0;
    for (int i = 1; i < this->n; ++i) {
        sum_2_n_Gamma += Gamma[i];
        sum_2_n_Epsilon += Epsilon[i];
    }
    for (int i = 0; i < this->n; ++i) {
        if (i == 0) {
            beta[i] = (1.0 - sum_2_n_Gamma / this->V) / (1.0 + sum_2_n_Epsilon);
            continue;
        }
        beta[i] = Gamma[i] / this->V + Epsilon[i] * beta[0];
    }
}

void PMIS::getDataFromFile() {
    FILE *fpo, *fps, *fpg, *fpw, *totalW;
    double valueO[this->n], valueS[this->n], valueG[this->n], valueW[this->n];

    fpo = fopen("../data/o.txt", "r");
    fps = fopen("../data/s.txt", "r");
    fpg = fopen("../data/g.txt", "r");
    fpw = fopen("../data/w.txt", "r");
    totalW  = fopen("../data/WTotal.txt", "r");

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

double PMIS::getTime() {
    double time = (this->m - 1) * (servers[0].getW() * alpha[0] * this->V + servers[0].getS()) +
                  servers[0].getW() * beta[0] * this->V + servers[0].getS();

    return time;
}

void PMIS::getOptimalModel() {
    this->optimalTime = getTime();

    // result return
    double result_return = 0.0, assign_weight;
    for (int i = 0; i < this->n; ++i) {
        assign_weight = (this->m - 1) * this->V * alpha[i] + this->V * beta[i];
        result_return += (servers[i].getO() + servers[i].getG() * this->theta * assign_weight);
        usingTimes[servers[i].getId()] += (servers[i].getG() * this->theta * assign_weight);
    }

    this->optimalTime += result_return;

    // restore time
    for (int i = 0; i < this->n; ++i) {
        usingTimes[servers[i].getId()] += ((servers[i].getS() + servers[i].getW() * alpha[i] * this->V) * (this->m - 1) +
                servers[i].getS() + servers[i].getW() * beta[i] * this->V);
    }
}

void PMIS::error(vector<int> &errorPlace, int errorInstallment) {
    // cal leftW
    this->leftW = 0.0;
    outputName += "_errorNum_" + to_string((int)errorPlace.size()) + "_installment_" + to_string(errorInstallment) + "_server";
    for (auto & i : errorPlace) {
        leftW += (this->V * alpha[i - 1] * (this->m - 1) + this->V * beta[i - 1]);
        outputName += '_' + to_string(i);
        usingTimes[i - 1] = 0;

    }
    this->W = this->leftW;

    // cal leftServer
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

    // cal best installment
    beforeInstallment = this->m;
    int errorBestInstallment = getBestInstallment();
    setM(errorBestInstallment);
    initValue();

    // cal best time
    this->startTime = this->optimalTime;
    getOptimalModel();
    this->optimalTime += this->startTime;

    this->m += beforeInstallment;
}

void PMIS::printResult() {
    FILE * fpResult;
    fpResult = fopen(("../output/" + outputName + ".txt").c_str(), "a+");

    if (fpResult == nullptr) {
        printf("The file result.txt can not be opened:\n");
        exit(-1);
    }

    fprintf(fpResult, "time: %lf : %lf + %lf\n", this->optimalTime, this->startTime, this->optimalTime - this->startTime);
    fprintf(fpResult, "installment: %d : %d + %d\n", this->m, this->beforeInstallment, this->m - this->beforeInstallment);
    fprintf(fpResult, "workload: %lf : %lf + %lf\n", this->WWithoutError, this->WWithoutError - this->W, this->W);

    // print usingTime
    fprintf(fpResult, "using time: (server id : using time)\n");
    this->usingRate = 0.0;
    for (int i = 0; i < this->numberWithoutError; ++i) {
        fprintf(fpResult, "%d : %lf\n", i, usingTimes[i]);
        usingRate += usingTimes[i] / (this->optimalTime * this->numberWithoutError);
    }
    fprintf(fpResult, "using rate: %lf\n", usingRate);

    fprintf(fpResult, "\n\n");

    fclose(fpResult);
}

int PMIS::getBestInstallment() {
    double best_time = INT_MAX;
    int installment = 0;
    for (int i = 2; i < 100; ++i) {
        setM(i);

        this->V = this->W / this->m;

        initAlpha();
        initBeta();

        bool flag = false;
        double count_alpha = 0, count_beta = 0;
        for (int j = 0; j < this->n; ++j) {
            count_alpha += alpha[j];
            count_beta += beta[j];
            if (alpha[j] < 0 || beta[j] < 0) {
                flag = true;
            }
        }
        if (flag) continue;
        if ((int)(((count_beta + 0.000005) * 100000) / 100000) != 1) continue;
        if ((int)(((count_alpha + 0.000005) * 100000) / 100000) != 1) continue;

        double time = getTime();
        if (best_time > time) {
            installment = i;
            best_time = time;
        }
    }

    return installment;
}

double PMIS::getOptimalTime() {
    return this->optimalTime;
}

double PMIS::getUsingRate() {
    return this->usingRate;
}
