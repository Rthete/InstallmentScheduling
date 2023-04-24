/*
 * @FilePath: \Installment_Scheduling\models\APMISRR.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-04-22 16:46:39
 * @LastEditTime: 2023-04-24 20:26:25
 */
#include "APMISRR.h"

APMISRR::APMISRR(int valueN, double valueTheta) :
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
        usingTimes(vector<double>(MAX_N, 0)) {}

// init value
void APMISRR::setW(double value) {
    this->W = value;
}

void APMISRR::setM(int value) {
    this->m = value;
}

void APMISRR::setLambda(double value) {
    this->lambda = value;
}

void APMISRR::initValue() {
    // 式(11)
    double temp = 0;
    for (int i = 0; i < this->n; ++i) {
        temp += 1 / servers[i].getW();
    }
    this->P = (this->m - this->lambda) / (m * (m - 1) * temp);

    // 式(12)
    for (int i = 0; i < this->n; ++i) {
        alpha[i] = P / servers[i].getW();
    }


    // 式(15)
    vector<double> a(this->n, 0);
    for (int i = 0; i < this->n; ++i) {
        a[i] = (servers[i].getW() + theta * servers[i].getG()) / servers[i+1].getW();
    }
    vector<double> c(this->n, 0);
    for (int i = 0; i < this->n; ++i) {
        c[i] = -((alpha[i] * theta * servers[i].getG() + alpha[i+1] * servers[i+1].getG()) / servers[i+1].getW());
    }

    // 式(17)
    vector<double> D(this->n, 0);
    for (int i = 2; i < this->n; ++i) {
        for (int j = 1; j < i - 1; ++j) {
            double multi_a = 1;
            for (int k = j + 1; k < i - 1; ++k) {
                multi_a *= a[k];
            }
            D[i] += c[j] * multi_a;
        }
    }

    // 式(20)
    double A = (servers[0].getW() + servers[1].getW() + theta * servers[1].getG()) / servers[0].getW();
    for (int i = 2; i < this->n; ++i) {
        double multi_a = 1;
        for (int j = 1; j < i - 1; ++j) {
            multi_a *= a[j];
        }
        A += multi_a * (1 + (theta * servers[i].getG()) / servers[0].getW());
    }

    // 式(21)
    double B = (alpha[1] * servers[1].getG()) / servers[0].getW();
    for (int i = 2; i < this->n; ++i) {
        B += D[i] * (1 + (theta * servers[i].getG()) / servers[0].getW());
    }

    // 式(19)
    beta[1] = (lambda / m - B) / A;

    // 式(18)
    for (int i = 2; i < this->n; ++i) {
        double multi_a = 1;
        for (int j = 1; j < i - 1; ++j) {
            multi_a *= a[j];
        }
        beta[i] = multi_a * beta[1] + D[i];
    }

    // 式(10)
    beta[0] = alpha[1] * servers[1].getG() + beta[1] * servers[1].getW();
    for (int i = 1; i < this->n; ++i) {
        beta[0] += beta[i] * theta * servers[i].getG();
    }
    beta[0] = beta[0] / servers[0].getW();

    double load_sum = 0;
    for (int i = 0; i < this->n; ++i) {
        // cout << "alpha" << i << "=" << alpha[i] << "\t";
        // cout << "beta" << i << "=" << beta[i] << endl;
        load_sum += alpha[i] * (m - 1);
        load_sum += beta[i];
    }
    // cout << "load_sum = " << load_sum << endl;
}

double APMISRR::getOptimalTime() {
    optimalTime = (m - 1) * P + beta[0] * servers[0].getW();
    cout << "optimalTime = " << optimalTime << endl;
    return optimalTime;
}

/**
 * @brief Read o & s & g & w & WTotal from "/data" and set every server.
 * 
 */
void APMISRR::getDataFromFile() {
    FILE *fpo, *fps, *fpg, *fpw, *totalW;
    double valueO[this->n], valueS[this->n], valueG[this->n], valueW[this->n];

    fpo = fopen("../data/o.txt", "r");
    fps = fopen("../data/s.txt", "r");
    fpg = fopen("../data/APMISRR/g.txt", "r");
    fpw = fopen("../data/APMISRR/w.txt", "r");
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
