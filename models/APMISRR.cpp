/*
 * @FilePath: \InstallmentScheduling\models\APMISRR.cpp
 * @Description:  
 * @Author: rthete
 * @Date: 2023-04-22 16:46:39
 * @LastEditTime: 2023-05-09 20:48:01
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

double APMISRR::getAlpha() {
    return alpha[1];
}

double APMISRR::getBeta() {
    return beta[1];
}

void APMISRR::initValue_cost() {
    // 每趟调度任务量
    this->V = (this->m - this->lambda) * this->W / (m * (m - 1));

    // 式(11)
    double temp = 0;
    for (int i = 0; i <= this->n; ++i) {
        temp += 1 / servers[i].getW();
    }
    double temp1 = 0;
    for (int i = 0; i <= this->n; ++i) {
        temp1 += (servers[i].getS() / servers[i].getW());
    }
    this->P = (this->V + temp1) / temp;
    cout << "P = " << P << endl;

    // 式(12)
    double sum = 0;
    for (int i = 0; i <= this->n; ++i) {
        alpha[i] = (P - servers[i].getS())/ (servers[i].getW() * this->V);
        // cout << "alpha[" << i << "] = " << alpha[i] << endl;
        sum += alpha[i];
        // cout << "start cost for computation: " << servers[i].getS() << endl;
        // cout << "compute time in first installment: " << alpha[i] * servers[i].getW() * this->V << endl;
        // cout << "start cost + compute time = " << servers[i].getS() + alpha[i] * servers[i].getW() * this->V << endl; // ==P
    }
    // cout << "sum = " << sum << endl; // alpha[i]之和==1

    temp = 0;
    for (int i = 1; i <= this->n; ++i) {
        temp += servers[i].getO() + alpha[i] * servers[i].getG() * this->V + servers[i].getO() + 
                alpha[i] * theta * servers[i].getG() * this->V;
    }

    // 内部调度间空闲时间
    double d = P - temp;
    cout << "d = " << d << endl;

    this->Vb = lambda / m * this->W;
    // 式(15)
    vector<double> a(this->n + 1, 0);
    for (int i = 0; i < this->n; ++i) {
        a[i] = (servers[i].getW() + theta * servers[i].getG()) / servers[i+1].getW();
    }
    vector<double> c(this->n + 1, 0);
    for (int i = 0; i < this->n; ++i) {
        c[i] = ((- alpha[i] * theta * servers[i].getG() - alpha[i+1] * servers[i+1].getG()) * this->V + 
                servers[i].getS() - servers[i].getO() - servers[i+1].getS() - servers[i+1].getO()) / 
                (servers[i+1].getW() * this->Vb);
    }

    // 式(17)
    vector<double> D(this->n + 1, 0);
    for (int i = 2; i <= this->n; ++i) {
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
    for (int i = 2; i <= this->n; ++i) {
        double multi_a = 1;
        for (int j = 1; j < i - 1; ++j) {
            multi_a *= a[j];
        }
        A += multi_a * (1 + (theta * servers[i].getG()) / servers[0].getW());
    }

    // 式(21)
    double B = (2 * servers[1].getO() + alpha[1] * servers[1].getG() * this->V + servers[1].getS() - servers[0].getS()) / (servers[0].getW() * this->Vb);
    for (int i = 2; i <= this->n; ++i) {
        B += D[i] * (1 + (theta * servers[i].getG()) / servers[0].getW()) + servers[i].getO() / (servers[0].getW() * this->Vb);
    }

    // 式(19)
    beta[1] = (1 - B) / A;

    // 式(18)
    for (int i = 2; i <= this->n; ++i) {
        double multi_a = 1;
        for (int j = 1; j < i - 1; ++j) {
            multi_a *= a[j];
        }
        beta[i] = multi_a * beta[1] + D[i];
    }

    // 式(10)
    beta[0] = servers[1].getO() + alpha[1] * servers[1].getG() * this->V + servers[1].getS() + beta[1] * servers[1].getW() * this->Vb;
    for (int i = 1; i <= this->n; ++i) {
        beta[0] += servers[i].getO() + beta[i] * theta * servers[i].getG() * this->Vb;
    }
    beta[0] = (beta[0] - servers[0].getS()) / (servers[0].getW() * this->Vb);
}

void APMISRR::initValue() {
    // 式(11)
    double temp = 0;
    for (int i = 0; i <= this->n; ++i) {
        temp += 1 / servers[i].getW();
    }
    this->P = (this->m - this->lambda) / (m * (m - 1) * temp);

    // 式(12)
    for (int i = 0; i <= this->n; ++i) {
        alpha[i] = P / servers[i].getW();
    }


    // 式(15)
    vector<double> a(this->n + 1, 0);
    for (int i = 0; i < this->n; ++i) {
        a[i] = (servers[i].getW() + theta * servers[i].getG()) / servers[i+1].getW();
    }
    vector<double> c(this->n + 1, 0);
    for (int i = 0; i < this->n; ++i) {
        c[i] = -((alpha[i] * theta * servers[i].getG() + alpha[i+1] * servers[i+1].getG()) / servers[i+1].getW());
    }

    // 式(17)
    vector<double> D(this->n + 1, 0);
    for (int i = 2; i <= this->n; ++i) {
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
    for (int i = 2; i <= this->n; ++i) {
        double multi_a = 1;
        for (int j = 1; j < i - 1; ++j) {
            multi_a *= a[j];
        }
        A += multi_a * (1 + (theta * servers[i].getG()) / servers[0].getW());
    }

    // 式(21)
    double B = (alpha[1] * servers[1].getG()) / servers[0].getW();
    for (int i = 2; i <= this->n; ++i) {
        B += D[i] * (1 + (theta * servers[i].getG()) / servers[0].getW());
    }

    // 式(19)
    beta[1] = (lambda / m - B) / A;

    // 式(18)
    for (int i = 2; i <= this->n; ++i) {
        double multi_a = 1;
        for (int j = 1; j < i - 1; ++j) {
            multi_a *= a[j];
        }
        beta[i] = multi_a * beta[1] + D[i];
    }

    // 式(10)
    beta[0] = alpha[1] * servers[1].getG() + beta[1] * servers[1].getW();
    for (int i = 1; i <= this->n; ++i) {
        beta[0] += beta[i] * theta * servers[i].getG();
    }
    beta[0] = beta[0] / servers[0].getW();

    
}

double APMISRR::getOptimalTime() {
    // 式(29)
    optimalTime = (m - 1) * P + beta[0] * servers[0].getW();
    cout << "optimalTime = " << optimalTime << endl;
    return optimalTime;
}

double APMISRR::getOptimalTime_cost() {
    // 式(29)
    optimalTime = (m - 1) * P + beta[0] * servers[0].getW() * this->V + servers[0].getO();
    cout << "optimalTime = " << optimalTime << endl;
    return optimalTime;
}

void APMISRR::isSchedulable() {
    int condition_1 = 0, condition_2 = 0, condition_3 = 0;
    
    // 条件1，推论4
    double right;
    right = P / servers[n - 1].getW() * 
            (theta * servers[n - 1].getG() / servers[n - 1].getW() + (1 + theta) * servers[n].getG() / servers[n].getW());
    if (beta[n - 1] >= right) {
        condition_1 = 1;
    }

    // 条件2，式(13)
    double sum = 0;
    for (int i = 1; i <= this->n; ++i) {
        // cout << servers[i].getG() << "\t" << servers[i].getW() << endl;
        sum += (servers[i].getG() / servers[i].getW());
    }
    if (sum <= 1 / (1 + theta)) {
        condition_2 = 1; 
    }
    // cout << "sum = " << sum << ", 1 / (1 + theta) = " << 1 / (1 + theta) << endl;

    if(alpha[1] >= beta[1]) {
        condition_3 = 1;
    }

    if (condition_1 && condition_2 && condition_3) {
        cout << "schedulable" << endl;
    } else {
        cout << "condition_1 = " << condition_1 << ", condition_2 = " << condition_2 << ", condition_3 = " << condition_3 << endl;
        cout << "not schedulable" << endl;
    }
}

int APMISRR::isSchedulable_cost() {

    // 打印alpha[i], beta[i]
    double load_sum = 0;
    for (int i = 0; i <=this->n; ++i) {
        cout << "alpha" << i << "=" << alpha[i] << "\t";
        cout << "beta" << i << "=" << beta[i] << endl;
        // load_sum += alpha[i] * (m - 1) * this->V;
        // load_sum += beta[i] * this->Vb;
    }
    // load_sum为1
    // cout << "load_sum = " << load_sum << endl;

    int condition1 = 1, condition2 = 1, condition3 = 1;

    // 负载分量必须为正数
    for (int i = 0; i <= this->n; ++i) {
        if(!(alpha[i] > 0 && beta[i] > 0)) {
            condition1 = 0;
            cout << "!!!!!!! not schedulable 1 !!!!!!!" << endl;
            break;
        }
    }

    // 引理8
    int num = this->n;
    if(servers[num].getO() + beta[num] * servers[num].getW() * this->V < 
        servers[num].getS() + alpha[num] * theta * servers[num].getG() * this->V + servers[num-1].getS() + 
        beta[num-1] * theta * servers[num-1].getG() * this->V) {
            condition2 = 0;
            cout << "!!!!!!! not schedulable 2 !!!!!!!" << endl;
        }

    // 处理器P_1的最后一趟负载分量是否能在传输窗口放下
    if(alpha[1] * this->V < beta[1] * this->Vb) {
        condition3 = 0;
        cout << "!!!!!!! not schedulable 3 !!!!!!!" << endl;
    }

    if(condition1 && condition2 && condition3) {
        cout << "```````````` schedulable ````````````" << endl;
        return 1;
    } else {
        return 0;
    }
}

/**
 * @brief Read o & s & g & w & WTotal from "/data" and set every server.
 * 
 */
void APMISRR::getDataFromFile() {
    FILE *fpo, *fps, *fpg, *fpw, *totalW;
    double valueO[this->n + 1], valueS[this->n + 1], valueG[this->n + 1], valueW[this->n + 1];

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

    for (int i = 0; i < this->n + 1; i++) {
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
