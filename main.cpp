#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

namespace sc21 {
#include "sc21.h"
}

unsigned int xor128() {
    static unsigned int x = 123456789, y = 362436069, z = 521288629,
                        w = 88675123;
    unsigned int t;
    t = (x ^ (x << 11));
    x = y;
    y = z;
    z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

// CとI_PROBを与えると誤差を返す関数。
double simulator(int C[][N_GROUP], double I_PROB[]) {
    double S[sc21::T + 1][N_GROUP] = {0}, I[sc21::T + 1][N_GROUP] = {0},
                       R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0;
    I[0][0] = 1.0;

    for (int t = 0; t < sc21::T; t++) {
        for (int i = 0; i < N_GROUP; i++) {
            double sum = 0;
            for (int j = 0; j < N_GROUP; j++) {
                sum += C[i][j] * I[t][j];
            }
            sum *= sc21::BETA2 * S[t][i];
            S[t + 1][i] = S[t][i] - sc21::BETA * S[t][i] * I[t][i] - sum;
            I[t + 1][i] = I[t][i] + sc21::BETA * S[t][i] * I[t][i] + sum -
                          sc21::GAMMA * I[t][i];
            R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[sc21::T][i] - I_PROB[i]) * (I[sc21::T][i] - I_PROB[i]);
    }

    return loss;
}

void print(double x[N_GROUP]) {
    for (int i = 0; i < N_GROUP; i++) {
        printf("%le ", x[i]);
    }
    printf("\n");
}

double backprop(double C[][N_GROUP], double I_PROB[],
                double C_back[][N_GROUP]) {
    double S[sc21::T + 1][N_GROUP] = {0}, I[sc21::T + 1][N_GROUP] = {0},
                       R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0;
    I[0][0] = 1.0;

    for (int t = 0; t < sc21::T; t++) {
        for (int i = 0; i < N_GROUP; i++) {
            double sum = 0;
            for (int j = 0; j < N_GROUP; j++) {
                sum += C[i][j] * I[t][j];
            }
            sum *= sc21::BETA2 * S[t][i];
            S[t + 1][i] = S[t][i] - sc21::BETA * S[t][i] * I[t][i] - sum;
            I[t + 1][i] = I[t][i] + sc21::BETA * S[t][i] * I[t][i] + sum -
                          sc21::GAMMA * I[t][i];
            R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }

        // print(S[t + 1]);
        // print(I[t + 1]);
        // print(R[t + 1]);
        // printf("\n");
    }
    // printf("\n");

    double S_back[sc21::T + 1][N_GROUP] = {0},
                            I_back[sc21::T + 1][N_GROUP] = {0},
                            R_back[sc21::T + 1][N_GROUP] = {0};

    std::fill(C_back[0], C_back[N_GROUP], 0.0);

    for (int i = 0; i < N_GROUP; i++) {
        I_back[sc21::T][i] = (I[sc21::T][i] - I_PROB[i]) * 2;
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[sc21::T][i] - I_PROB[i]) * (I[sc21::T][i] - I_PROB[i]);
    }

    for (int t = sc21::T - 1; t >= 0; t--) {
        for (int i = 0; i < N_GROUP; i++) {
            R_back[t][i] += R_back[t + 1][i];
            I_back[t][i] += sc21::GAMMA * R_back[t + 1][i] + I_back[t + 1][i] +
                            sc21::BETA * I_back[t + 1][i] * S[t][i] -
                            sc21::GAMMA * I_back[t + 1][i] -
                            sc21::BETA * S[t][i] * S_back[t + 1][i];
            S_back[t][i] += sc21::BETA * I_back[t + 1][i] * I[t][i] +
                            S_back[t + 1][i] -
                            sc21::BETA * I[t][i] * S_back[t + 1][i];
            double sum_back = I_back[t + 1][i] - S_back[t + 1][i];
            double sum = 0;
            for (int j = 0; j < N_GROUP; j++) {
                sum += C[i][j] * I[t][j];
            }

            S_back[t][i] += sc21::BETA2 * sum * sum_back;
            sum_back *= sc21::BETA2 * S[t][i];

            for (int j = 0; j < N_GROUP; j++) {
                I_back[t][j] += sum_back * C[i][j];
                C_back[i][j] += sum_back * I[t][j];
            }
        }
    }
    
    return loss;
}

void init_C(double C[][N_GROUP]) {
    // std::fill(C[0], C[N_GROUP], 0.0);
    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            if (i <= j) C[i][j] = C[j][i] = 0.5;
            if (j == N_GROUP - 1)
                printf("%le\n", C[i][j]);
            else
                printf("%le ", C[i][j]);
        }
    }
}

/*
class Adam:
    '''
    Adam (http://arxiv.org/abs/1412.6980v8)
    '''
    def __init__(self, lr=0.001, beta1=0.9, beta2=0.999):
        self.lr = lr
        self.beta1 = beta1
        self.beta2 = beta2
        self.iter = 0
        self.m = None
        self.v = None
    def update(self, params, grads):
        if self.m is None:
            self.m, self.v = [], []
            for param in params:
                self.m.append(np.zeros_like(param))
                self.v.append(np.zeros_like(param))
        self.iter += 1
        lr_t = self.lr * np.sqrt(1.0 - self.beta2**self.iter) / (1.0 -
self.beta1**self.iter)
        for i in range(len(params)):
            self.m[i] += (1 - self.beta1) * (grads[i] - self.m[i])
            self.v[i] += (1 - self.beta2) * (grads[i]**2 - self.v[i])
            params[i] -= lr_t * self.m[i] / (np.sqrt(self.v[i]) + 1e-7)
*/

double quant_loss_diff(double x) { return 2 * x * ((2 * x - 3) * x + 1); }

int main() {
    sc21::SC_input();

    double C[N_GROUP][N_GROUP], C_m[N_GROUP][N_GROUP] = {},
                                C_v[N_GROUP][N_GROUP] = {};
    init_C(C);
    const double lr = 1e-3, beta1 = 0.9, beta2 = 0.999;
    // const double quant_lambda

    const int q = 1000;
    for (int iter = 1; iter <= q; iter++) {
        double C_back[N_GROUP][N_GROUP];

        double loss = backprop(C, sc21::I_PROB, C_back);
        printf("%le\n", loss);

        const double lr_t = lr * std::sqrt(1.0 - std::pow(beta2, iter)) /
                            (1.0 - std::pow(beta1, iter));

        for (int i = 0; i < N_GROUP; i++) {
            for (int j = 0; j < N_GROUP; j++) {
                if (i < j) {
                    double grad = C_back[i][j] + C_back[j][i];

                    C_m[i][j] += (1 - beta1) * (grad - C_m[i][j]);
                    C_v[i][j] += (1 - beta2) * (grad * grad - C_v[i][j]);

                    C[i][j] -= lr_t * C_m[i][j] / (std::sqrt(C_v[i][j]) + 1e-7);
                    C[j][i] = C[i][j] = std::max(std::min(C[i][j], 1.0), 0.0);
                }
            }
        }
    }

    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            sc21::C[i][j] = std::round(C[i][j]);
        }
    }

    double loss = simulator(sc21::C, sc21::I_PROB);
    printf("%le\n", loss);

    sc21::SC_output();
}