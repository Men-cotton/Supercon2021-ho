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

double backprop(int C[][N_GROUP], double I_PROB[], double C_back[][N_GROUP]) {
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

void init_C(int C[][N_GROUP]) {
    std::fill(C[0], C[N_GROUP], 0);
    for (int i = 0; i < sc21::N_LINK; i++) {
        int a, b;
        while (true) {
            a = xor128() % N_GROUP;
            b = xor128() % N_GROUP;

            if (a < b && C[a][b] == 0) {
                break;
            }
        }

        C[a][b] = C[b][a] = 1;
    }
}

int main() {
    sc21::SC_input();

    int C[N_GROUP][N_GROUP];
    init_C(C);

    int q = 100;
    while (q--) {
        double C_back[N_GROUP][N_GROUP];

        double loss = backprop(C, sc21::I_PROB, C_back);
        printf("%le\n", loss);

        std::vector<std::tuple<double, int, int>> should_zero, should_one;

        for (int i = 0; i < N_GROUP; i++) {
            for (int j = 0; j < N_GROUP; j++) {
                if (!(i < j)) {
                    continue;
                }

                double grad = C_back[i][j] + C_back[j][i];
                if (grad < 0 && C[i][j] == 0) {
                    should_one.emplace_back(grad, i, j);
                } else if (grad > 0 && C[i][j] == 1) {
                    should_zero.emplace_back(-grad, i, j);
                }
            }
        }

        std::sort(should_zero.begin(), should_zero.end());
        std::sort(should_one.begin(), should_one.end());

        int swap_num =
            std::min({should_zero.size(), should_one.size(), (size_t)1});
        for (int s = 0; s < swap_num; s++) {
            {
                auto& [grad, i, j] = should_zero[s];
                C[i][j] = C[j][i] = 0;
            }
            {
                auto& [grad, i, j] = should_one[s];
                C[i][j] = C[j][i] = 1;
            }
        }
    }

    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            sc21::C[i][j] = C[i][j];
        }
    }

    sc21::SC_output();
}