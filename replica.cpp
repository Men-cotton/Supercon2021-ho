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
double simulator(const int C[][N_GROUP], const double I_PROB[]) {
    double S[2][N_GROUP] = {}, I[2][N_GROUP] = {};
    // R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0;
    I[0][0] = 1.0;

    for (int t = 0; t < sc21::T; t++) {
        std::fill(S[1], S[1] + N_GROUP, 0.0);
        std::fill(I[1], I[1] + N_GROUP, 0.0);

        for (int i = 0; i < N_GROUP; i++) {
            double sum = 0;
            for (int j = 0; j < N_GROUP; j++) {
                sum += C[i][j] * I[0][j];
            }
            sum *= sc21::BETA2 * S[0][i];
            S[1][i] = S[0][i] - sc21::BETA * S[0][i] * I[0][i] - sum;
            I[1][i] = I[0][i] + sc21::BETA * S[0][i] * I[0][i] + sum -
                      sc21::GAMMA * I[0][i];
            // R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
        std::swap(S[0], S[1]);
        std::swap(I[0], I[1]);
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[0][i] - I_PROB[i]) * (I[0][i] - I_PROB[i]);
    }

    return loss;
}

void modify(int C[][N_GROUP], const int change) {
    if (change <= 0) return;

    // int one[sc21::N_LINK][2];
    // int cnt = 0;
    // for (int i = 0; i < N_GROUP; i++) {
    //     for (int j = i + 1; j < N_GROUP; j++) {
    //         if (C[i][j] == 1) {
    //             one[cnt][0] = i;
    //             one[cnt][1] = j;
    //             cnt++;
    //         }
    //     }
    // }

    for (int i = 0; i < change; i++) {
        int a, b;
        while (true) {
            // int one_index = xor128() % sc21::N_LINK;
            // a = one[one_index][0];
            // b = one[one_index][1];

            // if (a != -1) {
            //     one[one_index][0] = -1;
            //     break;
            // }
            a = xor128() % N_GROUP;
            b = xor128() % N_GROUP;

            if (a < b && C[a][b] == 1) {
                break;
            }
        }

        int c, d;
        while (true) {
            c = xor128() % N_GROUP;
            d = xor128() % N_GROUP;

            if (c < d && C[c][d] == 0) {
                break;
            }
        }
        C[a][b] = 0;
        C[b][a] = 0;
        C[c][d] = 1;
        C[d][c] = 1;
    }
}

// 焼きなまし法
double sa(int C[][N_GROUP], double s_temp) {
    // double min = 100000.0;
    // auto start_t = std::chrono::system_clock::now();
    // double TIME_LIMIT = 0.05;
    double pre_score = simulator(C, sc21::I_PROB);
    // int epoch = 1;

    for (int _ = 0; _ < 100; _++) {
        // auto now_t = std::chrono::system_clock::now();
        // double e =
        //     std::chrono::duration_cast<std::chrono::seconds>(now_t - start_t)
        //         .count();
        // if (e > TIME_LIMIT) break;

        int new_state[N_GROUP][N_GROUP];

        for (int i = 0; i < N_GROUP; i++) {
            for (int j = 0; j < N_GROUP; j++) {
                new_state[i][j] = C[i][j];
            }
        }

        //だんだん変更量を減少させるように
        int change = 2;

        modify(new_state, change);

        double new_score = simulator(new_state, sc21::I_PROB);
        // min = std::min(min, new_score);

        double prob = exp((pre_score - new_score) / s_temp);

        // std::cout<<pre_score<<" "<<new_score<<" "<<temp<<" "<<prob<<"\n";

        double p = xor128() % 1280000;
        p /= 1280000;

        if (prob > p) {
            for (int i = 0; i < N_GROUP; i++) {
                for (int j = 0; j < N_GROUP; j++) {
                    C[i][j] = new_state[i][j];
                }
            }
            // printf("%lf\n", new_score);
            pre_score = new_score;
        }
    }
    // printf("%lf\n", pre_score);
    return pre_score;
}

void init_c(int C[][N_GROUP]) {
    std::vector<std::pair<double, int>> v;

    //初期化
    for (int i = 0; i < N_GROUP; i++) {
        //ここはsc21::I_PROB[i],iでもいいと思います
        v.emplace_back(sc21::I_PROB[i] / sc21::N[i], i);
        for (int j = 0; j < N_GROUP; j++) {
            C[i][j] = 0;
        }
    }

    std::sort(v.begin(), v.end());

    for (int i = 0; i < N_GROUP - 1; i++) {
        C[v[i].second][v[i + 1].second] = 1;
        C[v[i + 1].second][v[i].second] = 1;
    }
}

void initialize(int L[][N_GROUP][N_GROUP], double T[], int len) {
    double s[N_GROUP];
    memcpy(s, sc21::I_PROB, sizeof(sc21::I_PROB));
    std::sort(s, s + N_GROUP);
    double sum = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        sum += s[i];
    }

    //初期化

    for (int k = 0; k < len; k++) {
        init_c(L[k]);
    }
    std::cout << simulator(L[0], sc21::I_PROB) << "\n";

    for (int i = 0; i < len; i++) {
        T[i] = i * 1.5 + 0.001;
    }
}

int main() {
    sc21::SC_input();

    const double k = 0.1;
    const double eps = 1e-7;

    auto start_t = std::chrono::system_clock::now();
    double TIME_LIMIT = 90;
    const int L_num = 48;
    int L[L_num][N_GROUP][N_GROUP];
    double tmp[L_num], score[L_num];
    initialize(L, tmp, L_num);

    while (true) {
        // for (int i = 0; i < L_num; i++) {
        //     if (i != L_num - 1)
        //         printf("%lf ", tmp[i]);
        //     else
        //         printf("%lf\n", tmp[i]);
        // }

        auto now_t = std::chrono::system_clock::now();
        double e =
            std::chrono::duration_cast<std::chrono::seconds>(now_t - start_t)
                .count();
        if (e > TIME_LIMIT) break;

#pragma omp parallel for
        for (int i = 0; i < L_num; i++) {
            score[i] = sa(L[i], tmp[i]);
        }

        double min = score[0];
        for (int i = 0; i < L_num; i++) {
            min = std::min(min, score[i]);
        }

        printf("score: %lf\n", min);

        for (int i = 0; i < L_num - 1; i++) {
            double p =
                std::min(1.0, std::exp(((double)score[i] - score[i + 1]) *
                                       (1.0 / (k * tmp[i] + eps) -
                                        1.0 / (k * tmp[i + 1] + eps))));
            double q = (double)xor128() / INT32_MAX;
            // printf("exp: %lf\n", tmp[i]);
            // printf("p: %le\n", p);
            if (q > p) {
                std::swap(tmp[i], tmp[i + 1]);
            }
        }
    }

    double min = 1e9;
    int ind = -1;
    for (int i = 0; i < L_num; i++) {
        if (min > simulator(L[i], sc21::I_PROB)) {
            ind = i;
            min = simulator(L[i], sc21::I_PROB);
        }
    }

    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            sc21::C[i][j] = L[ind][i][j];
        }
    }

    //提出時にはここを消せ！！！！！！！！！！！
    std::cout << simulator(sc21::C, sc21::I_PROB) << "\n";

    sc21::SC_output();
}
