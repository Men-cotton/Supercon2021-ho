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
    double S[2][N_GROUP] = {0}, I[2][N_GROUP] = {0};
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

void modify(int C[][N_GROUP], int change) {
    if (change <= 0) return;
    for (int i = 0; i < change; i++) {
        int a = xor128() % N_GROUP, b = xor128() % N_GROUP,
            c = xor128() % N_GROUP, d = xor128() % N_GROUP;
        if (a < b && c < d && C[a][b] == 1 && C[c][d] == 0) {
            C[a][b] = 0;
            C[b][a] = 0;
            C[c][d] = 1;
            C[d][c] = 1;
            break;
        }
    }
}

const unsigned long long int cycle_per_sec = 2800000000;

unsigned long long int getCycle() {
    unsigned int low, high;
    __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
    return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
}

double getTime(unsigned long long int begin_cycle) {
    return (double)(getCycle() - begin_cycle) / cycle_per_sec;
}

// 焼きなまし法
void sa(int C[][N_GROUP]) {
    auto start_cycle = getCycle();
    double s[N_GROUP];
    double min = 100000.0;
    memcpy(s, sc21::I_PROB, sizeof(sc21::I_PROB));
    std::sort(s, s + N_GROUP);
    double sum = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        sum += s[i];
    }

    //初期化
    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            C[i][j] = 0;
        }
    }
    for (int i = 0; i < N_GROUP; i++) {
        C[i][(i + 1) % N_GROUP] = 1;
        C[(i + 1) % N_GROUP][i] = 1;
    }
    int cnt = 0;
    for (int i = 1; i < N_GROUP; i++) {
        if (cnt + N_GROUP >= sc21::N_LINK) {
            break;
        }
        if (s[i] < sum / (N_GROUP * 0.7) && C[0][i] == 0) {
            cnt++;
            C[0][i] = C[i][0] = 1;
        }
    }
    std::cout << cnt << "\n";

    for (int i = 0; i < sc21::N_LINK - N_GROUP - cnt; i++) {
        while (true) {
            int a = xor128() % N_GROUP, b = xor128() % N_GROUP;
            if ((a < b) && C[a][b] == 0) {
                C[a][b] = 1;
                C[b][a] = 1;
                break;
            }
        }
    }

    double start_temp = 50, end_temp = 0;
    double TIME_LIMIT = 90;
    double pre_score = simulator(C, sc21::I_PROB);
    int epoch = 1;

    while (true) {
        double now_t = getTime(start_cycle);
        if (now_t > TIME_LIMIT) break;

        int new_state[N_GROUP][N_GROUP];

        for (int i = 0; i < N_GROUP; i++) {
            for (int j = 0; j < N_GROUP; j++) {
                new_state[i][j] = C[i][j];
            }
        }

        //だんだん変更量を減少させるように
        int change = 1 + 70 - 70 * (now_t) / TIME_LIMIT;

        modify(new_state, change);

        double new_score = simulator(new_state, sc21::I_PROB);
        min = std::min(min, new_score);

        double temp =
            start_temp + (end_temp - start_temp) * (now_t) / TIME_LIMIT;
        double prob = exp((pre_score - new_score) / temp);

        // std::cout<<pre_score<<" "<<new_score<<" "<<temp<<" "<<prob<<"\n";

        long double p = xor128() % 1280000;
        p /= 1280000;

        if (prob > p * (now_t) / TIME_LIMIT) {
            for (int i = 0; i < N_GROUP; i++) {
                for (int j = 0; j < N_GROUP; j++) {
                    C[i][j] = new_state[i][j];
                }
            }
            pre_score = new_score;
        }
        if (epoch % 1000 == 0) {
            std::cout << simulator(sc21::C, sc21::I_PROB) << "\n";
        }
        epoch++;
    }
    std::cout << min << "\n";
}

int main() {
    sc21::SC_input();

    sa(sc21::C);

    //提出時にはここを消せ！！！！！！！！！！！
    std::cout << simulator(sc21::C, sc21::I_PROB) << "\n";

    sc21::SC_output();
}