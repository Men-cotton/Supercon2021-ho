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
double simulator_bool(const bool C[][N_GROUP], const double I_PROB[],
                      double LOSS[]) {
    double S[2][N_GROUP] = {}, I[2][N_GROUP] = {}, SUM[N_GROUP];
    // R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0;
    I[0][0] = 1.0;

    int C_LIST[sc21::N_LINK][2];
    int cnt = 0;
    for (int i = 0; i < N_GROUP; i++) {
        for (int j = i + 1; j < N_GROUP; j++) {
            if (C[i][j]) {
                C_LIST[cnt][0] = i;
                C_LIST[cnt][1] = j;
                cnt++;
            }
        }
    }

    for (int t = 0; t < sc21::T; t++) {
        std::fill(S[1], S[1] + N_GROUP, 0.0);
        std::fill(I[1], I[1] + N_GROUP, 0.0);
        std::fill(SUM, SUM + N_GROUP, 0.0);

        for (int i = 0; i < sc21::N_LINK; i++) {
            auto &[x, y] = C_LIST[i];
            SUM[x] += I[0][y];
            SUM[y] += I[0][x];
        }

        for (int i = 0; i < N_GROUP; i++) {
            double sum = SUM[i] * sc21::BETA2 * S[0][i];
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
        double l = std::abs(I[0][i] - I_PROB[i]);
        LOSS[i] = l;
        loss += l;
    }

    for (int i = 0; i < N_GROUP; i++) {
        int cnt = 1;
        for (int j = 0; j < N_GROUP; j++) {
            cnt += C[i][j];
        }
        LOSS[i] *= cnt;
    }

    return loss;
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

double loss_sum(bool C[][N_GROUP], double loss[N_GROUP]) {
    double sum = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        sum += loss[i];
    }
    return sum;
}

int weight_sample(double weight[N_GROUP], double sum) {
    double s = double(xor128() % 12800000) / 12800000 * sum;
    for (int i = 0; i < N_GROUP; i++) {
        if (s < weight[i]) {
            return i;
        }
        s -= weight[i];
    }
}

void modify(bool C[][N_GROUP], const int change, double loss[N_GROUP],
            double loss_sum) {
    if (change <= 0) return;

    for (int i = 0; i < change; i++) {
        int a, b;
        while (true) {
            a = xor128() % N_GROUP;
            b = weight_sample(loss, loss_sum);

            if (a != b && C[a][b] == 1) {
                break;
            }
        }
        int c;
        bool flag;
        while (true) {
            c = weight_sample(loss, loss_sum);
            flag = xor128() % 2;

            if (c != a) {
                if (flag) {
                    if (C[c][a] == 0) {
                        break;
                    }
                } else {
                    if (C[a][c] == 0) {
                        break;
                    }
                }
            }
        }
        C[a][b] = 0;
        C[b][a] = 0;
        C[c][a] = 1;
        C[a][c] = 1;
    }
}


// 焼きなまし法
double sa(bool C[][N_GROUP], double s_temp, bool best[][N_GROUP],
          double &min_score) {
    // double min = 100000.0;
    // auto start_t = std::chrono::system_clock::now();
    // double TIME_LIMIT = 0.05;
    double pre_loss[N_GROUP], new_loss[N_GROUP];
    double pre_score = simulator_bool(C, sc21::I_PROB, pre_loss);
    double pre_sample_sum = loss_sum(C, pre_loss);

    // int epoch = 1;

    int change = 1 + (int)s_temp / 20;

    for (int epoch = 0; epoch < 100; epoch++) {
        // auto now_t = std::chrono::system_clock::now();
        // double e =
        //     std::chrono::duration_cast<std::chrono::seconds>(now_t - start_t)
        //         .count();
        // if (e > TIME_LIMIT) break;

        bool new_state[N_GROUP][N_GROUP];

        for (int i = 0; i < N_GROUP; i++) {
            for (int j = 0; j < N_GROUP; j++) {
                new_state[i][j] = C[i][j];
            }
        }

        //だんだん変更量を減少させるように

        modify(new_state, change, pre_loss, pre_sample_sum);

        double new_score = simulator_bool(new_state, sc21::I_PROB, new_loss);
        // printf("%lf, %lf\n", pre_score, new_score);
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

            for (int i = 0; i < N_GROUP; i++) {
                pre_loss[i] = new_loss[i];
            }
            pre_sample_sum = loss_sum(C, pre_loss);

            if (min_score > pre_score) {
                min_score = pre_score;
                for (int i = 0; i < N_GROUP; i++) {
                    for (int j = 0; j < N_GROUP; j++) {
                        best[i][j] = C[i][j];
                    }
                }
            }
        }
    }

    // printf("%lf\n", pre_score);
    return pre_score;
}

void init_c(bool C[][N_GROUP], double s[N_GROUP], double sum) {
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
    // std::cout << cnt << "\n";

    for (int i = 0; i < sc21::N_LINK - (N_GROUP - 1) - cnt; i++) {
        while (true) {
            int a = xor128() % N_GROUP, b = xor128() % N_GROUP;
            if ((a < b) && C[a][b] == 0) {
                C[a][b] = 1;
                C[b][a] = 1;
                break;
            }
        }
    }
}

void initialize(bool L[][N_GROUP][N_GROUP], double T[], int len) {
    double s[N_GROUP];
    memcpy(s, sc21::I_PROB, sizeof(sc21::I_PROB));
    std::sort(s, s + N_GROUP);
    double sum = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        sum += s[i];
    }

    //初期化

    for (int k = 0; k < len; k++) {
        init_c(L[k], s, sum);
    }
    std::cout << simulator_bool(L[0], sc21::I_PROB) << "\n";

    double start_temp = 0.01, end_temp = 15.0;

    for (int i = 0; i < len; i++) {
        T[i] = start_temp + (end_temp - start_temp) * (double)i / (double)len;
    }
}

int main() {
    sc21::SC_input();

    const double k = 3.0;
    const double eps = 1e-7;

    auto start_t = std::chrono::system_clock::now();
    double TIME_LIMIT = 295;
    const int L_num = 48;
    bool L[L_num][N_GROUP][N_GROUP];
    bool BEST[L_num][N_GROUP][N_GROUP] = {};
    double tmp[L_num], score[L_num];
    double best_score[L_num];
    initialize(L, tmp, L_num);
    std::fill(best_score, best_score + L_num, 1e9);

    while (true) {
        // for (int i = 0; i < L_num; i++) {
        //     if (i != L_num - 1)
        //         printf("%lf ", tmp[i]);
        //     else
        //         printf("%lf\n", tmp[i]);
        // }

        auto now_t = std::chrono::system_clock::now();
        double e = std::chrono::duration_cast<std::chrono::milliseconds>(
                       now_t - start_t)
                       .count() /
                   1000.0;
        if (e > TIME_LIMIT) break;

        double min = best_score[0];
        for (int i = 0; i < L_num; i++) {
            min = std::min(min, best_score[i]);
        }

        for (int i = 0; i < L_num; i++) {
            best_score[i] = std::min(best_score[i], min + 0.01);
        }

        double temp_scale = (1 - e / TIME_LIMIT);
        // printf("%lf\n", temp_scale);

#pragma omp parallel for
        for (int i = 0; i < L_num; i++) {
            score[i] = sa(L[i], tmp[i] * temp_scale, BEST[i], best_score[i]);
        }

        min = best_score[0];
        for (int i = 0; i < L_num; i++) {
            min = std::min(min, best_score[i]);
        }

        printf("score: %lf\n", min);

        for (int ii = 0; ii < 5; ii++) {
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
    }

    double min = 1e9;
    int ind = -1;
    for (int i = 0; i < L_num; i++) {
        if (min > best_score[i]) {
            ind = i;
            min = best_score[i];
        }
    }

    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            sc21::C[i][j] = BEST[ind][i][j];
        }
    }

    //提出時にはここを消せ！！！！！！！！！！！
    std::cout << simulator(sc21::C, sc21::I_PROB) << "\n";

    sc21::SC_output();
}
