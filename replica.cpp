#include <bits/stdc++.h>
#include <cmath>
#include <omp.h>
#include <cstdio>
#include <cstdlib>

namespace sc21 {

#include "sc21.h"

}

inline unsigned int xor128() {
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 88675123;
    unsigned int t;
    t = (x ^ (x << 11));
    x = y;
    y = z;
    z = w;
    return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

inline double zero_one_random() {
    constexpr unsigned int MOD = 1 << 20;
    return (double) (xor128() & (MOD - 1)) / MOD;
}


// CとI_PROBを与えると誤差を返す関数。
double simulator_sa(const bool C[][N_GROUP], const double I_PROB[], double LOSS[], const bool abs = false) {
    double S[2][N_GROUP] = {}, I[2][N_GROUP] = {}, SUM[N_GROUP];
    // R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0, I[0][0] += 1.0;

    int C_LIST[sc21::N_LINK][2];
    {
        int cnt = 0;
        for (int i = 0; i < N_GROUP; i++) {
            for (int j = i + 1; j < N_GROUP; j++) {
                if (C[i][j]) {
                    C_LIST[cnt][0] = i, C_LIST[cnt][1] = j;
                    cnt++;
                }
            }
        }
    }

    for (int t = 0; t < sc21::T; t++) {
        std::fill(S[1], S[1] + N_GROUP, 0.0);
        std::fill(I[1], I[1] + N_GROUP, 0.0);
        std::fill(SUM, SUM + N_GROUP, 0.0);

        for (int i = 0; i < sc21::N_LINK; i++) {
            auto &[x, y] = C_LIST[i];
            SUM[x] += I[0][y], SUM[y] += I[0][x];
        }

        for (int i = 0; i < N_GROUP; i++) {
            double sum = SUM[i] * sc21::BETA2 * S[0][i];
            S[1][i] = S[0][i] - sc21::BETA * S[0][i] * I[0][i] - sum;
            I[1][i] = I[0][i] + sc21::BETA * S[0][i] * I[0][i] + sum - sc21::GAMMA * I[0][i];
            // R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
        std::swap(S[0], S[1]), std::swap(I[0], I[1]);
    }

    double loss = 0.0;
    if (abs) {
        for (int i = 0; i < N_GROUP; i++) {
            double l = std::abs(I[0][i] - I_PROB[i]);
            LOSS[i] = l;
            loss += l;
        }
    } else {
        for (int i = 0; i < N_GROUP; i++) {
            double l = std::abs(I[0][i] - I_PROB[i]);
            LOSS[i] = l;
            loss += l * l;
        }
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
double simulator(const bool C[][N_GROUP], const double I_PROB[]) {
    double S[2][N_GROUP] = {}, I[2][N_GROUP] = {};
    // R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0, I[0][0] += 1.0;

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
            I[1][i] = I[0][i] + sc21::BETA * S[0][i] * I[0][i] + sum - sc21::GAMMA * I[0][i];
            // R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
        std::swap(S[0], S[1]), std::swap(I[0], I[1]);
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[0][i] - I_PROB[i]) * (I[0][i] - I_PROB[i]);
    }

    return loss;
}

// CとI_PROBを与えると誤差を返す関数。
double simulator_final(const int C[][N_GROUP], const double I_PROB[]) {
    double S[2][N_GROUP] = {}, I[2][N_GROUP] = {};
    // R[sc21::T + 1][N_GROUP] = {0};
    for (int i = 0; i < N_GROUP; i++) {
        S[0][i] = sc21::N[i];
    }
    S[0][0] -= 1.0, I[0][0] += 1.0;

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
            I[1][i] = I[0][i] + sc21::BETA * S[0][i] * I[0][i] + sum - sc21::GAMMA * I[0][i];
            // R[t + 1][i] = R[t][i] + sc21::GAMMA * I[t][i];
        }
        std::swap(S[0], S[1]), std::swap(I[0], I[1]);
    }

    double loss = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        loss += (I[0][i] - I_PROB[i]) * (I[0][i] - I_PROB[i]);
    }

    return loss;
}

inline double loss_sum(const bool C[][N_GROUP], const double loss[N_GROUP]) {
    double sum = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        sum += loss[i];
    }
    return sum;
}

int weight_sample(double weight[N_GROUP], double sum) {
    double s = zero_one_random() * sum;
    for (int i = 0; i < N_GROUP; i++) {
        if (s < weight[i]) {
            return i;
        }
        s -= weight[i];
    }
    return N_GROUP - 1;
}

// 辺を change 本付け替える
void modify(bool C[][N_GROUP], const int change, double loss[N_GROUP], double loss_sum) {
    if (change <= 0) return;

    for (int i = 0; i < change; i++) {
        int base_vertex, pre_vertex;
        while (true) {
            base_vertex = xor128() % N_GROUP;
            pre_vertex = weight_sample(loss, loss_sum);

            if (base_vertex != pre_vertex && C[base_vertex][pre_vertex]) {
                break;
            }
        }

        int new_vertex;
        while (true) {
            new_vertex = weight_sample(loss, loss_sum);

            if (new_vertex != base_vertex) {
                if (!C[new_vertex][base_vertex] && !C[base_vertex][new_vertex]) {
                    break;
                }
            }
        }
        C[base_vertex][pre_vertex] = C[pre_vertex][base_vertex] = false;
        C[base_vertex][new_vertex] = C[new_vertex][base_vertex] = true;
    }
}

double yakinamashi(bool C[][N_GROUP], const double s_temp, bool best[][N_GROUP], double &min_score, const bool abs) {
    // double min = 100000.0;
    // auto start_t = std::chrono::system_clock::now();
    // double TIME_LIMIT = 0.05;
    double pre_loss[N_GROUP], new_loss[N_GROUP];
    double pre_score = simulator_sa(C, sc21::I_PROB, pre_loss, abs);
    double pre_sample_sum = loss_sum(C, pre_loss);

    // int epoch = 1;

    const int change = 1 + (int) s_temp / 20;

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

        double new_score = simulator_sa(new_state, sc21::I_PROB, new_loss, abs);
        // printf("%lf, %lf\n", pre_score, new_score);
        // min = std::min(min, new_score);

        double prob = exp((pre_score - new_score) / s_temp);

        // std::cout<<pre_score<<" "<<new_score<<" "<<temp<<" "<<prob<<"\n";

        if (prob > zero_one_random()) {
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

void init_c(bool C[][N_GROUP], const double s[N_GROUP], const double sum_I) {
    std::vector<std::pair<double, int>> ratio;

    //初期化
    for (int i = 0; i < N_GROUP; i++) {
        //ここはsc21::I_PROB[i],iでもいいと思います
        ratio.emplace_back(sc21::I_PROB[i] / sc21::N[i], i);
        for (int j = 0; j < N_GROUP; j++) {
            C[i][j] = false;
        }
    }

    std::sort(ratio.begin(), ratio.end());

    int cnt = 0;
    for (int i = 0; i < N_GROUP - 1; i++) {
        cnt++;
        int vertex1 = ratio[i].second, vertex2 = ratio[i + 1].second;
        C[vertex1][vertex2] = C[vertex2][vertex1] = true;
    }

    for (int i = 1; i < N_GROUP; i++) {
        if (cnt + 1 >= sc21::N_LINK) { // ランダム辺を一個は作る
            break;
        }
        if (s[i] < sum_I / (N_GROUP * 0.7) && !C[0][i]) {
            cnt++;
            C[0][i] = C[i][0] = true;
        }
    }
    // std::cout << cnt << "\n";

    for (int i = cnt; i < sc21::N_LINK; i++) { // ランダム辺生成
        while (true) {
            int a = xor128() % N_GROUP, b = xor128() % N_GROUP;
            if ((a < b) && !C[a][b]) {
                C[a][b] = C[b][a] = true;
                break;
            }
        }
    }
}

void initialize(bool L[][N_GROUP][N_GROUP], double T[], const int THREAD_NUM) {
    double s[N_GROUP];
    memcpy(s, sc21::I_PROB, sizeof(sc21::I_PROB));
    std::sort(s, s + N_GROUP);
    double I_sum = 0.0;
    for (int i = 0; i < N_GROUP; i++) {
        I_sum += s[i];
    }

    // グラフ初期化
    for (int thread = 0; thread < THREAD_NUM; thread++) {
        init_c(L[thread], s, I_sum);
    }
    std::cout << simulator(L[0], sc21::I_PROB) << "\n";

    { // 線形に高くなる温度を構築
        const double START_TEMP = 0.01, END_TEMP = 15.0;
        for (int i = 0; i < THREAD_NUM; i++) {
            T[i] = START_TEMP + (END_TEMP - START_TEMP) * (double) i / (double) THREAD_NUM;
        }
    }
}

int main() {
    sc21::SC_input();

    const double K = 3.0;
    const double EPS = 1e-7;
    const double TIME_LIMIT = 295;
    const int THREAD_NUM = 48;

    const auto start_t = std::chrono::system_clock::now();
    bool L[THREAD_NUM][N_GROUP][N_GROUP];
    bool BEST[THREAD_NUM][N_GROUP][N_GROUP] = {};
    double temperature[THREAD_NUM], score[THREAD_NUM];
    double best_score[THREAD_NUM];
    initialize(L, temperature, THREAD_NUM);
    std::fill(best_score, best_score + THREAD_NUM, 1e9);

    bool pre_abs = true;

    while (true) {
        // for (int i = 0; i < THREAD_NUM; i++) {
        //     if (i != THREAD_NUM - 1)
        //         printf("%lf ", temperature[i]);
        //     else
        //         printf("%lf\n", temperature[i]);
        // }

        double e;
        {
            auto now_t = std::chrono::system_clock::now();
            e = std::chrono::duration_cast<std::chrono::milliseconds>(now_t - start_t).count() / 1000.0;
        }
        if (e > TIME_LIMIT) break;

        {
            double current_min = best_score[0];
            for (int i = 0; i < THREAD_NUM; i++) {
                current_min = std::min(current_min, best_score[i]);
            }
            for (int i = 0; i < THREAD_NUM; i++) {
                best_score[i] = std::min(best_score[i], current_min + 0.01);
            }
        }

        double temp_scale = (1 - std::min(e / TIME_LIMIT * 1.03, 0.99999));
        // printf("%lf\n", temp_scale);

        bool abs = (e / TIME_LIMIT) < 0.9;

        if (pre_abs != abs) {
            std::fill(best_score, best_score + THREAD_NUM, 1e9);
        }
        pre_abs = abs;

        // 焼きなましを並列でやる
#pragma omp parallel for
        for (int i = 0; i < THREAD_NUM; i++) {
            score[i] = yakinamashi(L[i], temperature[i] * temp_scale, BEST[i], best_score[i], abs);
        }

        {
            double min = best_score[0];
            for (int i = 0; i < THREAD_NUM; i++) {
                min = std::min(min, best_score[i]);
            }
            printf("score: %lf\n", min);
        }

        for (int epoch = 0; epoch < 5; epoch++) {
            for (int thread = 0; thread < THREAD_NUM - 1; thread++) {
                double p;
                {
                    double score_diff = (double) score[thread] - score[thread + 1];
                    double hoge = 1.0 / (K * temperature[thread] + EPS) - 1.0 / (K * temperature[thread + 1] + EPS);
                    p = std::min(1.0, std::exp(score_diff * hoge));
                }
                double q = (double) xor128() / INT32_MAX;
                // printf("exp: %lf\n", temperature[thread]);
                // printf("p: %le\n", p);
                if (q > p) {
                    std::swap(temperature[thread], temperature[thread + 1]);
                }
            }
        }
    }

    int ind = -1;
    {
        double min = 1e9;
        for (int i = 0; i < THREAD_NUM; i++) {
            if (min > best_score[i]) {
                ind = i;
                min = best_score[i];
            }
        }
    }

    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            sc21::C[i][j] = BEST[ind][i][j];
        }
    }

    //提出時にはここを消せ！！！！！！！！！！！
    std::cout << simulator_final(sc21::C, sc21::I_PROB) << "\n";

    sc21::SC_output();
}
