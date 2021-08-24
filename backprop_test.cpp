#include "./main.cpp"

// CとI_PROBを与えると誤差を返す関数。
double simulator_double(double C[][N_GROUP], double I_PROB[]) {
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

// void diff(double C[][N_GROUP], double I_PROB[], double C_back[][N_GROUP]) {
//     std::fill(C_back[0], C_back[N_GROUP], 0.0);
//     double eps = 1e-3;

//     for (int i = 0; i < 1; i++) {
//         for (int j = 0; j < N_GROUP; j++) {
//             C[i][j] -= eps;
//             double from = simulator(C, I_PROB);
//             C[i][j] += 2 * eps;
//             double to = simulator(C, I_PROB);
//             C[i][j] -= eps;
//             C_back[i][j] = (to - from) / eps / 2;
//         }
//     }
// }
int main() {
    sc21::SC_input();

    double C_back_1[N_GROUP][N_GROUP], C_back_2[N_GROUP][N_GROUP];
    double C_double[N_GROUP][N_GROUP];
    for (int i = 0; i < N_GROUP; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            sc21::C[i][j] = (i + j) % 2;
            C_double[i][j] = (double)sc21::C[i][j];
        }
    }

    backprop(sc21::C, sc21::I_PROB, C_back_1);

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < N_GROUP; j++) {
            printf("(%d, %d) back: %le, diff: %le\n", i, j, C_back_1[i][j],
                   C_back_2[i][j]);
        }
    }
}