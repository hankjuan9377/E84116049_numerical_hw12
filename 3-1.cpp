#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const int Nr = 6;       // r: 0.5 to 1, step = 0.1
const int Nth = 6;      // theta: 0 to pi/3, step = pi/15
const double r0 = 0.5;
const double dr = 0.1;
const double dtheta = M_PI / 15;
const double tol = 1e-5;
const int max_iter = 10000;

int main() {
    double T[Nr][Nth] = {};  // T[i][j] = T(r_i, theta_j)
    double r[Nr], theta[Nth];

    // 初始化 r, theta
    for (int i = 0; i < Nr; ++i)
        r[i] = r0 + i * dr;
    for (int j = 0; j < Nth; ++j)
        theta[j] = j * dtheta;

    // 邊界條件
    for (int j = 0; j < Nth; ++j) {
        T[0][j] = 50;   // r = 0.5
        T[Nr-1][j] = 100; // r = 1
    }
    for (int i = 0; i < Nr; ++i) {
        T[i][0] = 0;              // theta = 0
        T[i][Nth-1] = 0;          // theta = pi/3
    }

    // Gauss-Seidel 迭代
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_err = 0.0;

        for (int i = 1; i < Nr - 1; ++i) {
            for (int j = 1; j < Nth - 1; ++j) {
                double ri = r[i];
                double term_r = (T[i+1][j] + T[i-1][j]) / (dr * dr);
                double term_rr = (T[i+1][j] - T[i-1][j]) / (2 * dr);
                double term_th = (T[i][j+1] + T[i][j-1]) / (dtheta * dtheta);

                double T_new = (
                    term_r +
                    (1.0 / ri) * term_rr +
                    (1.0 / (ri * ri)) * term_th
                ) / (2.0 / (dr * dr) + 2.0 / (ri * ri * dtheta * dtheta));

                double err = fabs(T_new - T[i][j]);
                T[i][j] = T_new;
                if (err > max_err) max_err = err;
            }
        }

        if (max_err < tol) {
            cout << "Converged in " << iter + 1 << " iterations.\n";
            break;
        }
    }

    // 輸出結果
    cout << fixed << setprecision(3);
    cout << "\nT(r, theta) 2D Grid (r: rows, theta: columns):\n";
    for (int i = 0; i < Nr; ++i) {
        for (int j = 0; j < Nth; ++j) {
            cout << setw(7) << T[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}

