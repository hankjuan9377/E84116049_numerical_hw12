#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

const double dr = 0.1;
const double dt = 0.5;
const double K = 0.1;
const int N_r = 6;         // 範圍 [0.5, 1.0]
const int N_t = 21;
const double r_min = 0.5;

const double alpha = 4 * K * dt / (dr * dr);  // α = 4KΔt / Δr2

int main() {
    vector<double> r(N_r);
    for (int i = 0; i < N_r; ++i)
        r[i] = r_min + i * dr;

    vector<vector<double> > T(N_t);
    for (int i = 0; i < N_t; ++i)
        T[i] = vector<double>(N_r, 0.0);

    // 初始條件
    for (int i = 0; i < N_r; ++i)
        T[0][i] = 200 * (r[i] - 0.5);

    // 時間迴圈
    for (int n = 0; n < N_t - 1; ++n) {
        double t = n * dt;

        // 三對角係數
        vector<double> a(N_r - 2), b(N_r - 2), c(N_r - 2), d(N_r - 2);

        for (int i = 1; i < N_r - 1; ++i) {
            double ri = r[i];
            double coeff1 = alpha;
            double coeff2 = alpha * dr / (2 * ri);

            a[i - 1] = -coeff1 + coeff2;
            b[i - 1] = 1 + 2 * coeff1;
            c[i - 1] = -coeff1 - coeff2;
            d[i - 1] = T[n][i];
        }

        // 修正 d 對邊界的影響
        // 左邊界: ?T/?r + 3T = 0 ? T? = T? / (1 + 3Δr)
        d[0] -= a[0] * (T[n + 1][0] = T[n + 1][1] / (1 + 3 * dr));
        // 右邊界: T = 100 + 40t
        T[n + 1][N_r - 1] = 100 + 40 * (t + dt);
        d[N_r - 3] -= c[N_r - 3] * T[n + 1][N_r - 1];

        // Thomas algorithm
        for (int i = 1; i < N_r - 2; ++i) {
            double m = a[i] / b[i - 1];
            b[i] -= m * c[i - 1];
            d[i] -= m * d[i - 1];
        }

        vector<double> x(N_r - 2);
        x[N_r - 3] = d[N_r - 3] / b[N_r - 3];
        for (int i = N_r - 4; i >= 0; --i)
            x[i] = (d[i] - c[i] * x[i + 1]) / b[i];

        // 更新 T[n+1]
        for (int i = 1; i < N_r - 1; ++i)
            T[n + 1][i] = x[i - 1];
    }

    // 輸出最後時間的結果
    cout << fixed << setprecision(4);
    cout << "T(r, t=10):" << endl;
    for (int j = 0; j < N_t; ++j) {
   		for (int i = 0; i < N_r; ++i) {
     	   cout << "r = " << r[i] << "t = " << 0+0.5*j << ", T = " << T[j][i] << endl;
    	}
	}

    return 0;
}

