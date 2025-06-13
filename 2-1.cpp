#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

// 參數
const double dr = 0.1;
const double dt = 0.5;
const double K = 0.1;
const int N_r = 6;         // 範圍 [0.5, 1.0]，共 6 點 (包含端點)
const int N_t = 21;        // t ? [0,10]，Δt=0.5 → 21 步
const double r_min = 0.5;
const double r_max = 1.0;

//計算 alpha

const double alpha = 20;

int main() {
    // 建立 r 座標
    vector<double> r(N_r);
    for (int i = 0; i < N_r; ++i)
        r[i] = r_min + i * dr;

    // 建立 T 矩陣 T[時間][空間]
    vector<vector<double> > T(N_t);
	for (int i = 0; i < N_t; ++i) {
   	 T[i] = vector<double>(N_r, 0.0);
	}


    // 初始條件 T(r,0) = 200(r - 0.5)
    for (int i = 0; i < N_r; ++i) {
        T[0][i] = 200 * (r[i] - 0.5);
    }

    // 時間迴圈
    for (int n = 0; n < N_t - 1; ++n) {
        double t = n * dt;
		
        // 更新內部節點 i = 1 到 N_r - 2
        for (int i = 1; i < N_r - 1; ++i) {
            double r_i = r[i];
            double term1 = T[n][i + 1] - 2 * T[n][i] + T[n][i - 1];
            double term2 = (T[n][i + 1] - T[n][i]) * (dr / r_i);
            T[n + 1][i] = T[n][i] + alpha * (term1 + term2);
        }

        // 邊界條件：右邊 T(1, t) = 100 + 40t
        T[n + 1][N_r - 1] = 100 + 40 * (t + dt);

        // 邊界條件：左邊 ?T/?r + 3T = 0 → 使用一階差分近似：
        // (T1 - T0)/dr + 3T0 = 0 → T0 = T1 / (1 + 3dr)
        T[n + 1][0] = T[n][0] * (1 - 3 * dr);
    }

    // 輸出結果（只輸出最後一個時間步）
    cout << fixed << setprecision(4);
    cout << "T(r, t=10):" << endl;
    for (int i = 0; i < N_r; ++i) {
        cout << "r = " << r[i] << ", T = " << T[N_t - 1][i] << endl;
    }

    return 0;
}

