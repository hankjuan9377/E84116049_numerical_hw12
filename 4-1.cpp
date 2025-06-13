#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const int Nx = 11;        // x: 0.0 to 1.0, step = 0.1
const int Nt = 21;        // 模擬 t=0 ~ 2.0, step = 0.1
const double dx = 0.1;
const double dt = 0.1;
const double pi = acos(-1.0);
const double lambda2 = (dt * dt) / (dx * dx); // λ^2 = 1

double initial_u(double x) {
    return cos(2 * pi * x);
}

double initial_ut(double x) {
    return 2 * pi * sin(2 * pi * x);
}

int main() {
    double u_prev[Nx] = {}; // p^0
    double u_curr[Nx] = {}; // p^1
    double u_next[Nx] = {}; // p^{n+1}
    double x[Nx];

    for (int i = 0; i < Nx; ++i)
        x[i] = i * dx;

    // 初始狀態 t = 0
    for (int i = 0; i < Nx; ++i)
        u_prev[i] = initial_u(x[i]);

    // 套用邊界 t=0
    u_prev[0] = 1.0;
    u_prev[Nx - 1] = 2.0;

    // 使用 Taylor 展開得 t = dt 的值（u_curr）
    for (int i = 1; i < Nx - 1; ++i) {
        double utt = (u_prev[i+1] - 2*u_prev[i] + u_prev[i-1]) / (dx*dx);
        u_curr[i] = u_prev[i] + dt * initial_ut(x[i]) + 0.5 * dt * dt * utt;
    }
    // 邊界值
    u_curr[0] = 1.0;
    u_curr[Nx - 1] = 2.0;

    // 輸出 t = 0, t = dt
    cout << fixed << setprecision(5);
    cout << "t = 0.0\n";
    for (int i = 0; i < Nx; ++i)
        cout << "x = " << setw(4) << x[i] << ", p = " << u_prev[i] << endl;

    cout << "\nt = 0.1\n";
    for (int i = 0; i < Nx; ++i)
        cout << "x = " << setw(4) << x[i] << ", p = " << u_curr[i] << endl;

    // 往後時間層迭代
    for (int n = 2; n <= 10; ++n) {
        for (int i = 1; i < Nx - 1; ++i) {
            u_next[i] = 2 * u_curr[i] - u_prev[i] +
                        lambda2 * (u_curr[i+1] - 2 * u_curr[i] + u_curr[i-1]);
        }
        u_next[0] = 1.0;
        u_next[Nx - 1] = 2.0;

        // 輸出每個時間層
        cout << "\nt = " << n * dt << endl;
        for (int i = 0; i < Nx; ++i)
            cout << "x = " << setw(4) << x[i] << ", p = " << u_next[i] << endl;

        // 更新時間層
        for (int i = 0; i < Nx; ++i) {
            u_prev[i] = u_curr[i];
            u_curr[i] = u_next[i];
        }
    }

    return 0;
}

