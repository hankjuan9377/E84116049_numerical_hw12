#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const int Nx = 11; // x: 0 to pi ¡÷ 0.1pi*10
const int Ny = 6;  // y: 0 to pi/2 ¡÷ 0.1pi*5
const double h = 0.1 * M_PI;
const double tolerance = 1e-6;
const int max_iter = 10000;

double f(double x, double y) {
    return x * y;
}

int main() {
    double u[Nx][Ny] = {}; // 2D grid for u(x,y)
    double x[Nx], y[Ny];

    // Generate grid
    for (int i = 0; i < Nx; ++i) x[i] = i * h;
    for (int j = 0; j < Ny; ++j) y[j] = j * h;

    // Apply boundary conditions
    for (int j = 0; j < Ny; ++j) {
        u[0][j] = cos(y[j]);          // u(0, y) = cos(y)
        u[Nx - 1][j] = -cos(y[j]);    // u(pi, y) = -cos(y)
    }
    for (int i = 0; i < Nx; ++i) {
        u[i][0] = cos(x[i]);          // u(x, 0) = cos(x)
        u[i][Ny - 1] = 0.0;           // u(x, pi/2) = 0
    }

    // Gauss-Seidel iteration
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_error = 0.0;
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double u_old = u[i][j];
                u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - h * h * f(x[i], y[j]));
                double error = fabs(u[i][j] - u_old);
                if (error > max_error) max_error = error;
            }
        }
        if (max_error < tolerance) {
            cout << "Converged in " << iter + 1 << " iterations.\n";
            break;
        }
    }

    // Output result table
    cout << fixed << setprecision(6);
    cout << "u(x, y) values:\n";
    for (int j = Ny - 1; j >= 0; --j) {
        for (int i = 0; i < Nx; ++i) {
            cout << setw(10) << u[i][j] << " ";
        }
        cout << "\n";
    }

    return 0;
}

