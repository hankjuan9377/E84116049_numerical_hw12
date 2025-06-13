#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

// �Ѽ�
const double dr = 0.1;
const double dt = 0.5;
const double K = 0.1;
const int N_r = 6;         // �d�� [0.5, 1.0]�A�@ 6 �I (�]�t���I)
const int N_t = 21;        // t ? [0,10]�A�Gt=0.5 �� 21 �B
const double r_min = 0.5;
const double r_max = 1.0;

//�p�� alpha

const double alpha = 20;

int main() {
    // �إ� r �y��
    vector<double> r(N_r);
    for (int i = 0; i < N_r; ++i)
        r[i] = r_min + i * dr;

    // �إ� T �x�} T[�ɶ�][�Ŷ�]
    vector<vector<double> > T(N_t);
	for (int i = 0; i < N_t; ++i) {
   	 T[i] = vector<double>(N_r, 0.0);
	}


    // ��l���� T(r,0) = 200(r - 0.5)
    for (int i = 0; i < N_r; ++i) {
        T[0][i] = 200 * (r[i] - 0.5);
    }

    // �ɶ��j��
    for (int n = 0; n < N_t - 1; ++n) {
        double t = n * dt;
		
        // ��s�����`�I i = 1 �� N_r - 2
        for (int i = 1; i < N_r - 1; ++i) {
            double r_i = r[i];
            double term1 = T[n][i + 1] - 2 * T[n][i] + T[n][i - 1];
            double term2 = (T[n][i + 1] - T[n][i]) * (dr / r_i);
            T[n + 1][i] = T[n][i] + alpha * (term1 + term2);
        }

        // ��ɱ���G�k�� T(1, t) = 100 + 40t
        T[n + 1][N_r - 1] = 100 + 40 * (t + dt);

        // ��ɱ���G���� ?T/?r + 3T = 0 �� �ϥΤ@���t������G
        // (T1 - T0)/dr + 3T0 = 0 �� T0 = T1 / (1 + 3dr)
        T[n + 1][0] = T[n][0] * (1 - 3 * dr);
    }

    // ��X���G�]�u��X�̫�@�Ӯɶ��B�^
    cout << fixed << setprecision(4);
    cout << "T(r, t=10):" << endl;
    for (int i = 0; i < N_r; ++i) {
        cout << "r = " << r[i] << ", T = " << T[N_t - 1][i] << endl;
    }

    return 0;
}

