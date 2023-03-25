#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

int factorial(int n) {
    int f = 1;
    for (int i = 2; i <= n; i++) {
        f *= i;
    }
    return f;
}


double Wigner_d(int j, int m1, int m2, double theta) {
    double d = 0;
    double s1 = 1.0, s2 = 1.0, s3 = 1.0, s4 = 1.0;
    double a = cos(theta / 2);
    double b = sin(theta / 2);

    for (int k = 0; k <= j; k++) {
        if ((j - m1 - k) < 0 || (j + m2 - k) < 0 || (m1 - m2 + k) < 0) {
            continue;
        }

        s1 = pow(-1.0, k);
        s2 = sqrt(factorial(j + m1) * factorial(j - m1) * factorial(j + m2) * factorial(j - m2));
        s3 = pow(b, j - m1 - k);
        s4 = pow(a, j + m2 - k);

        d += s1 / (factorial(j - m1 - k) * factorial(m1 - m2 + k) * factorial(j + m2 - k)) * s2 * s3 * s4;
    }

    return d;
}

vector<vector<double>> rotation_matrix(int l, double theta) {
    vector<vector<double>> R((l + 1) * (l + 1), vector<double>((l + 1) * (l + 1), 0.0));

    for (int i = 0; i <= l; i++) {
        for (int m1 = -i; m1 <= i; m1++) {
            int index1 = i * i + m1 + i;
            for (int j = 0; j <= l; j++) {
                for (int m2 = -j; m2 <= j; m2++) {
                    int index2 = j * j + m2 + j;
                    double d = Wigner_d(j, m2, m1, theta);
                    R[index1][index2] = d;
                }
            }
        }
    }

    return R;
}


int main() {
    int l = 2;
    double theta = M_PI / 4.0;
    vector<vector<double>> R = rotation_matrix(l, theta);

    cout << "Rotation matrix R for l = " << l << " and theta = " << theta << ":" << endl;
    for (int i = 0; i < (l + 1) * (l + 1); i++) {
        for (int j = 0; j < (l + 1) * (l + 1); j++) {
            cout << R[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}

