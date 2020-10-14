#include <iostream>
#include <iomanip>

double k(double r) {
    //return r * r;
    return r;
}

double u(double r) {
    //return 15;
    return r * r;
}

double func(double r) {
    //return 15;
    return r * r - 6 * r;
}

double r(double i, double h) {
    return i * h;
}

double q(double r) {
    return 1;
}

int main() {
    double R = 5; // [0, 5] - r
    double X = 1;
    int N = 10;
    double h = R / N;
    double h2 = h / 2;
    double nu = 75; // 15
    // ri = i * h; fi = func(ri)
    // Av = F, размерность - N + 1

    for (int j = 0; j < 19; j++) {
        h = R / N;
        h2 = h / 2;
        auto F = new double[N + 1];
        F[0] = 1.0 / 2 * h * h2 / 2 * func(0);
        for (int i = 1; i < N + 1; ++i) {
            F[i] = h * i * h * func(i * h);
        }
        F[N] = h2 * N * h * func(N * h) + N * h * nu;

        // Зададимся тремя векторами, определяющими трёхдиагональную матрицу А
        auto A = new double[N + 1];
        A[0] = 0;
        for (int i = 1; i < N + 1; ++i) {
            A[i] = -(i - 1.0 / 2) * h * k((i - 1.0 / 2) * h) / h;
        }

        auto B = new double[N + 1];
        B[N] = 0;
        for (int i = 0; i < N + 1; ++i) {
            B[i] = -(i + 1.0 / 2) * h * k((i + 1.0 / 2) * h) / h;
        }

        auto C = new double[N + 1];
        C[0] = r(1.0 / 2, h) * k(r(1.0 / 2, h)) / h + r(1.0 / 2, h) / 2 * h2 * q(r(0, h));
        for (int i = 1; i < N; ++i) {
            C[i] = r(i + 1.0 / 2, h) * k(r(i + 1.0 / 2, h)) / h + r(i - 1.0 / 2, h) * k(r(i - 1.0 / 2, h)) / h +
                   h * r(i, h) * q(r(i, h));
        }
        C[N] = r(N - 1.0 / 2, h) * k(r(N - 1.0 / 2, h)) / h + h2 * r(N, h) * q(r(N, h)) + r(N, h) * X;

        // прогонка
        auto alpha = new double[N + 1];
        auto beta = new double[N + 1];

        alpha[1] = -B[0] / C[0];
        beta[1] = F[0] / C[0];
        for (int i = 1; i < N; ++i) {
            alpha[i + 1] = -B[i] / (A[i] * alpha[i] + C[i]);
            beta[i + 1] = (F[i] - A[i] * beta[i]) / (A[i] * alpha[i] + C[i]);
        }

        auto v = new double[N + 1];
        v[N] = (F[N] - A[N] * beta[N]) / (A[N] * alpha[N] + C[N]);

        for (int i = N - 1; i > -1; --i) {
            v[i] = alpha[i + 1] * v[i + 1] + beta[i + 1];
        }

        //std::cout << "  r\t     v[i]" << "\t\t" << "u(r)" << "\n";
        double err = 0;
        for (int i = 0; i < N + 1; ++i) {
            if (err < std::abs(v[i] - u(r(i, h)))) {
                err = std::abs(v[i] - u(r(i, h)));
            }
            //std::cout << std::setw(4) << r(i, h) << std::setw(16) << v[i] << std::setw(16) << u(r(i, h)) << std::setw(16) << v[i] - u(r(i, h)) << "\n";
        }
        std::cout << std::setw(10) << N << "\t";
        std::cout << std::setw(15) << err << "\n";
        N *= 2;

        delete[] A;
        delete[] B;
        delete[] C;
        delete[] F;
        delete[] v;
        delete[] alpha;
        delete[] beta;
    }

    return 0;
}
