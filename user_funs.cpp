#include"user_funs.h"
#include <functional>

matrix ff0T(matrix x, matrix ud1, matrix ud2) {
    matrix y(1, 1);
    y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
    return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
    int n = get_len(Y[0]);
    double teta_max = Y[1](0, 0);
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1));
    Y[0].~matrix();
    Y[1].~matrix();
    return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(2, 1);
    double m = 1, l = 0.5, b = 0.5, g = 9.81;
    double I = m * pow(l, 2);
    dY(0) = Y(1);
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
    return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    return {-cos(0.1 * m2d(x)) * exp(-pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow(0.1 * m2d(x), 2)};
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    matrix Y0(3, 1);
    Y0(0) = 5.0;
    Y0(1) = 1.0;
    Y0(2) = 20.0;
    double t0 = 0.0;
    double dt = 1.0;
    int timesteps = 2000;

    matrix *Y = solve_ode(df1R, t0, dt, timesteps, Y0, x(0), x);

    double Tmax = Y[1](0, 2);
    for (int i = 1; i <= timesteps; ++i) {
        if (Y[1](i, 2) > Tmax) {
            Tmax = Y[1](i, 2);
        }
    }
    delete[] Y;

    return {fabs(Tmax - 50.0)};
}

matrix df1R(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix derivatives(3, 1);
    double a = 0.98, b = 0.63, g = 9.81, PA = 0.5, PB = 1, DB = 0.00365665;
    double F_in = 0.01, T_in = 20, TA = 90.0;
    double DA = m2d(ud1(0));

    double FA_out, FB_out;

    if (Y(0) > 0)
        FA_out = a * b * DA * sqrt(2 * g * Y(0) / PA);
    else
        FA_out = 0;

    if (Y(1) > 0)
        FB_out = a * b * DB * sqrt(2 * g * Y(1) / PB);
    else
        FB_out = 0;

    derivatives(0) = -FA_out;
    derivatives(1) = FA_out + F_in - FB_out;
    derivatives(2) = FA_out / Y(1) * (TA - Y(2)) + F_in / Y(1) * (T_in - Y(2));

    return derivatives;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    matrix y(1, 1);
    y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
    return y;
}

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2) {
    double alpha = Y(0);
    double omega = Y(1);

    double k1 = ud2(0);
    double k2 = ud2(1);

    double m_arm = 1.0;
    double m_weight = 5.0;
    double l = 1.0;
    double g = 9.81;
    double b = 0.5;

    double I = (m_arm + m_weight) * pow(l, 2);

    double M = k1 * (M_PI - alpha) + k2 * (0.0 - omega);

    matrix dY(2, 1);
    dY(0) = omega;
    dY(1) = (M - m_arm * g * l * sin(alpha) - b * omega) / I;

    return dY;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    double k1 = x(0);
    double k2 = x(1);

    double m = 1.0, l = 1.0, g = 9.81, b = 0.5;
    double alpha_target = M_PI;
    double omega_target = 0.0;

    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix control_params(2, 1);
    control_params(0) = k1;
    control_params(1) = k2;

    matrix *result = solve_ode(df2R, 0.0, 0.1, 100.0, Y0, matrix(), control_params);

    double Q = 0.0;
    for (int i = 0; i < get_len(result[0]); ++i) {
        double alpha = result[1](i, 0);
        double omega = result[1](i, 1);
        double M = k1 * (alpha_target - alpha) + k2 * (omega_target - omega);
        Q += 10 * pow(alpha_target - alpha, 2) + pow(omega_target - omega, 2) + pow(M, 2);
    }

    delete[] result;

    return {Q};
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    // Wyci�ganie wsp�rz�dnych x1 i x2 z macierzy x
    double x1 = x(0, 0);
    double x2 = x(1, 0);

    // Obliczanie mianownika: sqrt((x1 / pi)^2 + (x2 / pi)^2)
    double denominator = sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2));


    // Obliczanie warto�ci funkcji celu
    double result = sin(M_PI * denominator) / (M_PI * denominator);

    // Zwracanie wyniku jako macierz 1x1
    return matrix(result);
}

bool g1(matrix x1) {
    return (-x1() + 1 <= 0);
}

bool g2(matrix x2) {
    return (-x2() + 1 <= 0);
}

bool g3(matrix x1, matrix x2, double alpha) {
    return (sqrt(pow(x1(), 2) + pow(x2(), 2)) - alpha <= 0);
}

// matrix fR3(matrix x, matrix ud1, matrix ud2) {
//     matrix y;
//     matrix Y0(4, new double[4]{ 0, x(0), 100, 0 });
//     matrix* Y = solve_ode(fd3T, 0, 0.01, 7, Y0, ud1, x(1));
//     int n = get_len(Y[0]);
//     int i0 = 0, i50 = 0;
//
//     for (int i = 0; i < n; i++) {
//         if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
//             i50 = i;
//         if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
//             i0 = i;
//     }
//
//     y = -Y[1](i0, 0);
//
//     if (abs(x(0)) - 10 > 0)
//         y = y + ud2 * pow(abs(x(0)) - 10, 2);
//     if (abs(x(1)) - 25 > 0)
//         y = y + ud2 * pow(abs(x(1)) - 25, 2);
//     if (abs(Y[1](i50, 0) - 5) - 1 > 0)
//         y = y + ud2 * pow(abs(Y[1](i50, 0) - 5) - 1, 2);
//
//     return y;
// }

matrix fff3T(matrix x, matrix ud1, matrix ud2) {
    double argument = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
    matrix y = sin(argument) / argument;

    // Definicja ogranicze�
    std::vector<std::function<double(matrix)>> constraints = {
        [](matrix x) { return -x(0) + 1; },           // g1: -x1 + 1 <= 0
        [](matrix x) { return -x(1) + 1; },           // g2: -x2 + 1 <= 0
        [&ud1](matrix x) { return norm(x) - ud1(0); } // g3: ||x|| - a <= 0
    };

    // Zewn�trzna kara
    if (ud2(1) > 1) {
        for (const auto &g : constraints) {
            double violation = std::max(0.0, g(x));
            y = y + ud2(0) * pow(violation, 2);
        }
    }
    // Wewn�trzna kara
    else {
        for (const auto &g : constraints) {
            double violation = g(x);
            if (violation >= 0) {
                return matrix(1, 1, 1e10); // Punkt poza obszarem dopuszczalnym
            }
            y = y - ud2(0) / (-violation);
        }
    }

    return y;
}