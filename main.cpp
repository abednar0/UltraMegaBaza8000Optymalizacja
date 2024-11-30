#include "opt_alg.h"
#include<iomanip>

void lab0();

void lab1();

void lab2();

void lab3();

void lab4();

void lab5();

void lab6();

int main() {
    try {
        lab3();
    } catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    system("pause");
    return 0;
}

void lab0() {
    //Funkcja testowa
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Wahadlo
    Nmax = 1000;
    epsilon = 1e-2;
    lb = 0;
    ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix();
    Y[1].~matrix();
}

void lab1() {
    try {
        // OPTYMALIZACJE

        string FILE_PATH = R"(C:\Users\Ala\Desktop\optyma\Optimalization.csv)";

        fstream csvFile;
        const int numOfElements = 100;
        double minVal = -15, maxVal = 15;
        set<double> uniqueNums;

        default_random_engine randomEngine(random_device{}());
        uniform_real_distribution<double> realDistribution(minVal, maxVal);

        while (uniqueNums.size() < numOfElements) {
            uniqueNums.insert(realDistribution(randomEngine));
        }

        for (double x0: uniqueNums) {
            // Expansion
            double d = 4.25, alpha = 1.00043;
            int Nmax = 1000;

            unique_ptr<double[]> p(expansion(ff1T, x0, d, alpha, Nmax));

            if (!p) {
                cout << "Expansion returned null for x0: " << x0 << endl;
                continue;
            }

            try {
                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error writing expansion results: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Fibonacci
            try {
                solution Fibonacci = fib(ff1T, p[0], p[1], 0.0001);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Fibonacci.x) << "," << m2d(Fibonacci.y) << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error in Fibonacci: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Lagrange
            try {
                solution Lagrange = lag(ff1T, p[0], p[1], 0.0001, 1e-09, Nmax);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Lagrange.x) << "," << m2d(Lagrange.y) << "," << solution::f_calls << "\n";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error in Lagrange: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();
        }

        // SYMULACJA
        double DA0 = 0.001, d = 0.001, alpha = 2.257, epsilon = 1e-10, gamma = 1e-10;
        int Nmax = 10000;

        unique_ptr<double[]> range(expansion(ff1R, DA0, d, alpha, Nmax));

        if (!range) {
            throw runtime_error("Expansion for simulation returned null");
        }

        solution::clear_calls();

        solution FibonacciSolution = fib(ff1R, range[0], range[1], epsilon);
        cout << "Optimal DA value for Fibonacci's method: " << FibonacciSolution.x << "\n";
        cout << "Value of result function: " << FibonacciSolution.y << "\n";
        cout << "Number of calls of result function: " << solution::f_calls << "\n\n";
        solution::clear_calls();

        solution LagrangeSolution = lag(ff1R, range[0], range[1], epsilon, gamma, Nmax);
        cout << "Optimal DA value for Lagrange's method: " << LagrangeSolution.x << "\n";
        cout << "Value of result function: " << LagrangeSolution.y << "\n";
        cout << "Number of calls of result function: " << solution::f_calls << "\n";
        solution::clear_calls();

        double Pa = 0.5, Va0 = 5, Vb0 = 1, Tb0 = 20, t_0 = 0, t_step = 1, t_end = 2000;
        matrix Y0 = matrix(3, 1, Va0);
        Y0(1) = Vb0;
        Y0(2) = Tb0;

        double Da = m2d(FibonacciSolution.x);
        unique_ptr<matrix[]> FibonacciSimulation(solve_ode(df1R, t_0, t_step, t_end, Y0, Da, Pa));

        FILE_PATH =
                R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\FibonacciSimulation.csv)";
        csvFile.open(FILE_PATH, ios::app);
        if (!csvFile.is_open()) {
            throw runtime_error("Could not open file for Fibonacci simulation");
        }
        csvFile << FibonacciSimulation[1] << "\n";
        csvFile.close();

        solution::clear_calls();

        Da = m2d(LagrangeSolution.x);
        unique_ptr<matrix[]> LagrangeSimulation(solve_ode(df1R, t_0, t_step, t_end, Y0, Da, Pa));

        FILE_PATH =
                R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\LagrangeSimulation.csv)";
        csvFile.open(FILE_PATH, ios::app);
        if (!csvFile.is_open()) {
            throw runtime_error("Could not open file for Lagrange simulation");
        }
        csvFile << LagrangeSimulation[1] << "\n";
        csvFile.close();

        solution::clear_calls();
    } catch (const exception &e) {
        cout << "Fatal error in lab1: " << e.what() << endl;
    } catch (...) {
        cout << "Unknown error in lab1" << endl;
    }
}


void lab2() {
    try {
        std::string FILE_PATH = R"(C:\Users\Ala\Desktop\optyma\Optimalization.csv)";
        std::ofstream csvFile(FILE_PATH, std::ios::out);
        if (!csvFile.is_open()) throw std::runtime_error("Could not open file");

        double minVal = -1.0, maxVal = 1.0;
        std::default_random_engine randomEngine(std::random_device{}());
        std::uniform_real_distribution<double> realDistribution(minVal, maxVal);

        double epsilon = 1e-6;
        int Nmax = 1000;
        matrix ud1, ud2;
        std::vector<double> step_sizes = {0.44, 0.22, 0.11};


            matrix x0(2, 1);
            x0(0, 0) = realDistribution(randomEngine);
            x0(1, 0) = realDistribution(randomEngine);


                solution Xopt_HJ = HJ(ff2T, x0, 0.22, 0.45, epsilon, Nmax, ud1, ud2);
                solution::clear_calls();

                matrix s0(2, 1, 0.22);
                solution Xopt_Rosen = Rosen(ff2T, x0, s0, 1.22, 0.55, epsilon, Nmax, ud1, ud2);
                solution::clear_calls();



        csvFile.close();
    } catch (const std::exception &e) {
        std::cerr << "Fatal error in lab2: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error in lab2" << std::endl;
    }
}

matrix generate_random_point() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist_x1(-M_PI, M_PI); // Zakres dla x1
    std::uniform_real_distribution<> dist_x2(-M_PI, M_PI); // Zakres dla x2

    double vals[2] = {dist_x1(gen), dist_x2(gen)};
    return matrix(2, vals); // Macierz 2x1
}

void lab3() {
    try {
        // Ustawienia dla testu
        matrix x0 = matrix(2, 1, 3.0);  // Punkt pocz¹tkowy
        double penalty_start = 1.0;     // Pocz¹tkowy wspó³czynnik kary
        double penalty_scale_ext = 2.0; // Skalowanie dla zewnêtrznej funkcji kary
        double penalty_scale_int = 0.5; // Skalowanie dla wewnêtrznej funkcji kary
        double epsilon = 1e-3;          // Kryterium dok³adnoœci
        int Nmax = 10000;               // Maksymalna liczba iteracji

        // Parametr ograniczenia (wartoœæ a)
        matrix a_value = matrix(1, 1, 4.4934);

        // Optymalizacja zewnêtrzna
        solution external_opt = pen(fff3T, x0, penalty_start, penalty_scale_ext, epsilon, Nmax, a_value, matrix(2, new double[2]{1.0, 2.0}));
        std::cout << "Zewnêtrzna kara:" << std::endl;
        std::cout << "x_opt: " << external_opt.x(0) << ", " << external_opt.x(1) << std::endl;
        std::cout << "y_opt: " << external_opt.y(0) << std::endl;
        std::cout << "Liczba wywo³añ funkcji celu: " << solution::f_calls << std::endl;

        solution::clear_calls();

        // Optymalizacja wewnêtrzna
        solution internal_opt = pen(fff3T, x0, penalty_start, penalty_scale_int, epsilon, Nmax, a_value, matrix(2, new double[2]{1.0, 0.5}));
        std::cout << "Wewnêtrzna kara:" << std::endl;
        std::cout << "x_opt: " << internal_opt.x(0) << ", " << internal_opt.x(1) << std::endl;
        std::cout << "y_opt: " << internal_opt.y(0) << std::endl;
        std::cout << "Liczba wywo³añ funkcji celu: " << solution::f_calls << std::endl;

    }
    catch (const std::exception &e) {
        std::cerr << "Fatal error in lab3: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error in lab3" << std::endl;
    }
}

void lab4() {
}

void lab5() {
}

void lab6() {
}
