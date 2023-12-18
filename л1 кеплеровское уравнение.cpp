#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

const double epsilon = pow(10, -4);
const double pi{ 3.1415926535 };
const double Mercury_weight = 0.33 * pow(10, 24);
const double gravity_const = 6.67430 * pow(10, -20);
const double a = 7689.7;
const double e = 0.656722629;
const double mu = gravity_const * Mercury_weight;
const double T = 2 * pi * sqrt(pow(a, 3) / mu);
const double n = sqrt(mu / (a * a * a));
//double calculateFunction(double M, double e, double E);
double M(int t);
double iteration_method(double M0, double E, int i);
//void method_of_half_devision(double M0, double e, double epsilon, double pi, std::vector<double> M, int i);
//void golden_ratio_method(double M, double e, double epsilon, double pi);
//void newton_method(double M, double e, double epsilon);

int main()
{
    int i = 0;
    
    //6.28284
    int v = iteration_method(i * n, i * n, i);
    std::cout << "end " << v << " " << T * n + e * sin(v);
    /*method_of_half_devision(M0, e, epsilon, pi, std::vector<double> M, int i);
    golden_ratio_method(M, e, epsilon, pi);
    newton_method(M, e, epsilon);*/
}

// Функция для вычисления значения функции F(E) = M + e*sin(E) - E
//double calculateFunction(double M, double e, double E) {
//    return M + e * sin(E) - E;
//}

double iteration_method(double M0, double E0, int i) {
    while(i < T) {
        double E = M0 + e * sin(E0);
        //std::cout << M0 << "\t" << E << "\n";
        std::ofstream chartE("chartE.txt", std::ios::app);
        if (chartE.is_open()) {
            chartE << i << "\t" << M0 << "\t" << E << "\n";
        }
        chartE.close();

        if (abs(E - E0) < epsilon and E != 0) {
            return E;
        }
        else {
            ++i;
            M0 = i * n;
            E0 = E;
        }
    }
}

//void method_of_half_devision(double M0, double e, double epsilon, double pi, std::vector<double> M, int i) {
//    double a = -3.0 * pi;
//    double b = 2.0 * pi;
//
//    while (fabs(b - a) >= epsilon) {
//        double c = (a + b) / 2.0;
//        double functionvalue = c - e * sin(c) - M0;
//    
//        if (functionvalue == 0.0)
//            break;
//        else if (functionvalue * (a - b) < 0)
//            b = c;
//        else
//            a = c;
//        M0 = M[i];
//        i++;
//    }
//    std::cout << (a + b) / 2.0 << "\n";
//}
//
//void golden_ratio_method(double M, double e, double epsilon, double pi) {
//    const double goldenRatio = (1 + sqrt(5)) / 2; // Золотое сечение
//    
//        double a = 0;
//        double b = 2 * pi; // Пределы для эксцентрической аномалии
//        double c = b - (b - a) / goldenRatio;
//        double d = a + (b - a) / goldenRatio;
//    
//        while (abs(b - a) > 2 * epsilon) {
//            if (calculateFunction(M, e, c) * calculateFunction(M, e, d) < 0) {
//                b = d;
//            }
//            else {
//                a = c;
//            }
//    
//            c = b - (b - a) / goldenRatio;
//            d = a + (b - a) / goldenRatio;
//        }
//    
//        // Найденная эксцентрическая аномалия
//        double E = (b + a) / 2;
//        std::cout << E << "\n";
//}
//
//void newton_method(double M, double e, double epsilon) {
//    // Начальное приближение для эксцентрической аномалии
//    double E = M;
//
//    // Максимальное количество итераций
//    int maxIterations = 1000;
//
//    // Итерации методом Ньютона
//    for (int i = 0; i < maxIterations; ++i) {
//        double f = E - e * sin(E) - M;
//        double f_prime = 1.0 - e * cos(E);
//
//        // Метод Ньютона
//        double delta = f / f_prime;
//        E = E - delta;
//
//        // Проверка условия сходимости
//        if (fabs(delta) < epsilon || fabs(f) < epsilon) {
//            break;
//        }
//    }
//
//    std::cout << E << "\n";
//}
