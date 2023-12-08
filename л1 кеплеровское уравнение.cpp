#include <iostream>
#include <cmath>
#include <vector>

//double calculateFunction(double M, double e, double E);
double M(int T, int i, double n);
double iteration_method(double M0, double epsilon, double E, double e, int T, int i, double n);
//void method_of_half_devision(double M, double e, double epsilon, double pi);
//void golden_ratio_method(double M, double e, double epsilon, double pi);
//void newton_method(double M, double e, double epsilon);

int main()
{
    double epsilon = pow(10, -8);
    double pi{ 3.1415926535 };
    double Mercury_weight = 0.33 * pow(10, 24);
    double gravity_const = 6.67430 * pow(10, -20);
    double a = 7689.7;
    double e = 0.656722629;
    double M0 = 0;
    double E0 = 0;
    int i = 0;
    double mu = gravity_const * Mercury_weight;
    int T = 2*pi*sqrt(pow(a, 3)/mu);
    double n = sqrt(mu / (a * a * a));

    iteration_method(M0, epsilon, E0, e, T, i, n);
    /*method_of_half_devision(M, e, epsilon, pi);
    golden_ratio_method(M, e, epsilon, pi);
    newton_method(M, e, epsilon);*/
}

double M(int T, int i, double n) {
    std::vector<double> M(T);
    for (int t = 0; t < T; ++t) M[t] = t * n;
    return M[i];
}

// Функция для вычисления значения функции F(E) = M + e*sin(E) - E
//double calculateFunction(double M, double e, double E) {
//    return M + e * sin(E) - E;
//}

double iteration_method(double M0, double epsilon, double E0, double e, int T, int i, double n) {
    double E = M0 + e * sin(E0);
    std::cout << E << "\n";
    if (abs(E - E0) < epsilon) {
        return E;
    }
    else {
        M0 = M(T, i, n);
        i++;
        return iteration_method(M0, epsilon, E, e, T, i, n);
    }
}

//void method_of_half_devision(double M, double e, double epsilon, double pi) {
//    double a = -3.0 * pi;
//    double b = 2.0 * pi;
//
//    while (fabs(b - a) >= epsilon) {
//        double c = (a + b) / 2.0;
//        double functionvalue = c - e * sin(c) - M;
//    
//        if (functionvalue == 0.0)
//            break;
//        else if (functionvalue * (a - b) < 0)
//            b = c;
//        else
//            a = c;
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
