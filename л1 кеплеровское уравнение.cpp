#include <iostream>
#include <cmath>

double iteration_method(double M, double epsilon, double E, double e);
void method_of_half_devision(double M, double e, double epsilon, double pi);
void golden_ratio_method(double M, double e, double accuracy);
void newton_method();

int main()
{
    double epsilon = pow(10, -8);
    double pi{ 3.1415926535 };
    double Mercury_weight = 0.33 * pow(10, 24);
    double gravity_const = 6.67430 * pow(10, -20);
    double a = 7689.7;
    double e = 0.656722629;
    double T = 12 * 60 * 60;
    double M = 0;
    
    double mu = gravity_const * Mercury_weight;
    double n = sqrt(mu / (a * a * a));
    for (int t = 0; t < T; ++t) {
        M = (t - 0) * n;
    }
    iteration_method(M, epsilon, M, e);
}

double iteration_method(double M, double epsilon, double E0, double e) {
    double E = M + e * sin(E0);
    std::cout << E << "\n";
    if (abs(E - E0) < epsilon) {
        std::cout << E;
        return E;
    }
    else {
        return iteration_method(M, epsilon, E, e);
    }
}

void method_of_half_devision(double M, double e, double epsilon, double pi) {
    double a = -3.0 * pi;
    double b = 2.0 * pi;

    while (fabs(b - a) >= epsilon) {
        double c = (a + b) / 2.0;
        double functionvalue = c - e * sin(c) - M;
    
        if (functionvalue == 0.0)
            break;
        else if (functionvalue * (a - b) < 0)
            b = c;
        else
            a = c;
    }
}

void golden_ratio_method(double M, double e, double accuracy) {
    const double goldenRatio = (1 + sqrt(5)) / 2; // Золотое сечение
    //
    //    double a = 0;
    //    double b = 2 * M_PI; // Пределы для эксцентрической аномалии
    //    double c = b - (b - a) / goldenRatio;
    //    double d = a + (b - a) / goldenRatio;
    //
    //    while (abs(b - a) > 2 * accuracy) {
    //        if (calculateFunction(M, e, c) * calculateFunction(M, e, d) < 0) {
    //            b = d;
    //        }
    //        else {
    //            a = c;
    //        }
    //
    //        c = b - (b - a) / goldenRatio;
    //        d = a + (b - a) / goldenRatio;
    //    }
    //
    //    // Найденная эксцентрическая аномалия
    //    double E = (b + a) / 2;
    //    return E;
}

void newton_method(double meanAnomaly, double eccentricity, double tolerance) {
    // Начальное приближение для эксцентрической аномалии
//    double E = meanAnomaly;
//
//    // Максимальное количество итераций
//    int maxIterations = 1000;
//
//    // Итерации методом Ньютона
//    for (int i = 0; i < maxIterations; ++i) {
//        double f = E - eccentricity * sin(E) - meanAnomaly;
//        double f_prime = 1.0 - eccentricity * cos(E);
//
//        // Метод Ньютона
//        double delta = f / f_prime;
//        E = E - delta;
//
//        // Проверка условия сходимости
//        if (fabs(delta) < tolerance || fabs(f) < tolerance) {
//            break;
//        }
//    }
//
//    return E;
}


//int main() {
//    // Заданные параметры орбиты
//    double eccentricity = 0.957;
//
//    // Итерации для каждого момента времени от 0 до 28080
//    for (int time = 0; time <= 28080; ++time) {
//        // Приведение момента времени к радианам
//        double meanAnomalyRad = (2.0 * PI * time) / 28080.0;
//
//        // Вычисление эксцентрической аномалии с точностью от 10^(-3) до 10^(-10)
//        for (int i = 3; i <= 10; ++i) {
//            double precision = pow(10, -i);
//            double eccentricAnomaly = calculateEccentricAnomaly(meanAnomalyRad, eccentricity, precision);
//            double calculatedMeanAnomaly = meanAnomaly(eccentricity, eccentricAnomaly);
//
//            std::cout << "Time: " << time << ", Precision: 1e-" << i
//                << ", Eccentric Anomaly: " << eccentricAnomaly
//                << ", Calculated Mean Anomaly: " << calculatedMeanAnomaly << std::endl;
//        }
//    }
//
//    return 0;
//}



//// Функция для вычисления значения функции F(E) = M + e*sin(E) - E
//double calculateFunction(double M, double e, double E) {
//    return M + e * sin(E) - E;
//}
//
//// Функция для вычисления эксцентрической аномалии методом золотого сечения
//double calculateEccentricAnomaly(double M, double e, double accuracy) {
//    const double goldenRatio = (1 + sqrt(5)) / 2; // Золотое сечение
//
//    double a = 0;
//    double b = 2 * M_PI; // Пределы для эксцентрической аномалии
//    double c = b - (b - a) / goldenRatio;
//    double d = a + (b - a) / goldenRatio;
//
//    while (abs(b - a) > 2 * accuracy) {
//        if (calculateFunction(M, e, c) * calculateFunction(M, e, d) < 0) {
//            b = d;
//        }
//        else {
//            a = c;
//        }
//
//        c = b - (b - a) / goldenRatio;
//        d = a + (b - a) / goldenRatio;
//    }
//
//    // Найденная эксцентрическая аномалия
//    double E = (b + a) / 2;
//    return E;
//}


//double eccentricAnomaly(double meanAnomaly, double eccentricity, double tolerance) {
//    // Начальное приближение для эксцентрической аномалии
//    double E = meanAnomaly;
//
//    // Максимальное количество итераций
//    int maxIterations = 1000;
//
//    // Итерации методом Ньютона
//    for (int i = 0; i < maxIterations; ++i) {
//        double f = E - eccentricity * sin(E) - meanAnomaly;
//        double f_prime = 1.0 - eccentricity * cos(E);
//
//        // Метод Ньютона
//        double delta = f / f_prime;
//        E = E - delta;
//
//        // Проверка условия сходимости
//        if (fabs(delta) < tolerance || fabs(f) < tolerance) {
//            break;
//        }
//    }
//
//    return E;
//}
//
//int main() {
//    // Параметры орбиты
//    double meanAnomaly = 6.7895;
//    double eccentricity = 0.957;
//
//    // Задайте значения meanAnomaly и eccentricity, например:
//    // meanAnomaly = 30.0; // в градусах
//    // eccentricity = 0.5;
//
//    // Перевод в радианы
//    meanAnomaly = meanAnomaly * M_PI / 180.0;
//
//    // Вычисление эксцентрической аномалии с точностью от 10^(-3) до 10^(-10)
//    for (int i = 3; i <= 10; ++i) {
//        double tolerance = pow(10, -i);
//        double eccentricAnomalyResult = eccentricAnomaly(meanAnomaly, eccentricity, tolerance);
//
//        std::cout << "Tolerance: " << tolerance << ", Eccentric Anomaly: " << eccentricAnomalyResult << std::endl;
//    }
//
//    return 0;
//}