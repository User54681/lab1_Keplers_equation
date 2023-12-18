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
double excentricToTrue(double E);
double iteration_method(double M0, double E, int i);
double radFind(double p, double teta);
double speedRad(double p, double teta);
double speedN(double p, double teta);
double fullSpeed(double p, double teta);
//void method_of_half_devision(double M0, double e, double epsilon, double pi, std::vector<double> M, int i);
//void golden_ratio_method(double M, double e, double epsilon, double pi);
//void newton_method(double M, double e, double epsilon);

int main()
{
    int i = 0;
    
    //iteration_method(i * n, i * n, i);

    for (int j = 0; j <= 360; j++) {
        double p = a * (1 - e * e);
        double teta = j * pi / 180;
        std::ofstream chart("chart.txt", std::ios::app);
        if (chart.is_open()) {
            chart << j*T/360  << "\t" << radFind(p, teta) << "\t" << speedRad(p, teta) << "\t" << speedN(p, teta) << "\t" << fullSpeed(p, teta) << "\n";
        }
        chart.close();
    }

    std::cout << "end";
    /*method_of_half_devision(M0, e, epsilon, pi, std::vector<double> M, int i);
    golden_ratio_method(M, e, epsilon, pi);
    newton_method(M, e, epsilon);*/
}

double radFind(double p, double teta)
{
    return p / (1 + e * cos(teta));
}

double speedRad(double p, double teta)
{
    return sqrt(mu / p) * e * sin(teta);
}

double speedN(double p, double teta)
{
    return sqrt(mu / p) * (1 + e * cos(teta));
}

double fullSpeed(double p, double teta)
{
    double rad = speedRad(p, teta);
    double n = speedN(p, teta);
    return sqrt(rad * rad + n * n);
}

double excentricToTrue(double E)
{
    if (atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2 > 0)
    {
        return atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;
    }
    else
    {
        return atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2 + 2 * pi;
    }
}

double iteration_method(double M0, double E0, int i) {
    while(i < T) {
        double E = M0 + e * sin(E0);
        double v = excentricToTrue(E);
        std::ofstream chartE("chartE.txt", std::ios::app);
        if (chartE.is_open()) {
            chartE << i << "\t" << M0 << "\t" << E << "\t" << v << "\n";
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
