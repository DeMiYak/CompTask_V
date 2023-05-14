//
// Created by Пользователь on 08.05.2023.
//
#include"Task5.1.h"

namespace ublas = boost::numeric::ublas;

ublas::vector<double> weightFunctionMoment(Integrator weightIntegrator, int nodeNum) {
    int nodeACT = 2 * nodeNum;
    ublas::vector<double> vectorMoment(nodeACT);
    for (int i = 0; i < nodeACT; ++i) {
        vectorMoment(i) = weightIntegrator.JSimpson().value;
        weightIntegrator.ModifyFormula("x");
    }
    return vectorMoment;

}

//ublas::vector<double> makeNodeVector(int nodeNum) {
//    std::vector<double> tempNodeVector;
//    ublas::vector<double> ublasNodeVector(nodeNum);
//
//    for (int i = 0; i < nodeNum; ++i) {
//        double tempNodeElement;
//        while (std::cin >> tempNodeElement) {
//            if (std::find(tempNodeVector.begin(), tempNodeVector.end(), tempNodeElement) == tempNodeVector.end()) {
//                tempNodeVector.push_back(tempNodeElement);
//                break;
//            } else { std::wcout << L"Ошибка: совпадающие числа. Пожалуйста, используйте другое число."; }
//        }
//    }
//
//    std::sort(tempNodeVector.begin(), tempNodeVector.end());
//
//    for (int i = 0; i < nodeNum; ++i) {
//        ublasNodeVector(i) = tempNodeVector[i];
//    }
//
//    return ublasNodeVector;
//}

ublas::vector<double> FindPolynomialCoefficient(const ublas::vector<double> weightFuncMoment) {
    const unsigned size = weightFuncMoment.size() / 2;
    ublas::matrix<double> matrixInverseMoment(size, size);
    ublas::vector<double> weightMomentHalf(size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrixInverseMoment(i, j) = weightFuncMoment(i + j);
        }
        weightMomentHalf(i) = -weightFuncMoment(i + size);
    }

    InvertMatrix(matrixInverseMoment, matrixInverseMoment);

    ublas::vector<double> polynomialCoefficient = prod(matrixInverseMoment, weightMomentHalf);

    return polynomialCoefficient;
}

Formula BuildQFHADA(const ublas::vector<double> &polynomialCoefficient) {
    int degree = polynomialCoefficient.size();
    Formula QFHADA = "x^" + to_string(degree);

    degree--;

    while (degree >= 0) {
        Formula temp = to_string(polynomialCoefficient(degree)) + "*x^" + to_string(degree);
        QFHADA = QFHADA + temp;
        degree--;
    }
    return QFHADA;
}

double CauchyBound(const ublas::vector<double> &vector) {
    double max = 0;
    for (int i = 0; i < vector.size(); ++i) {
        if (max < abs(vector(i))) max = abs(vector(i));
    }
    return 1 + max;
}

double secf(double beg, double end, Formula formula) {
    return beg - formula.Evaluate(beg) * (end - beg) / (formula.Evaluate(end) - formula.Evaluate(beg));
}

ublas::vector<double> FindPolynomialRoot(const ublas::vector<double> &polynomialCoefficient, Formula &QFHADA) {
    int rootNum = polynomialCoefficient.size();
    double foundRootBound = CauchyBound(polynomialCoefficient);
    int segmentNum = 1e5;
    double epsilon = 1e-12;
    double step = 2 * foundRootBound / segmentNum;
    double x = -foundRootBound;
    std::vector<double> polynomialRootStd;
    ublas::vector<double> polynomialRootUblasVector(rootNum);

    while (x < foundRootBound) {
        double a = x, b = x + step;
        double c = secf(a, b, QFHADA);
        if (QFHADA.Evaluate(a) * QFHADA.Evaluate(b) <= 0) {
            while (abs(c - b) >= epsilon) {
                a = b;
                b = c;
                c = secf(a, c, QFHADA);
            }
            polynomialRootStd.push_back(c);
        }
        x += step;
    }

    if (polynomialRootStd.size() == rootNum) {
        for (int i = 0; i < rootNum; ++i) {
            polynomialRootUblasVector(i) = polynomialRootStd[i];
        }
    } else { std::wcout << L"Ошибка: число получившихся корней меньше степени многочлена"; }

    return polynomialRootUblasVector;
}

ublas::vector<double> FindInterpolationCoefficient(const ublas::vector<double> &polynomialRootVector,
                                                   const ublas::vector<double> &weightFuncMoment) {
    int halfSize = weightFuncMoment.size() / 2;
    ublas::vector<double> weightVectorFirstHalf(halfSize);
    ublas::matrix<double> matrixWandermond(halfSize, halfSize);

    for (int i = 0; i < halfSize; ++i) {
        for (int j = 0; j < halfSize; ++j) {
            matrixWandermond(i, j) = pow(polynomialRootVector(j), i);
        }
        weightVectorFirstHalf(i) = weightFuncMoment(i);
    }

    InvertMatrix(matrixWandermond, matrixWandermond);

    ublas::vector<double>interpolationCoefficient = prod(matrixWandermond, weightVectorFirstHalf);

    return interpolationCoefficient;
}

double InterpolationQF(const ublas::vector<double> &interpolationCoefficient, const ublas::vector<double> &polynomialRootVector, Formula &QFHADA){
    int size = polynomialRootVector.size();
    ublas::vector<double> rootValueVector(size);

    for (int i = 0; i < size; ++i) {
        rootValueVector(i) = QFHADA.Evaluate(polynomialRootVector(i));
    }

    return dot_product(interpolationCoefficient, rootValueVector);
}

void Task5_1() {
    wcout << L"Задание 5.1." << endl
          << L"Приближённое вычисление интегралов при помощи квадратурных формул\nНаивысшей Алгебраической Степени Точности (КФ НАСТ)";
    wcout << endl << endl << L"Вариант 7";
    wcout << endl << L"Функции: " << endl << "f(x) = sin(x)" << endl << "p(x) = |x - 0.5|" << endl << endl;

    double startingPoint;
    double endingPoint;
    int nodeNum;
    int segmentNum = 1e4;

    wcout << L"Введите значения:" << endl;

    wcout << L"A - начало интегрирования: ";
    cin >> startingPoint;
    wcout << L"B - конец интегрирования: ";
    cin >> endingPoint;
    wcout << L"N - число узлов: ";
    cin >> nodeNum;

    Formula weightFormula = Formula("[x - 0,5]");
    Integrator weightIntegrator = Integrator(weightFormula, startingPoint, endingPoint, segmentNum);

    ublas::vector<double> weightFuncMoment = weightFunctionMoment(weightIntegrator, nodeNum);

    Formula funcFormula = Formula("sin(x)");

    Formula resultFormula = weightFormula * funcFormula;
    Integrator integratorResult(resultFormula, startingPoint, endingPoint, segmentNum);

    Formula formulaCompare = Formula("sin(x)*sgn(x - 0,5) + ((1 - 2*x)*cos(x)*sgn(x - 0,5))/2");
    double valueCompare = formulaCompare.Evaluate(endingPoint) - formulaCompare.Evaluate(startingPoint);


}