//
// Created by Пользователь on 08.05.2023.
//
#include"Task5.3.h"

double StretchedGaussValue(InterpolateQF IQFGauss, const double &startingPoint, const double &endingPoint) {
    if (IQFGauss._startingPoint != startingPoint || IQFGauss._endingPoint != endingPoint) {
        int size = IQFGauss._nodeNum;
        IQFGauss._startingPoint = startingPoint;
        IQFGauss._endingPoint = endingPoint;

        IQFGauss.RewriteIntegratorParameters(startingPoint, endingPoint, IQFGauss._segmentNum);

        double coefOne = (IQFGauss._endingPoint - IQFGauss._startingPoint) / 2;
        double coefTwo = coefOne + IQFGauss._startingPoint;

        for (int i = 0; i < size; ++i) {
            IQFGauss._polynomialRoot(i) = coefOne * IQFGauss._polynomialRoot(i) + coefTwo;
            IQFGauss._interpolationCoefficient(i) *= coefOne;
        }
    }
    return IQFGauss.IQFValue();
}

double CQFGauss(const InterpolateQF &IQFGauss, const double &startingPoint, const double &endingPoint,
                const int &CQFsegmentNum) {
    double step = (endingPoint - startingPoint) / CQFsegmentNum;
    double x = startingPoint;
    double CQFsum = 0;
    while (x < endingPoint) {
        CQFsum += StretchedGaussValue(IQFGauss, x, x + step);
        x += step;
    }
    return CQFsum;
}

void Task5_3() {
    wcout << L"Задание 5.3." << endl
          << L"Приближённое вычисление интеграла при помощи составной КФ Гаусса" << endl;
    wcout << endl << L"Вариант 7" << endl << L"Функции: " << endl;
    wcout << L"f(x) = sin(x)" << endl << "p(x) = |x - 0,5|";

    double startingPoint;
    double endingPoint;
    int nodeNum;
    int CQFsegmentNum;
    int segmentNum = 1e4;


    wcout << endl << L"Введите значения:" << endl;

    wcout << L"A - начало интегрирования: ";
    cin >> startingPoint;
    wcout << L"B - конец интегрирования: ";
    cin >> endingPoint;
    wcout << L"N - число узлов: ";
    cin >> nodeNum;
    wcout << L"m - число разбиений для составной КФ Гаусса: ";
    cin >> CQFsegmentNum;

    Formula GaussFormula = Formula("sin(x)*[x - 0,5]");
    Formula WeightFormula = Formula("1");

    wcout << endl << L"КФ Гаусса" << endl;

    Formula LegendrePolynomial = BuildLegendrePolynomial(nodeNum);

    InterpolateQF IQFGauss(WeightFormula, GaussFormula, LegendrePolynomial, segmentNum, nodeNum);

    double compoundQFGaussValue = CQFGauss(IQFGauss, startingPoint, endingPoint, CQFsegmentNum);
    double result = IQFPreciseResult(IQFGauss);

    IQFGauss.PrintData();

    Compare(compoundQFGaussValue, result);
}