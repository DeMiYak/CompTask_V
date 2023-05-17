//
// Created by Пользователь on 08.05.2023.
//
#include"Task5.3.h"

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

    double compoundQFGaussValue = IQFGauss.CQFGauss(startingPoint, endingPoint, CQFsegmentNum);

    IQFGauss.RewriteIntegratorParameters(startingPoint, endingPoint, segmentNum);

    double result = IQFPreciseResult(IQFGauss);

    IQFGauss.PrintData();

    Compare(compoundQFGaussValue, result);
}