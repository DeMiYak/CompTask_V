//
// Created by Пользователь on 08.05.2023.
//
#include"Task5.1.h"
#include "InterpolateQF.h"


Formula BuildLegendrePolynomial(const int &nodeNum = 0) {
    if (nodeNum <= 0) return {"1"};
    if (nodeNum == 1) return {"x"};
    else {
        Formula p0("1");
        Formula p1("x");
        Formula p2;
        for (int i = 2; i <= nodeNum; ++i) {
            p2 = to_string(double(2 * i - 1) / i) * p1 * "x" - to_string(double(i - 1) / i) * p0;
            p0 = p1;
            p1 = p2;
        }

        return p2;
    }

}

void Task5_2() {
    wcout << L"Задание 5.2." << endl
          << L"КФ Гаусса, Мелера, их узлы и коэффициценты\nВычисление интегралов с их помощью" << endl;
    wcout << endl << L"Вариант 7" << endl << L"Функции: " << endl;
    wcout << L"Гаусс: f(x) = (x + 0.8)/(x^2 + 1.2)^(1/2)" << endl << L"Мелер: f(x) = e^(2x)*x^2" << endl << "p(x) = 1";

    double startingPoint;
    double endingPoint;
    int nodeNum;
    int segmentNum = 1e4;

    wcout << endl << L"Введите значения:" << endl;

    wcout << L"A - начало интегрирования: ";
    cin >> startingPoint;
    wcout << L"B - конец интегрирования: ";
    cin >> endingPoint;
    wcout << L"N - число узлов: ";
    cin >> nodeNum;

    Formula GaussFormula = Formula("(x+0,8)/(x^2+1,2)^(1/2)");
    Formula MehlerFormula = Formula("e^(2*x)*x^2");
    Formula WeightFormula = Formula("1");

    Formula MehlerFormulaNew = MehlerFormula / "(1 - x^2)^(1/2)";

    wcout << endl << L"КФ Гаусса" << endl;

    Formula LegendrePolynomial = BuildLegendrePolynomial(nodeNum);

    InterpolateQF IQFGauss(WeightFormula, GaussFormula, LegendrePolynomial, segmentNum, nodeNum);

    IQFGauss.StretchIntegratorSegment(startingPoint, endingPoint);

    double resultGauss = IQFPreciseResult(IQFGauss);
    double GaussValue = IQFGauss.IQFValue();

    IQFGauss.PrintData();

    Compare(GaussValue, resultGauss);


    wcout << endl << L"КФ Мелера" << endl;

    InterpolateQF IQFMehler(WeightFormula, MehlerFormula, nodeNum);
    IQFMehler.RewriteIntegratorParameters(-1, 1, segmentNum);

    IQFMehler.MakeMehlerValue();

    double IQFMehlerValue = IQFMehler.IQFValue();

    Integrator MehlerIntegrator(IQFMehler.GetFormulaIntegrator().GetFormula() / "(1 - x^2)^(1/2)", -1, 1, segmentNum);

    IQFMehler.ModifyFormulaIntegrator(MehlerIntegrator);

    double resultMehler = IQFPreciseResult(IQFMehler);

    IQFMehler.PrintData();

    Compare(IQFMehlerValue, resultMehler);
}