//
// Created by Пользователь on 08.05.2023.
//
#include"Task5.1.h"
#include"InterpolateQF.h"


namespace ublas = boost::numeric::ublas;

double IQFPreciseResult(InterpolateQF &IQF) {
    Integrator resultIntegrator = Integrator(
            IQF.GetWeightIntegrator().GetFormula() * IQF.GetFormulaIntegrator().GetFormula(),
            IQF.GetFormulaIntegrator().GetStartingPoint(), IQF.GetFormulaIntegrator().GetEndingPoint(),
            IQF.GetFormulaIntegrator().GetSegmentNum());

    if (resultIntegrator.JSimpson().value == INFINITY)
        return resultIntegrator.JMiddleRectangle().value;

    return resultIntegrator.JSimpson().value;
}

bool CheckPrecisionIQF(InterpolateQF &IQF) {
    int nodeACT = 2 * IQF.GetNodeNum() - 1;
    string checkPrecisionString("x^" + to_string(nodeACT));
    ublas::vector<double> polynomialRoot = IQF.GetPolynomialRoot();
    ublas::vector<double> interpolationCoefficient = IQF.GetInterpolationCoefficient();
    Integrator weightIntegrator = IQF.GetWeightIntegrator();

    weightIntegrator.ModifyFormulaMultiply(checkPrecisionString);

    double sumCompare = 0;
    for (int i = 0; i < polynomialRoot.size(); ++i) {
        sumCompare += interpolationCoefficient(i) * pow(polynomialRoot(i), nodeACT);
    }

    double weight = weightIntegrator.JSimpson().value;

    cout << "weight: " << weight << "\nsumCompare: " << sumCompare << endl;

    if (abs(weight - sumCompare) < 1e-2)
        return true;

    return false;
}

void Compare(const double &IQFValue, const double &result) {
    double absoluteError = abs(IQFValue - result);
    double relativeError = absoluteError / abs(result) * 100;
    wcout << endl << endl << L"ИКФ: " << IQFValue << L"\nИнтеграл: " << result << L"\nАбсолютная погрешность: "
          << absoluteError;
    wcout << L"\nОтносительная погрешность: " << relativeError << "%" << endl;

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

    Formula funcFormula = Formula("sin(x)");
    Formula weightFormula = Formula("[x - 0,5]");

    Integrator resultIntegrator = Integrator(weightFormula * funcFormula, startingPoint, endingPoint, segmentNum);
    double result = resultIntegrator.JSimpson().value;

    InterpolateQF IQF(weightFormula, funcFormula, startingPoint, endingPoint, segmentNum, nodeNum);

    IQF.PrintData();

    Compare(IQF.IQFValue(), result);
}