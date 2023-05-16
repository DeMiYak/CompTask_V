//
// Created by Пользователь on 08.05.2023.
//
#include"Task5.1.h"

Formula BuildLegendrePolynomial(const int &nodeNum = 0){
    if(nodeNum <= 0) return {"1"};
    if(nodeNum == 1) return {"x"};
    else{
        Formula p0("1");
        Formula p1("x");
        Formula p2;
        for (int i = 2; i <= nodeNum ; ++i) {
            p2 = to_string(double(2*i - 1)/i)*p1*"x"-to_string(double(i-1)/i)*p0;
            p0 = p1;
            p1 = p2;
        }

        return p2;
    }

}

void MakeGaussCoefficient(InterpolateQF &IQFGauss){
    int size = IQFGauss._nodeNum;
    ublas::vector<double> interpolationCoefficient(size);
    Integrator QFHADPIntegrator(IQFGauss._QFHADP, -1, 1, IQFGauss._segmentNum);

    for (int i = 0; i < size; ++i) {
        double rootValue = IQFGauss._polynomialRoot(i);
        interpolationCoefficient(i) = 2/((1 - pow(rootValue, 2))*pow(QFHADPIntegrator.FirstDerivative(rootValue), 2));
    }
    IQFGauss._interpolationCoefficient = interpolationCoefficient;
}

void MakeMehlerValue(InterpolateQF &IQFMehler){
    int nodeNum = IQFMehler._nodeNum;
    const double pi = acos(-1);
    const double interpolationCoefficientConst = pi/nodeNum;
    for (int i = 0; i < nodeNum; ++i) {
        IQFMehler._interpolationCoefficient(i) = interpolationCoefficientConst;
        IQFMehler._polynomialRoot(i) = cos(double(2*i+1)/(2*nodeNum)*pi);
    }
}

void StretchIntegratorSegment(InterpolateQF &IQFGauss, const double &startingPoint, const double &endingPoint){
    if(IQFGauss._startingPoint != startingPoint || IQFGauss._endingPoint != endingPoint){
        int size = IQFGauss._nodeNum;
        IQFGauss._startingPoint = startingPoint;
        IQFGauss._endingPoint = endingPoint;

        IQFGauss.RewriteIntegratorParameters(startingPoint, endingPoint, IQFGauss._segmentNum);

        double coefOne = (IQFGauss._endingPoint - IQFGauss._startingPoint)/2;
        double coefTwo = coefOne + IQFGauss._startingPoint;

        for (int i = 0; i < size; ++i) {
            IQFGauss._polynomialRoot(i) = coefOne*IQFGauss._polynomialRoot(i) + coefTwo;
            IQFGauss._interpolationCoefficient(i) *= coefOne;
        }
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

    Formula MehlerFormulaNew = MehlerFormula/"(1 - x^2)^(1/2)";

    cout << MehlerFormulaNew.GetFormula();

    wcout << endl << L"КФ Гаусса" << endl;

    Formula LegendrePolynomial = BuildLegendrePolynomial(nodeNum);

    InterpolateQF IQFGauss(WeightFormula, GaussFormula, LegendrePolynomial, segmentNum, nodeNum);

    MakeGaussCoefficient(IQFGauss);

    StretchIntegratorSegment(IQFGauss, startingPoint, endingPoint);

    double resultGauss = IQFPreciselResult(IQFGauss);
    double GaussValue = IQFGauss.IQFValue();

    IQFGauss.PrintData();

    Compare(GaussValue, resultGauss);


    wcout << endl << L"КФ Мелера" << endl;

    InterpolateQF IQFMehler(WeightFormula, MehlerFormula, nodeNum);
    IQFMehler.RewriteIntegratorParameters(-1, 1, segmentNum);

    MakeMehlerValue(IQFMehler);

    double IQFMehlerValue = IQFMehler.IQFValue();

    Integrator MehlerIntegrator(IQFMehler.GetFormulaIntegrator().GetFormula() / "(1 - x^2)^(1/2)", -1, 1, segmentNum);

    IQFMehler.ModifyFormulaIntegrator(MehlerIntegrator);

    double resultMehler = IQFPreciselResult(IQFMehler);

    IQFMehler.PrintData();

    Compare(IQFMehlerValue, resultMehler);
}