//
// Created by Пользователь on 08.05.2023.
//

#include"InterpolateQF.h"

namespace ublas = boost::numeric::ublas;

double InterpolateQF::StretchedGaussValue(InterpolateQF IQFGauss, const double &startingPoint, const double &endingPoint) {
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

double InterpolateQF::CQFGauss(const double &startingPoint, const double &endingPoint,
                               const int &CQFsegmentNum) {
    double step = (endingPoint - startingPoint) / CQFsegmentNum;
    double x = startingPoint;
    double CQFsum = 0;
    while (x < endingPoint) {
        CQFsum += StretchedGaussValue(*this, x, x + step);
        x += step;
    }
    return CQFsum;
}

void InterpolateQF::MakeMehlerValue() {
    int nodeNum = _nodeNum;
    const double pi = acos(-1);
    const double interpolationCoefficientConst = pi / nodeNum;
    for (int i = 0; i < nodeNum; ++i) {
        _interpolationCoefficient(i) = interpolationCoefficientConst;
        _polynomialRoot(i) = cos(double(2 * i + 1) / (2 * nodeNum) * pi);
    }
}

void InterpolateQF::StretchIntegratorSegment(const double &startingPoint, const double &endingPoint) {
    if (_startingPoint != startingPoint || _endingPoint != endingPoint) {
        int size = _nodeNum;
        _startingPoint = startingPoint;
        _endingPoint = endingPoint;

        RewriteIntegratorParameters(startingPoint, endingPoint, _segmentNum);

        double coefOne = (_endingPoint - _startingPoint) / 2;
        double coefTwo = coefOne + _startingPoint;

        for (int i = 0; i < size; ++i) {
            _polynomialRoot(i) = coefOne * _polynomialRoot(i) + coefTwo;
            _interpolationCoefficient(i) *= coefOne;
        }
    }
}


void InterpolateQF::MakeGaussCoefficient() {
    int size = _nodeNum;
    ublas::vector<double> interpolationCoefficient(size);
    Integrator QFHADPIntegrator(_QFHADP, -1, 1, _segmentNum);

    for (int i = 0; i < size; ++i) {
        double rootValue = _polynomialRoot(i);
        interpolationCoefficient(i) =
                2 / ((1 - pow(rootValue, 2)) * pow(QFHADPIntegrator.FirstDerivative(rootValue), 2));
    }
    _interpolationCoefficient = interpolationCoefficient;
}

void InterpolateQF::MakeWeightMoment(Integrator weightIntegrator) {
    int nodeACT = 2 * _nodeNum;
    ublas::vector<double> vectorMoment(nodeACT);
    for (int i = 0; i < nodeACT; ++i) {
        vectorMoment(i) = weightIntegrator.JSimpson().value;
        weightIntegrator.ModifyFormulaMultiply("x");
    }
    _weightFuncMoment = vectorMoment;

}

void InterpolateQF::MakePolynomialCoefficient() {
    const unsigned size = _weightFuncMoment.size() / 2;
    ublas::matrix<double> matrixInverseMoment(size, size);
    ublas::vector<double> weightMomentHalf(size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrixInverseMoment(i, j) = _weightFuncMoment(i + j);
//            cout << matrixInverseMoment(i,j) << "   ";
        }
        weightMomentHalf(i) = -_weightFuncMoment(i + size);
//        cout << endl;
    }
    InvertMatrix(matrixInverseMoment, matrixInverseMoment);

    ublas::vector<double> polynomialCoefficient = prod(matrixInverseMoment, weightMomentHalf);

    _polynomialCoefficient = polynomialCoefficient;
}

void  InterpolateQF::MakeQFHADP() {
    int degree = _polynomialCoefficient.size();
    string QFHADP = "x^" + to_string(degree);

    degree--;

    while (degree >= 0) {
        double coefficient = _polynomialCoefficient(degree);
        string temp = to_string(coefficient) + "*x^" + to_string(degree);
        if (coefficient >= 0) QFHADP = QFHADP + "+" + temp;
        else QFHADP = QFHADP + temp;
        degree--;
    }
    _QFHADP = Formula(QFHADP);
}

double InterpolateQF::CauchyBound(const ublas::vector<double> &vectorUblas) {
    double max = 0;
    for (int i = 0; i < vectorUblas.size(); ++i) {
        if (max < fabs(vectorUblas(i))) max = fabs(vectorUblas(i));
    }
    return 1 + max;
}

void InterpolateQF::MakePolynomialRoot() {
    int rootNum = _polynomialCoefficient.size();
    double foundRootBound = CauchyBound(_polynomialCoefficient);
    int segmentNum = _segmentNum;
    if(segmentNum%2 == 0) segmentNum += 1;
    double epsilon = 1e-12;
    double step = 2 * foundRootBound / segmentNum;
    double x = -foundRootBound;
    std::vector<double> polynomialRootStd;
    ublas::vector<double> polynomialRootUblasVector(rootNum);

    while (x < foundRootBound) {
        double a = x, b = x + step, c;

        if (_QFHADP.Evaluate(a) * _QFHADP.Evaluate(b) < 0) {
            while (b - a >= 2 * epsilon) {
                c = (a + b) / 2;
                if (_QFHADP.Evaluate(b) * _QFHADP.Evaluate(c) < 0) {
                    a = c;
                } else b = c;
            }
            polynomialRootStd.push_back((a+b)/2);
        }

        x += step;
    }

    if (polynomialRootStd.size() == rootNum) {
        for (int i = 0; i < rootNum; ++i) {
            polynomialRootUblasVector(i) = polynomialRootStd[i];
        }
    } else { std::wcout << L"Ошибка: число получившихся корней меньше степени многочлена"; }

    _polynomialRoot = polynomialRootUblasVector;
}

void InterpolateQF::MakeInterpolationCoefficient() {
    int halfSize = _weightFuncMoment.size() / 2;
    ublas::vector<double> weightVectorFirstHalf(halfSize);
    ublas::matrix<double> matrixWandermond(halfSize, halfSize);

    for (int i = 0; i < halfSize; ++i) {
        for (int j = 0; j < halfSize; ++j) {
            matrixWandermond(i, j) = pow(_polynomialRoot(j), i);
        }
        weightVectorFirstHalf(i) = _weightFuncMoment(i);
    }

    InvertMatrix(matrixWandermond, matrixWandermond);

    ublas::vector<double> interpolationCoefficient = prod(matrixWandermond, weightVectorFirstHalf);

    _interpolationCoefficient = interpolationCoefficient;
}

double InterpolateQF::IQFValue() {
    int size = _polynomialRoot.size();
    ublas::vector<double> rootValueVector(size, 0);

    for (int i = 0; i < size; ++i) {
        rootValueVector(i) = _formulaIntegrator.GetFormula().Evaluate(_polynomialRoot(i));
    }

    double sum = 0;

    for (int i = 0; i < size; ++i) {
        sum += _interpolationCoefficient(i) * rootValueVector(i);
    }

    return sum;
}

void InterpolateQF::PrintData() {
    std::vector<ublas::vector<double>> Data = {GetWeightFuncMoment(), GetPolynomialRoot(), GetInterpolationCoefficient()};
    std::vector<const wchar_t *> DataName = {L"Моменты весовой функции: ", L"Узлы интерполирования: ", L"Коэффициенты КФ НАСТ: "};
    for (int i = 0; i < 3; ++i) {
        wcout << DataName[i] << endl;
        for(const auto &element : Data[i])
            cout << element << "\t";
        cout << endl << endl;
    }
    wcout << L"Найденный ортогональный многочлен: " << endl;
    cout << _QFHADP.GetFormula();
}

InterpolateQF::InterpolateQF() = default;

InterpolateQF::InterpolateQF(const Integrator &weightIntegrator, const Integrator &formulaIntegrator) {
    _weightIntegrator = weightIntegrator;
    _formulaIntegrator = formulaIntegrator;
}

InterpolateQF::InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula) {
    _weightIntegrator = Integrator(weightFormula);
    _formulaIntegrator = Integrator(formulaFormula);
}

InterpolateQF::InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula, const int &nodeNum) {
    _weightIntegrator = Integrator(weightFormula);
    _formulaIntegrator = Integrator(formulaFormula);
    _nodeNum = nodeNum;
    _polynomialRoot = ublas::vector<double>(nodeNum, 0);
    _interpolationCoefficient = ublas::vector<double>(nodeNum, 0);
    _polynomialCoefficient = ublas::vector<double>(nodeNum, 0);
}

InterpolateQF::InterpolateQF(const string &weightString, const string &formulaString) {
    _weightIntegrator = Integrator(Formula(weightString));
    _formulaIntegrator = Integrator(Formula(formulaString));
}

InterpolateQF::InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula, const double &startingPoint,
                                           const double &endingPoint, const int &segmentNum, const int &nodeNum) {
    _weightIntegrator = Integrator(weightFormula);
    _formulaIntegrator = Integrator(formulaFormula);
    if (GoodInsert(startingPoint, endingPoint, segmentNum, nodeNum)){
        _nodeNum = nodeNum;
        RewriteIntegratorParameters(startingPoint, endingPoint, segmentNum);
    }

    MakeWeightMoment(_weightIntegrator);
    MakePolynomialCoefficient();
    MakeQFHADP();
    MakePolynomialRoot();
    MakeInterpolationCoefficient();
}

InterpolateQF::InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula, const Formula &QFHADP, const int &segmentNum, const int &nodeNum){
    _weightIntegrator = Integrator(weightFormula);
    _formulaIntegrator = Integrator(formulaFormula);
    _startingPoint = -1;
    _endingPoint = 1;
    RewriteSegmentNum(segmentNum);
    RewriteNodeNum(nodeNum);

    _polynomialCoefficient = ublas::vector<double>(nodeNum, 0);
    _QFHADP = QFHADP;
    MakePolynomialRoot();
    MakeGaussCoefficient();
}

void InterpolateQF::RewriteIntegratorParameters(const double &startingPoint, const double &endingPoint, const int &segmentNum){
    _startingPoint = startingPoint;
    _endingPoint = endingPoint;
    _segmentNum = segmentNum;
    _weightIntegrator.SetValues(_startingPoint, _endingPoint, _segmentNum);
    _formulaIntegrator.SetValues(_startingPoint, _endingPoint, _segmentNum);

}