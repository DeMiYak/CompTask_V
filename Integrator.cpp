#include <functional>
#include"Integrator.h"


using namespace std;

bool Integrator::ValidateValues(double startingPoint, double endingPoint, int segmentNum) {
    bool condOne = startingPoint > endingPoint;
    bool condTwo = segmentNum <= 0;
    if (condOne || condTwo) {
        if (condOne)
            wcout << L"Ошибка: конец меньше начала" << endl;
        if (condTwo)
            wcout << L"Ошибка: неположительное число разбиений" << endl;
        cout << endl;
        return false;
    }
    return true;

}


void Integrator::SetValues(double startingPoint, double endingPoint, int segmentNum) {
    if (ValidateValues(startingPoint, endingPoint, segmentNum)) {
        _startingPoint = startingPoint;
        _endingPoint = endingPoint;
        _segmentNum = segmentNum;
    }
}


MethodResult Integrator::JLeftRectangle(bool includeError) {
    double integralResult = 0;
    double step = (_endingPoint - _startingPoint) / _segmentNum;

    for (int stepNum = 0; stepNum < _segmentNum; stepNum++) {
        double x = _startingPoint + step * stepNum;
        integralResult += _formula.Evaluate(x);
    }

    if (includeError) {
        double moduleMax = findModuleMax(1);
        double errorCoefficient = pow(_endingPoint - _startingPoint, 2) / (2 * _segmentNum);
        return {step * integralResult, moduleMax * errorCoefficient};
    }
    return step*integralResult;
}

MethodResult Integrator::JRightRectangle(bool includeError) {
    double integralResult = 0;
    double step = (_endingPoint - _startingPoint) / _segmentNum;

    for (int stepNum = 0; stepNum < _segmentNum; stepNum++) {
        double x = _startingPoint + step * (stepNum + 0.5);
        integralResult += _formula.Evaluate(x);
    }

    if (includeError) {
        double moduleMax = findModuleMax(1);
        double errorCoefficient = pow(_endingPoint - _startingPoint, 2) / (2 * _segmentNum);
        return {step * integralResult, moduleMax * errorCoefficient};
    }
    return step*integralResult;
}

MethodResult Integrator::JMiddleRectangle(bool includeError) {
    double integralResult = 0;
    double step = (_endingPoint - _startingPoint) / _segmentNum;

    for (int stepNum = 0; stepNum < _segmentNum; stepNum++) {
        double x = _startingPoint + step * (stepNum + 0.5);
        integralResult += _formula.Evaluate(x);
    }

    if (includeError) {
        double moduleMax = findModuleMax(2);
        double errorCoefficient = pow(_endingPoint - _startingPoint, 3) / (24 * pow(_segmentNum, 2));
        return {step * integralResult, moduleMax * errorCoefficient};
    }
    return step*integralResult;
}

MethodResult Integrator::JTrapezoid(bool includeError) {
    double integralResult = 0;
    double step = (_endingPoint - _startingPoint) / _segmentNum;
    double x = _startingPoint;
    for (int stepNum = 0; stepNum < _segmentNum; stepNum++) {
        double y = _startingPoint + step * (stepNum + 1);
        integralResult += (_formula.Evaluate(x) + _formula.Evaluate(y)) / 2;

        x = y;
    }
    if(includeError) {
        double moduleMax = findModuleMax(2);
        double errorCoefficient = pow(_endingPoint - _startingPoint, 3) / (12 * pow(_segmentNum, 2));
        return {step * integralResult, moduleMax * errorCoefficient};
    }
    return step*integralResult;
}

MethodResult Integrator::JSimpson(bool includeError) {
    double integralResult = 0;
    double step = (_endingPoint - _startingPoint) / _segmentNum;
    double x = _startingPoint;
    for (int stepNum = 0; stepNum < _segmentNum; stepNum++) {
        double y = _startingPoint + step * (stepNum + 1);
        double betweenXY = (x + y) / 2;
        integralResult += (_formula.Evaluate(x) + 4 * _formula.Evaluate(betweenXY) + _formula.Evaluate(y)) / 6;

        x = y;
    }
    if(includeError){
        double moduleMax = findModuleMax(4);
        double errorCoefficient = pow(_endingPoint - _startingPoint, 5) / (2880 * pow(_segmentNum, 4));
        return {step * integralResult, moduleMax*errorCoefficient};
    }
    return step*integralResult;
}

MethodResult Integrator::JThreeEighths() {
    double integralResult = 0;
    double step = (_endingPoint - _startingPoint) / _segmentNum;
    double x = _startingPoint;

    for (int stepNum = 0; stepNum < _segmentNum; stepNum++) {
        double xOneThird = x + step / 3;
        double y = _startingPoint + step * (stepNum + 1);
        double yOneThird = y - step / 3;
        integralResult += (_formula.Evaluate(x) + 3 * (_formula.Evaluate(xOneThird) + _formula.Evaluate(yOneThird)) +
                           _formula.Evaluate(y)) / 8;

        x = y;
    }

    return step * integralResult;
}

double Integrator::ZeroDerivative(double x, double step){
    return _formula.Evaluate(x);
}

double Integrator::FirstDerivative (double x, double step) {
    if (step == -1)
        step = (_endingPoint - _startingPoint) / _segmentNum;
    return (_formula.Evaluate(x + step) - _formula.Evaluate(x - step)) / (2 * step);
}

double Integrator::SecondDerivative(double x, double step) {
    if (step == -1)
        step = (_endingPoint - _startingPoint) / _segmentNum;
    return (_formula.Evaluate(x + step) - 2 * _formula.Evaluate(x) + _formula.Evaluate(x - step)) / (step * step);
}

double Integrator::FourthDerivative(double x, double step) {
    if (step == -1)
        step = (_endingPoint - _startingPoint) / _segmentNum;
    return (_formula.Evaluate(x + 2 * step) - 4 * _formula.Evaluate(x + step) + 6 * _formula.Evaluate(x) -
            4 * _formula.Evaluate(x - step) + _formula.Evaluate(x - 2 * step)) / pow(step, 4);
}

double Integrator::findModuleMax(int derivativeIndex = -1) {
    int n = 10e4;
    double step = (_endingPoint - _startingPoint) / n;
    double Xmax = -1;
    switch (derivativeIndex) {
        case 1:
            for (int i = 0; i < n; ++i) {
                double x = _startingPoint + step * i;
                Xmax = max(Xmax, fabs(FirstDerivative(x, step)));
            }
            break;
        case 2:
            for (int i = 0; i < n; ++i) {
                double x = _startingPoint + step * i;
                Xmax = max(Xmax, fabs(SecondDerivative(x, step)));
            }
            break;
        case 4:
            for (int i = 0; i < n; ++i) {
                double x = _startingPoint + step * i;
                Xmax = max(Xmax, fabs(FourthDerivative(x, step)));
            }
            break;
        default:
            for (int i = 0; i < n; ++i) {
                double x = _startingPoint + step * i;
                Xmax = max(Xmax, fabs(_formula.Evaluate(x)));
            }
            break;
    }
    return Xmax;
}

//double Integrator::FirstDerivative(){
//    return pow
//}
//
//double Integrator::SecondDerivative(){
//    return pow(max(fabs(_endingPoint), fabs(_startingPoint)), 2)/2;
//}
//
//double Integrator::FourthDerivative(){
//    return 1;
//}

//double Integrator::FindErrorOfLeftAndRightRectangles(){
//
//    return pow(_endingPoint - _startingPoint, 2)/(2*_segmentNum);
//}
//
//double Integrator::FindErrorOfMiddleRectangles(){
//    return pow(_endingPoint - _startingPoint, 3)/(24*pow(_segmentNum, 2));
//}
//
//double Integrator::FindErrorOfTrapezoid(){
//
//}
//
//double Integrator::FindErrorOfSimpson(){
//
//}
//
//double Integrator::FindErrorOfThreeEighths(){
//    return pow(_endingPoint - _startingPoint, 5)/(pow(_segmentNum, 4)*2880);
//}
/** @brief
 * Интегрирует строку, заменяя элементарные функции на их интегралы
 *
 * Следует добавить:
 * проверку на пробел вначале и в конце строки
 *
 * @return
 */

//Formula Integrator::IntegrateFormula() {
//    string formulaString = _formula;
//    string delimiter = " -+";
//    vector<string> words{};
//
//    size_t position = 0;
//    while(formulaString[position++] == ' '){}
//    formulaString.erase(0, position);
//
//    while((position = formulaString.find(delimiter))!=string::npos){
//        words.push_back(formulaString.substr(0, position));
//        while(delimiter.find(formulaString[position++]) != string::npos){}
//        formulaString.erase(0, position);
//    }
//
//    for(const auto &word : words){
//        size_t position = 0;
//        if((position = word.find("x")) != string::npos){
//
//        }
//    }
//    return Formula(formulaString);
//}

