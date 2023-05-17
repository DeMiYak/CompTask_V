//
// Created by Пользователь on 08.05.2023.
//

#ifndef COMPTASK_V_ROOTFINDER_H
#define COMPTASK_V_ROOTFINDER_H

#include "RPN.h"
#include"Integrator.h"
#include"invert_matrix.cpp"

class InterpolateQF {
public:
    InterpolateQF();

    InterpolateQF(const Integrator &weightIntegrator, const Integrator &formulaIntegrator);

    InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula);

    InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula, const int &nodeNum);

    InterpolateQF(const string &weightString, const string &formulaString);

    InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula, const double &startingPoint,
                  const double &endingPoint, const int &segmentNum, const int &nodeNum);

    InterpolateQF(const Formula &weightFormula, const Formula &formulaFormula, const Formula &QFHADP,
                  const int &segmentNum, const int &nodeNum);

    void MakeGaussCoefficient();

    void StretchIntegratorSegment(const double &startingPoint, const double &endingPoint);

    void MakeMehlerValue();

    double CQFGauss(const double &startingPoint, const double &endingPoint,
                    const int &CQFsegmentNum);


    void RewriteIntegratorParameters(const double &startingPoint, const double &endingPoint, const int &segmentNum);

    void RewriteNodeNum(const int &nodeNum) {
        _nodeNum = nodeNum;
        _weightFuncMoment.resize(2 * _nodeNum);
        _polynomialCoefficient.resize(_nodeNum);
        _polynomialRoot.resize(_nodeNum);
        _interpolationCoefficient.resize(_nodeNum);
    }

    void RewriteSegmentNum(const int &segmentNum) {
        if (segmentNum > 0) {
            _segmentNum = segmentNum;
            _weightIntegrator.RewriteSegmentNum(_segmentNum);
            _formulaIntegrator.RewriteSegmentNum(_segmentNum);
        }
    }

    double IQFValue();

    void PrintData();

    ublas::vector<double> GetWeightFuncMoment() { return _weightFuncMoment; }

    ublas::vector<double> GetPolynomialCoefficient() { return _polynomialCoefficient; }

    ublas::vector<double> GetPolynomialRoot() { return _polynomialRoot; }

    ublas::vector<double> GetInterpolationCoefficient() { return _interpolationCoefficient; }

    Integrator GetWeightIntegrator() { return _weightIntegrator; }

    Integrator GetFormulaIntegrator() { return _formulaIntegrator; }

    void ModifyFormulaIntegrator(const Integrator &formulaIntegrator) {
        _formulaIntegrator = formulaIntegrator;
        if (!(formulaIntegrator.GetSegmentNum() == _segmentNum &&
              formulaIntegrator.GetStartingPoint() == _startingPoint &&
              formulaIntegrator.GetEndingPoint() == _endingPoint))
            _formulaIntegrator.SetValues(_startingPoint, _endingPoint, _segmentNum);
    }

    int GetNodeNum() { return _nodeNum; }

private:

    Integrator _weightIntegrator;

    Integrator _formulaIntegrator;

    Formula _QFHADP;

    double _startingPoint = 0;

    double _endingPoint = 1;

    int _segmentNum = 10;

    int _nodeNum = 2;

    ublas::vector<double> _weightFuncMoment = ublas::vector<double>(2 * _nodeNum, 0);

    ublas::vector<double> _polynomialCoefficient = ublas::vector<double>(_nodeNum, 0);

    ublas::vector<double> _polynomialRoot = ublas::vector<double>(_nodeNum, 0);

    ublas::vector<double> _interpolationCoefficient = ublas::vector<double>(_nodeNum, 0);

    bool GoodInsert(const double &startingPoint, const double &endingPoint, const int &segmentNum, const int &nodeNum) {
        if ((segmentNum > 0) && (nodeNum > 1) && (startingPoint <= endingPoint))
            return true;
        return false;
    }

    void MakeWeightMoment(Integrator weightIntegrator);

    void MakePolynomialCoefficient();

    void MakeQFHADP();

    void MakePolynomialRoot();

    double CauchyBound(const ublas::vector<double> &vectorUblas);

    void MakeInterpolationCoefficient();

    double StretchedGaussValue(InterpolateQF IQFGauss, const double &startingPoint, const double &endingPoint);
};

#endif //COMPTASK_V_ROOTFINDER_H
