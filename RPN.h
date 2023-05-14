#ifndef RPN_H_INCLUDED
#define RPN_H_INCLUDED

#include<vector>
#include<string>
#include<cmath>
#include<cstdlib>
#include<stack>
#include<iostream>
#include<iomanip>

using namespace std;

class Formula {
public:

    Formula() { _formula = "0"; }

    Formula(string formula) { _formula = formula; }

    Formula(const Formula &formula) {
        _formula = formula._formula;
    };

    Formula operator=(const Formula& formula){
        _formula = formula._formula;
        return *this;
    }

    double Evaluate(double x);

    Formula operator*(const Formula &formulaTwo) {
        return *this * formulaTwo._formula;
    };

    Formula operator*(const string &formulaTwo) {
        return Formula('(' + this->_formula + ")*(" + formulaTwo + ')');
    };

    friend Formula operator*(const string &formulaOne, const Formula &formulaTwo) {
        return Formula('(' + formulaOne + ")*(" + formulaTwo._formula + ')');
    };

    Formula operator+(const Formula &formulaTwo){
        return *this + formulaTwo._formula;
    };

    Formula operator+(const string &formulaTwo){
        return Formula(this->_formula + "+" + formulaTwo);
    };

    friend Formula operator+(const string &formulaOne, const Formula &formulaTwo){
        return Formula(formulaOne + "+" + formulaTwo._formula);
    }

    Formula operator-(const Formula &formulaTwo){
        return *this - formulaTwo._formula;
    }

    Formula operator-(const string &formulaTwo){
        return Formula(this->_formula + "-(" + formulaTwo + ")");
    }

    friend Formula operator-(const string &formulaOne, const Formula &formulaTwo){
        return Formula(formulaOne + "-(" + formulaTwo._formula + ")");
    }

    string GetFormula() const{ return _formula; }

private:

    string _formula;

    const string _operator = "~+*^- sin cos ln tan asin acos atan sgn";

    bool delimiter(char c);

    bool isOperator(string c);

    int priority(string _operator);

    void ProcessOperator(stack<double> &_stack, string _operator);

    string ReplaceVariableWithDouble(double x);

    double CalculateFromString(const string &_string);
};

double sgn(double x);

#endif // RPN_H_INCLUDED
