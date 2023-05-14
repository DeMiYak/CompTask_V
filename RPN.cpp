#include"RPN.h"

using namespace std;

double sgn(double x) {
    if (x > 0) return 1;
    else if (x == 0) return 0;
    else return -1;
}

bool Formula::delimiter(char c) {
    return c == ' ';
}

bool Formula::isOperator(string c) {
    return _operator.find(c) != string::npos;
}

int Formula::priority(string operatorValue) {
    return
            operatorValue == "+" || operatorValue == "-" ? 1 :
            operatorValue == "*" || operatorValue == "/" ? 2 :
            operatorValue == "^" ? 3 :
            operatorValue == "~" || operatorValue == "sin" ||
            operatorValue == "cos" || operatorValue == "ln" ||
            operatorValue == "tan" || operatorValue == "asin" ||
            operatorValue == "acos" || operatorValue == "atan" ||
            operatorValue == "sgn" ? 4 :
            -1;
}

void Formula::ProcessOperator(stack<double> &valueStack, string operatorValue) {
    double rightOperand = valueStack.top();
    valueStack.pop();

    if (operatorValue == "~")
        valueStack.push(-rightOperand);

    else if (operatorValue == "sin")
        valueStack.push(sin(rightOperand));

    else if (operatorValue == "cos")
        valueStack.push(cos(rightOperand));

    else if (operatorValue == "ln")
        valueStack.push(log(rightOperand));

    else if (operatorValue == "tan")
        valueStack.push(tan(rightOperand));

    else if (operatorValue == "atan")
        valueStack.push(atan(rightOperand));

    else if (operatorValue == "asin")
        valueStack.push(asin(rightOperand));

    else if (operatorValue == "acos")
        valueStack.push(acos(rightOperand));

    else if(operatorValue == "sgn")
        valueStack.push(sgn(rightOperand));

    else {
        double leftOperand = valueStack.top();
        valueStack.pop();
        switch (operatorValue[0]) {
            case '+':
                valueStack.push(leftOperand + rightOperand);
                break;
            case '-':
                valueStack.push(leftOperand - rightOperand);
                break;
            case '*':
                valueStack.push(leftOperand * rightOperand);
                break;
            case '/':
                valueStack.push(leftOperand / rightOperand);
                break;
            case '^':
                valueStack.push(pow(leftOperand, rightOperand));
                break;
        }
    }
}

string Formula::ReplaceVariableWithDouble(double x) {
    string convertedVariable = to_string(x);
    string formulaCopy = _formula;
    for (size_t i = 0; i < formulaCopy.length(); i++) {
        if (formulaCopy[i] == 'x') {
            formulaCopy.insert(i, convertedVariable);
            if (formulaCopy[i] == '-') formulaCopy[i] = '~';
            i += convertedVariable.length();
            formulaCopy.erase(i, 1);
        } else if (formulaCopy[i] == 'e') {
            double e = 2.71828182845904;
            string stringE = to_string(e);
            formulaCopy.insert(i, stringE);
            i += stringE.length();
            formulaCopy.erase(i, 1);
        } else if (formulaCopy[i] == 'p') {
            double pi = acos(-1);
            string stringPi = to_string(pi);
            formulaCopy.insert(i, stringPi);
            i += stringPi.length();
            formulaCopy.erase(i, 2);

        }
    }
    return formulaCopy;
}

double Formula::CalculateFromString(const string &_string) {
    stack<double> valueStack;
    stack<string> operatorStack;
    for (size_t i = 0; i < _string.length(); i++) {

        if (!delimiter(_string[i])) {

            if (_string[i] == '(') {
                operatorStack.emplace("(");
            } else if (_string[i] == ')') {

                while (operatorStack.top() != "(")
                    ProcessOperator(valueStack, operatorStack.top()),
                            operatorStack.pop();

                operatorStack.pop();

            } else if (_string[i] == '[') {
                operatorStack.emplace("[");
            } else if (_string[i] == ']') {

                while (operatorStack.top() != "[")
                    ProcessOperator(valueStack, operatorStack.top()),
                            operatorStack.pop();

                operatorStack.pop();

                double temporaryValue = valueStack.top();
                valueStack.pop();
                valueStack.push(abs(temporaryValue));

            } else {
                string currentOperator(1, _string[i]);
                string compareLetter = "cslat";
                string compareOperation = "~-+*/^";
                if (compareLetter.find(currentOperator) != string::npos ||
                    compareOperation.find(currentOperator) != string::npos) {
                    if (compareLetter.find(currentOperator) != string::npos) {
                        i++;
                        while (_string[i] != '(') {
                            currentOperator += _string[i];
                            i++;
                        }

                    } else {
                        i++;
                    }

                    while (!operatorStack.empty() && priority(operatorStack.top()) >= priority(currentOperator))
                        ProcessOperator(valueStack, operatorStack.top()), operatorStack.pop();

                    operatorStack.emplace(currentOperator);
                    i--;
                } else {
                    string operandString(0, _string[i]);

                    while (i < _string.length() && '0' <= _string[i] && _string[i] <= '9')
                        operandString += _string[i++];

                    if (_string[i] == '.' || _string[i] == ',') operandString += _string[i++];

                    while (i < _string.length() && isalnum(_string[i])) {
                        operandString += _string[i++];
                    }

                    i--;

                    valueStack.push(atof(operandString.c_str()));
                }
            }
        }

    }


    while (!operatorStack.empty())
        ProcessOperator(valueStack, operatorStack.top()), operatorStack.pop();
    return valueStack.top();
}


double Formula::Evaluate(double x) {
    string stringModified = ReplaceVariableWithDouble(x);

    return CalculateFromString(stringModified);
}
