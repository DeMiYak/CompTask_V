#include"test.h"
#include"Task5.1.h"
#include"Task5.2.h"
#include"Task5.3.h"
#include<time.h>

namespace ublas = boost::numeric::ublas;

void testFormula() {
    Formula testOperator[] = {Formula("sin(pi)"),
                              Formula("cos(pi)"),
                              Formula("asin(1)"),
                              Formula("acos(~1)"),
                              Formula("atan(3^(0,5)/3)"),
                              Formula("tan(pi/4)"),
                              Formula("e^pi"),
                              Formula("ln(pi)"),
                              Formula("x + x^2 + 3"),
                              Formula("sgn(x)"),
                              Formula("~(x)"),
                              Formula("(cos(x)+x*(2+x^(~(x))))"),
                              Formula("ln([cos(x)+x*(2+x^(~(x)))])")};
    double x = -1;
    for (int i = 0; i < 13; ++i)
        cout << testOperator[i].Evaluate(x) << endl;
}

int main() {
    _wsetlocale(LC_ALL, L"russian");
    wcout.precision(16);
    cout.precision(16);

//    string sampleText = "[x-0,5]";
//    double x;
//
//    std::cin >> x;
//
//    Formula sampleFormula = Formula(sampleText);
//
//
//    cout << sampleFormula.Evaluate(x);
//
//    testAll();

//    Task5_1();
    clock_t timeStart = clock();
    Task5_2();
    wcout << L"Программа завершилась за " << (double)(clock() - timeStart)/CLOCKS_PER_SEC << L" (секунды)";
    return 0;
}
