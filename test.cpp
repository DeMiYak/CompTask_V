#include"test.h"

using namespace boost::numeric::ublas;

void testOne() {
    std::cout << "Test One: " << std::endl;
    matrix<std::complex<double> > m(3, 3);
    for (unsigned i = 0; i < m.size1(); ++i)
        for (unsigned j = 0; j < m.size2(); ++j)
            m(i, j) = std::complex<double>(3 * i + j, 3 * i + j);

    std::cout << -m << std::endl;
    std::cout << conj(m) << std::endl;
    std::cout << real(m) << std::endl;
    std::cout << imag(m) << std::endl;
    std::cout << trans(m) << std::endl;
    std::cout << herm(m) << std::endl << std::endl;
}

void testTwo() {
    std::cout << "Test Two: " << std::endl;

    matrix<double> m1(3, 3), m2(3, 3);
    for (unsigned i = 0; i < std::min(m1.size1(), m2.size1()); ++i)
        for (unsigned j = 0; j < std::min(m1.size2(), m2.size2()); ++j)
            m1(i, j) = m2(i, j) = 3 * i + j;

    std::cout << m1 + m2 << std::endl;
    std::cout << m1 - m2 << std::endl << std::endl;
}

void testThree() {
    std::cout << "Test Three: " << std::endl;

    matrix<double> m(3, 3);
    vector<double> v(3);
    for (unsigned i = 0; i < std::min(m.size1(), v.size()); ++i) {
        for (unsigned j = 0; j < m.size2(); ++j)
            m(i, j) = 3 * i + j;
        v(i) = i;
    }
    std::cout << m << std::endl;
    std::cout << v << std::endl;
    std::cout << prod(m, v) << std::endl;
    std::cout << prod(v, m) << std::endl;
    std::cout << prod(m, m) << std::endl;
}

void testFour(){
    std::cout << "Test Four: " << std::endl;

    std::vector<std::vector<double>> inputTemporary = {{1, 1, 1, 1}, {2, 3, 4, 5}, {4, 9, 16, 25}, {8, 27, 64, 125}};
    matrix<double> input(4,4);
    matrix<double> inputInverse(4, 4);
    for (unsigned i = 0; i < input.size1(); ++i) {
        for (unsigned j = 0; j < input.size2(); ++j)
            input(i, j) = inputTemporary[i][j];
    }
    InvertMatrix<double>(input, inputInverse);
    InvertMatrix(input, input);
    std::cout << input << std::endl << inputInverse << std::endl << prod(input, inputInverse) << std::endl << prod(inputInverse, input);
}

void testAll () {
    testOne();

    testTwo();

    testThree();

    testFour();
}