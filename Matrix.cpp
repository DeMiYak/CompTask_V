//#include"Matrix.h"
//
//void Matrix::BuildMatrix() {
//    wcout << endl << L"Введите значения в матрице размера " << _dimVertical << "x" << _dimHorizontal;
//    wcout << L" слева-направо, сверху-вниз" << endl;
//
//    for (int i = 0; i < _dimVertical; ++i) {
//        for (int j = 0; j < _dimHorizontal; ++j) {
//            cin >> _matrix[i][j];
//        }
//    }
//}
//
//Matrix Matrix::operator*(const Matrix &matrix){
//    Matrix returnMatrix;
//    unsigned int horizontalSizeMatrix = matrix.GetDimHorizontal();
//
//    if(_dimHorizontal == matrix._dimVertical){
//        returnMatrix._matrix.resize(_dimVertical);
//        for (int i = 0; i < _dimVertical; ++i) {
//            for (int j = 0; j < horizontalSizeMatrix; ++j) {
//                double sum = 0;
//                for (int k = 0; k < _dimHorizontal; ++k) {
//                    sum += _matrix[i][k] + matrix._matrix[k][j];
//                }
//                returnMatrix._matrix[i].push_back(sum);
//            }
//        }
//    }
//
//    return returnMatrix;
//}
//
//Matrix Matrix::operator+(const Matrix &matrix){
//    Matrix returnMatrix;
//
//    if(_dimVertical == matrix._dimVertical && _dimHorizontal == matrix._dimHorizontal){
//        returnMatrix._matrix.resize(_dimVertical);
//        for (int i = 0; i < _dimVertical; ++i) {
//            for (int j = 0; j < _dimHorizontal; ++j) {
//                returnMatrix._matrix[i].push_back(_matrix[i][j] + matrix._matrix[i][j]);
//            }
//        }
//    }
//
//    return returnMatrix;
//}
//
//Matrix Matrix::operator-(const Matrix &matrix){
//    Matrix returnMatrix;
//
//    if(_dimVertical == matrix._dimVertical && _dimHorizontal == matrix._dimHorizontal){
//        returnMatrix._matrix.resize(_dimVertical);
//        for (int i = 0; i < _dimVertical; ++i) {
//            for (int j = 0; j < _dimHorizontal; ++j) {
//                returnMatrix._matrix[i].push_back(_matrix[i][j] - matrix._matrix[i][j]);
//            }
//        }
//    }
//
//    return returnMatrix;
//}