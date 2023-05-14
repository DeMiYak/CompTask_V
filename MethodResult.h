//
// Created by Пользователь on 07.05.2023.
//

#ifndef COMPTASKIV_METHODRESULT_H
#define COMPTASKIV_METHODRESULT_H

struct MethodResult {
    MethodResult(double value, double error) {
        this->value = value;
        this->error = error;
    }

    MethodResult(double value){
        this->value = value;
    }

    double value;
    double error = NAN;
};

#endif //COMPTASKIV_METHODRESULT_H
