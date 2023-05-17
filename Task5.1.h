//
// Created by Пользователь on 08.05.2023.
//

#ifndef COMPTASK_V_TASK5_1_H
#define COMPTASK_V_TASK5_1_H

#include"Integrator.h"
#include"RPN.h"
#include"invert_matrix.cpp"
#include"InterpolateQF.h"

void Task5_1();

double IQFPreciseResult(InterpolateQF &IQF);

bool CheckPrecisionIQF(InterpolateQF &IQF);

void Compare(const double &IQFValue, const double &result);

#endif //COMPTASK_V_TASK5_1_H
