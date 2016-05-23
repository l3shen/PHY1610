// highest_value.h
#ifndef HIGHEST_VALUE_H
#define HIGHEST_VALUE_H
#include <rarray>
#include <iostream>
rarray<double,1> Highest_Value(rarray<double,1> Coefficients, rarray<int,1> &File_Index, int N_Highest);
void Print(rarray<int,1> File_Index, rarray<double,1> HighestValues, int N_Max);
#endif
