//
// Created by Tiger on 5/4/2022.
//

#ifndef UNTITLED_COMMON_FUNCTIONS_H
#define UNTITLED_COMMON_FUNCTIONS_H

#include "dynpd.h"

vector<bool> row_has_nan(const RowMatrixXd &x);

RowMatrixXd common_inv(Ref < RowMatrixXd > );

double standard_normalCDF(double);

bool areConsecutive(int*, int);
//void tokenize(string const &str, const char* delim,
//            std::vector<std::string> &out);

vector<string> splitString(string, char);
bool variable_exists(string,  vector<string>&);
int getIndex(vector<string> &, string );

void lag(Ref<RowMatrixXd> mat, Ref<RowMatrixXd> lagged, int N, int lag_number, double fill);

RowMatrixXd get_first_diff_table(Ref<RowMatrixXd> ori_arr, int N);

RowMatrixXd	get_fod_table(Ref<RowMatrixXd> ori_arr, int N);

void saveData(string fileName, RowMatrixXd  matrix);

#endif //UNTITLED_COMMON_FUNCTIONS_H
