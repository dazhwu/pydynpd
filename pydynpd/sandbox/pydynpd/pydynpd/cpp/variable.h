//
// Created by Tiger on 5/9/2022.
//

#ifndef PYDYNPD_VARIABLE_H
#define PYDYNPD_VARIABLE_H

#include "dynpd.h"
#include "List_Variables.h"

class regular_variable{
public:
	string name;
	int lag;
	regular_variable(string, int);
};


class gmm_var{
public:
	string name;
	int min_lag;
	int max_lag;
	int lag;
	gmm_var(string, int, int, int);
};

vector<regular_variable> list_to_dep_indep_iv (List_Variables temp_list);
std::tuple< vector<gmm_var>, vector<gmm_var> > list_to_Gmm(List_Variables temp_list);

#endif //PYDYNPD_VARIABLE_H
