//
// Created by Tiger on 5/9/2022.
//

#include "variable.h"

regular_variable::regular_variable(string _name, int _lag) {
	name = _name;
	lag = _lag;
}

gmm_var::gmm_var(string _name, int _min, int _max, int _lag) {
	name = _name;
	min_lag = _min;
	max_lag = _max;
	lag = _lag;
}


vector<regular_variable> list_to_dep_indep_iv (List_Variables temp_list){
    vector<regular_variable> tbr;
    for(std::size_t i=0, max=temp_list.names.size(); i<max; ++i){
        string var_name = temp_list.names[i];
        vector<int> lags = temp_list.lags[i];
        int num_lags = lags.size();
        for (int j=lags[0]; j<lags[0] + num_lags; ++j){
            regular_variable new_var = regular_variable(var_name, j);
            tbr.push_back(new_var);

        }


    }

    return (tbr);
}

std::tuple< vector<gmm_var>, vector<gmm_var> > list_to_Gmm(List_Variables temp_list){
    vector<gmm_var> tbr_D,     tbr_L;
     for(std::size_t i=0, max=temp_list.names.size(); i<max; ++i){
        string var_name = temp_list.names[i];
        vector<int>  lags = temp_list.lags[i];
        int num_lags = lags.size();

        int min_lag=lags[0];

        int max_lag=min_lag + num_lags-1;
        gmm_var new_var = gmm_var(var_name, min_lag, max_lag, 0);
        tbr_D.push_back(new_var);

        int Lmin_lag = min_lag - 1;
        if (Lmin_lag <0)
            Lmin_lag =0;
        gmm_var temp_var = gmm_var(var_name, Lmin_lag, min_lag, 0); // # the 3rd argument (min_lag) not used
        tbr_L.push_back(temp_var);


        }

    return std::make_tuple(tbr_D, tbr_L);
    }
