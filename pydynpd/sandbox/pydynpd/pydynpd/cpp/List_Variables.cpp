#include "List_Variables.h"

List_Variables::List_Variables(){
	
}


void List_Variables::append(string name, vector<int> new_lag_vector, bool min_adj_lag, bool max_adj_lag, vector<string> cols){
	
	if (!variable_exists(name, cols))
		throw std::invalid_argument("Variable: " + name + " does not exist");

	std::vector<string>::iterator existing_variable =std::find(names.begin(), names.end(), name);
	if (existing_variable==names.end()){
		names.push_back(name);
		lags.push_back(new_lag_vector);
		adjustable_min_lags.push_back(min_adj_lag);
		adjustable_max_lags.push_back(max_adj_lag);
	}else{
		int index=existing_variable-names.begin();
		lags[index].insert(lags[index].end(), new_lag_vector.begin(), new_lag_vector.end());
		if (min_adj_lag)
			adjustable_min_lags[index] = min_adj_lag;
		
		if (max_adj_lag)
			adjustable_max_lags[index]=max_adj_lag;
	}
}

string List_Variables::purge(){
	std::vector< vector<int> >::iterator vec_lags;
	for (vec_lags=lags.begin(); vec_lags!=lags.end(); vec_lags++){
		int index=vec_lags-lags.begin();
		std::sort( lags[index].begin(), lags[index].end() );
		lags[index].erase( std::unique( lags[index].begin(), lags[index].end() ), lags[index].end() );
		if (!areConsecutive(&lags[index][0], lags[index].size())){
			return "variable " + names[index] + " has gaps";
		}
	}
	return "";
}


