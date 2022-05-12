#pragma once

#include <algorithm>
#include <climits>
#include <numeric>
#include <regex>
#include <string>
#include <vector>

#include "Common_Functions.h"
#include "List_Variables.h"
#include "info.h"
using std::regex;
using std::smatch;
using std::string;
using std::vector;

class Command {
public:
  int largest_T;

  string command_str;
  vector<string> cols;

  Command(int, string, vector<string>);

  void parse_command();
  void parse_dep_indep();
  void parse_gmm_iv();
  // bool variable_exists(string);
  bool parse_spaced_vars_range(string, List_Variables &);
  bool parse_spaced_vars_single(string, List_Variables &);
  bool parse_spaced_vars_auto(string, List_Variables &);
  void parse_spaced_vars(std::vector<std::string> list_vars,
                         List_Variables &dest_list);
  void parse_gmmStyle(vector<string> &);
  void parse_endo_pred(vector<string> &);
  void process_GMM(vector<string> &, int, int, string);
  void parse_endo_pred_general(vector<string> &, string, string, string, int);
  void parse_IV(vector<string> &);
  void parse_options();
};


std::tuple<List_Variables, List_Variables, List_Variables, List_Variables,
           vector<string>, struct model_options> process_command(int T_MAX, string commandstr, vector<string> df_col_names);
