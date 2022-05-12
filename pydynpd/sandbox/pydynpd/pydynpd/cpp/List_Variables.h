#pragma once
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "dynpd.h"
#include "Common_Functions.h"

using std::string;
using std::vector;

class List_Variables {
public:
  vector<string> names;
  vector<vector<int>> lags;
  vector<int> min_lags;
  vector<int> max_lags;
  vector<bool> adjustable_min_lags;
  vector<bool> adjustable_max_lags;
  List_Variables();
  // List_Variables(const vector<string> &);
  void append(string, vector<int>, bool, bool, vector<string> );
  string purge();
};
