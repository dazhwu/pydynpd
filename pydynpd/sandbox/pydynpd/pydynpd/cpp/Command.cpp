#include "Command.h"

List_Variables temp_part1_list;
List_Variables temp_iv_list;
List_Variables LGMM_list;
List_Variables DGMM_list;
string part_1, part_2, part_3;
struct model_options options;


std::tuple<List_Variables, List_Variables, List_Variables, List_Variables, vector<string>, struct model_options> process_command(int T_MAX, string commandstr, vector<string> df_col_names)
{

	temp_iv_list = List_Variables();
	temp_part1_list = List_Variables();
	LGMM_list = List_Variables();
	DGMM_list = List_Variables();
	options = model_options();
	part_1 = "";
	part_2 = "";
	part_3 = "";

	Command new_command = Command(T_MAX, commandstr, df_col_names);

	vector<string> combined{part_1, part_2, part_3};


	return std::make_tuple(temp_part1_list, temp_iv_list, DGMM_list, LGMM_list,
						   combined, options);
}

Command::Command(int T_MAX, string commandstr, vector<string> df_col_names)
{
	command_str = commandstr;
	largest_T = T_MAX;
	cols = df_col_names;

	parse_command();
}

void Command::parse_command()
{

	vector<string> parts = splitString(command_str, '|');

	if (parts.size() <= 1)
		throw std::invalid_argument("There should be at least two parts in command string");

	if (parts.size() > 3)
		throw std::invalid_argument("too many parts in command string");

	if (parts.size() == 3)
	{
		part_3 = parts[2];
	}
	else
	{
		part_3 = "";
		// self.options = options_info()
	}

	part_1 = parts[0];

	parse_dep_indep();

	part_2 = parts[1];
	parse_gmm_iv();
}

void Command::parse_dep_indep()
{

	vector<string> list_vars = splitString(part_1, ' ');
	parse_spaced_vars(list_vars, temp_part1_list);
}

void Command::parse_gmm_iv()
{
	vector<string> matching_parts;

	parse_gmmStyle(matching_parts);
	parse_endo_pred(matching_parts);
	parse_IV(matching_parts);

	//	part2_cpy = part_2
	//	for part in matching_parts :
	// part2_cpy = part2_cpy.replace(part, "")

	//	if len(part2_cpy.strip()) > 0:
	// print(part2_cpy.strip() + ': invalid GMM or IV statement')
	//	exit();
}

bool Command::parse_spaced_vars_range(string var, List_Variables &dest_list)
{
	string lag_range = "^L[(]([0-9]{1,})[:]([0-9]{1,})[)][.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$";
	regex pt_range(lag_range);
	smatch match_groups_multiple;
	bool isMatch = regex_match(var, match_groups_multiple, pt_range);
	if (isMatch)
	{
		int LB = std::stoi(match_groups_multiple[1]);
		int UB = std::stoi(match_groups_multiple[2]);
		if (LB > UB)
		{
			int temp = LB;
			LB = UB;
			UB = temp;
		}
		string name = match_groups_multiple[3];

		// if (!variable_exists(name))
		//     throw std::invalid_argument(var + ": variable " + name + " is invalid");

		vector<int> ivec(UB - LB + 1);
		std::iota(ivec.begin(), ivec.end(), LB);
		dest_list.append(name, ivec, false, false, cols);
		return true;
	}
	else
	{
		return false;
	}
}

bool Command::parse_spaced_vars_single(string var, List_Variables &dest_list)
{
	string lag_single = "^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$";
	regex pt_single(lag_single);
	smatch match_groups_single;
	bool isMatch = regex_match(var, match_groups_single, pt_single);
	if (isMatch)
	{
		int lag = stoi(match_groups_single[1]);
		string name = match_groups_single[2];

		vector<int> ivec(1);
		ivec[0] = lag;
		dest_list.append(name, ivec, false, false, cols);
		return true;
	}
	else
		return false;
}

bool Command::parse_spaced_vars_auto(string var, List_Variables &dest_list)
{
	string lag_auto = "^L[(]([0-9]{1,})[:]([?])[)][.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$";
	regex pt_auto(lag_auto);
	smatch match_groups_auto;
	bool isMatch = regex_match(var, match_groups_auto, pt_auto);
	if (isMatch)
	{
		int LB = stoi(match_groups_auto[1]);
		string name = match_groups_auto[3];

		options.beginner = true;
		vector<int> ivec(1);
		ivec[0] = LB;
		dest_list.append(name, ivec, false, true, cols);
		return true;
	}
	else
		return false;
}

// bool Command::variable_exists(string name) {
//     if (std::find(cols.begin(), cols.end(), name) == cols.end()) {
//         return false;
//     }
//     return true;
// }
void Command::parse_spaced_vars(vector<string> list_vars, List_Variables &dest_list)
{

	bool isMatch;
	for (string var : list_vars)
	{

		isMatch = parse_spaced_vars_range(var, dest_list);
		if (!isMatch)
		{
			isMatch = parse_spaced_vars_single(var, dest_list);

			if (!isMatch)
			{
				isMatch = parse_spaced_vars_single(var, dest_list);
				if (!isMatch)
					isMatch = parse_spaced_vars_auto(var, dest_list);
				if (!isMatch)
				{
					vector<int> ivec(1);
					ivec[0] = 0;
					dest_list.append(var, ivec, false, false, cols);
				}
			}
		}
	}
}

void Command::parse_gmmStyle(vector<string> &matching_parts)
{
	smatch match;
	string temp = part_2;
	regex r("gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{0,}[:][ ]{0,}(?:(?:[.])|(?:[0-9]{1,}))[ ]{0,}[)]");
	int min_lag, max_lag;
	// int i = 1;
	while (regex_search(temp, match, r))
	{
		string part = match.str(0);
		matching_parts.push_back(part);
		smatch match_groups_multiple;
		regex prog_1("^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{0,}[:][ ]{0,}((?:[.])|(?:[0-9]{1,}))[ ]{0,}[)]$");
		regex_match(part, match_groups_multiple, prog_1);

		vector<string> vars = splitString(match_groups_multiple[1], ' ');

		min_lag = stoi(match_groups_multiple[2]);
		if (match_groups_multiple[3] == '.')
			max_lag = largest_T;
		else
			max_lag = stoi(match_groups_multiple[3]);

		process_GMM(vars, min_lag, max_lag, part);

		temp = match.suffix().str();
	}
}

void Command::parse_endo_pred_general(vector<string> &matching_parts, string temp, string pattern, string group_str, int min_lag)
{
	smatch match;
	regex r(pattern);

	while (regex_search(temp, match, r))
	{
		string part = match.str(0);
		matching_parts.push_back(part);
		smatch match_groups_multiple;
		regex prog_1(group_str);
		regex_match(part, match_groups_multiple, prog_1);

		vector<string> vars = splitString(match_groups_multiple[1], ' ');
		process_GMM(vars, min_lag, largest_T, part);
		temp = match.suffix().str();
	}
}

void Command::parse_endo_pred(vector<string> &matching_parts)
{

	string temp = part_2;
	parse_endo_pred_general(matching_parts, temp, "endo[(][a-zA-Z_0-9 ]{1,}[)]", "^endo[(]([a-zA-Z_0-9 ]{1,})[)]$", 2);
	parse_endo_pred_general(matching_parts, temp, "pred[(][a-zA-Z_0-9 ]{1,}[)]", "^pred[(]([a-zA-Z_0-9 ]{1,})[)]$", 2);
}

void Command::process_GMM(vector<string> &vars, int min_lag, int max_lag, string part)
{
	if (min_lag > max_lag)
		throw std::invalid_argument(part + ": minimum lag cannot be greater than maximum lag");

	if (min_lag < 0)
		throw std::invalid_argument(part + ": lags must be non-negative");

	if (vars.size() == 0)
		throw std::invalid_argument(part + ": no variable is included");

	for (string var : vars)
	{
		if (!variable_exists(var, cols))
			throw std::invalid_argument(part + ": " + var + " does not exist");
		if (variable_exists(var, DGMM_list.names))
			throw std::invalid_argument(part + ": " + var + " cannot be declared in part 2 for twice or more");

		vector<int> ivec(max_lag - min_lag + 1);
		std::iota(ivec.begin(), ivec.end(), min_lag);

		DGMM_list.append(var, ivec, false, false, cols);
		int Lmin_lag;
		if (min_lag - 1 < 0)
			Lmin_lag = 0;
		else
			Lmin_lag = min_lag - 1;

		vector<int> Livec(max_lag - Lmin_lag + 1);
		std::iota(Livec.begin(), Livec.end(), Lmin_lag);
		LGMM_list.append(var, Livec, false, false, cols);
	}
}

void Command::parse_IV(vector<string> &matching_parts)
{
	smatch match;
	string temp = part_2;
	regex r("iv[(].{1,}[)]");

	while (regex_search(temp, match, r))
	{
		string part = match.str(0);
		matching_parts.push_back(part);
		smatch match_groups_multiple;
		regex prog_1("^iv[(](.{1,})[)]$");
		regex_match(part, match_groups_multiple, prog_1);

		vector<string> vars = splitString(match_groups_multiple[1], ' ');

		parse_spaced_vars(vars, temp_iv_list);

		temp = match.suffix().str();
	}
}

void Command::parse_options()
{
	vector<string> list_options = splitString(part_3, ' ');

	for (string option : list_options)
	{
		if (option == "onestep")
			options.steps = 1;

		else if (option == "iterated")
			options.steps = 1000;

		else if (option == "nolevel")
			options.level = false;

		else if (option == "hqic")
			options.mmsc = "hqic";

		else if (option == "fod")
			options.transformation = "fod";

		else if (option == "timedumm")
			options.timedumm = true;

		else if (option == "collapse")
			options.collapse = true;

		else
			throw std::invalid_argument(option + ": is not a valid option");
	}
}
