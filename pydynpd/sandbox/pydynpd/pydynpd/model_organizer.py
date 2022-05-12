from pydynpd.command import command
from pydynpd.panel_data import panel_data
from pydynpd.variable import regular_variable


class list_models:
    def __init__(self):
        self.list_variables = []
        self.list_command_str = []


class model_oranizer(object):
    def __init__(self, user_command: command, pdata: panel_data):

        self._temp_list_indep = []
        self._temp_list_command = []

        self.ending_time = pdata.T
        first_list = []
        self._temp_list_indep.append(first_list)
        self._temp_list_command.append('')

        the_list = user_command._temp_part1_list
        for i in range(len(the_list.names)):
            var_name = the_list.names[i]
            lags = the_list.lags[i]
            num_lags = len(lags)
            last_lag = lags[0] + num_lags - 1
            for j in range(lags[0], last_lag + 1):
                new_var = regular_variable(var_name, j)
                self._add_item(new_var)

            if the_list.adjustable_max_lags[i] == True:
                self._explode_model(var_name, last_lag)

        self.models = list_models()
        for i in range(len(self._temp_list_indep)):
            new_variables = {}
            new_variables['dep_indep'] = self._temp_list_indep[i]
            new_variables['Dgmm'] = user_command.variables['Dgmm']
            new_variables['Lgmm'] = user_command.variables['Lgmm']
            new_variables['iv'] = user_command.variables['iv']
            self.models.list_variables.append(new_variables)
            self.models.list_command_str = self._temp_list_command

    def _add_item(self, new_var):
        for i in range(len(self._temp_list_indep)):

            self._temp_list_indep[i].append(new_var)
            if new_var.lag == 0:
                new_str = ' ' + new_var.name + ' '
            else:
                new_str = ' L' + str(new_var.lag) + '.' + new_var.name + ' '

            self._temp_list_command[i] += new_str

    def _explode_model(self, var_name, last_lag):
        new_list_variables = []
        new_list_command_str = []
        for i in range(len(self._temp_list_indep)):
            model = self._temp_list_indep[i]
            command_str = self._temp_list_command[i]

            for j in range(last_lag + 1, self.ending_time + 1):
                new_model = model.copy()
                new_str = str(command_str)

                for k in range(last_lag + 1, j + 1):
                    new_var = regular_variable(var_name, k)
                    new_model.append(new_var)
                    new_str += ' L' + str(k) + '.' + var_name + ' '

                new_list_variables.append(new_model)
                new_list_command_str.append(new_str)

        self._temp_list_indep += new_list_variables
        self._temp_list_command += new_list_command_str

    def _translate_to_command_str(self, variables, options):
        pass

    def check_model():
        pass
