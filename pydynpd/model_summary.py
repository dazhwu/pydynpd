import scipy.stats as st
from prettytable import PrettyTable


class model_summary(object):


    def print_summary(self, model):


        if model.options.steps == 2:
            str_steps = 'two-step '
        elif model.options.steps == 1:
            str_steps = 'one-step '
        else:
            str_steps = str(model.options.steps) + '-step '

        if model.options.level:
            str_gmm = 'system GMM'
        else:
            str_gmm = 'difference GMM'

        print('Dynamic panel-data estimation, ' + str_steps + str_gmm)
        print(self.basic_information(model))
        print(self.regression_table(model))
        print(self.test_results(model))

    def basic_information(self, model):
        basic_table = PrettyTable()
        basic_table.field_names = ["    ", "   ", "  "]
        basic_table.border = False
        basic_table.header = False
        basic_table.align = 'l'
        basic_table.add_row(['Group variable: ' + model.pdata._individual, ' ', 'Number of obs = ' + str(model.num_obs)])
        basic_table.add_row(['Time variable: ' + model.pdata._time, ' ', 'Number of groups = ' + str(model.N)])
        basic_table.add_row(['Number of instruments = ' + str(model.z_information.num_instr), ' ', ''])

        return (basic_table.get_string())

    def test_results(self, model):

        str_toprint = 'Hansen test of overid. restrictions: chi(' + str(model.hansen.df) + ') = ' + '{:.3f}'.format(
            model.hansen.test_value)
        str_toprint = str_toprint + ' Prob > Chi2 = ' + '{:.3f}'.format(model.hansen.p_value) + '\n'

        for i in range(len(model.AR_list)):
            AR = model.AR_list[i]
            P = st.norm.sf(abs(AR)) * 2
            str_toprint = str_toprint + 'Arellano-Bond test for AR(' + str(
                i + 1) + ') in first differences: z = ' + "{:.2f}".format(AR) + ' Pr > z =' + '{:.3f}'.format(P) + '\n'
        return (str_toprint)

    def regression_table(self, model):
        step=len(model.step_results)
        regression_result = model.step_results[step - 1]
        beta = regression_result.beta
        std_err = regression_result.std_err

        variables=model.variables

        dep_name = variables['dep_indep'][0].name
        var_names = []
        for i in range(1, len(variables['dep_indep'])):
            var_name = variables['dep_indep'][i].name
            var_lag = variables['dep_indep'][i].lag
            if (var_lag) >= 1:
                var_name = 'L' + str(var_lag) + '.' + var_name
            var_names.append(var_name)

        if model.options.level:
            var_names.append('_con')

        num_indep = len(var_names)
        r_table = PrettyTable()

        r_table.field_names = [dep_name, "coef.", "Corrected Std. Err.", "z", "P>|z|"]

        r_table.float_format = '.7'
        # , "z", "P>|z|", "[95% Conf. Interval]" ]
        for i in range(num_indep):
            var_name = var_names[i]
            coeff = beta[i, 0]
            stderr = std_err[i]
            z = coeff / stderr
            p = st.norm.sf(abs(z)) * 2
            r_table.add_row([var_name, coeff, stderr, z, p])

        return r_table.get_string()