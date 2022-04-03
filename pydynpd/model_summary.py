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

        to_print = []
        # to_print.append(model.command_str)
        to_print.append('Dynamic panel-data estimation, ' + str_steps + str_gmm)
        to_print.append(self.basic_information(model))
        to_print.append(self.regression_table(model))
        to_print.append(self.test_results(model))
        for line in to_print:
            print(line)

    def basic_information(self, model):
        basic_table = PrettyTable()
        basic_table.field_names = ["    ", "   ", "  "]
        basic_table.border = False
        basic_table.header = False
        basic_table.align = 'l'
        basic_table.add_row(
            ['Group variable: ' + model.pdata._individual, ' ', 'Number of obs = ' + str(model.num_obs)])
        basic_table.add_row(['Time variable: ' + model.pdata._time, ' ', 'Min obs per group: ' + str(model.min_obs)])
        basic_table.add_row(['Number of instruments = ' + str(model.z_information.num_instr), ' ',
                             'Max obs per group: ' + str(model.max_obs)])
        basic_table.add_row(
            ['Number of groups = ' + str(model.N), ' ', 'Avg obs per group: ' + '{0:.2f}'.format(model.avg_obs)])

        return (basic_table.get_string())

    def test_results(self, model):

        str_toprint = 'Hansen test of overid. restrictions: chi(' + str(model.hansen.df) + ') = ' + '{:.3f}'.format(
            model.hansen.test_value)
        str_toprint = str_toprint + ' Prob > Chi2 = ' + '{:.3f}'.format(model.hansen.p_value) + '\n'

        for i in range(len(model.AR_list)):
            the_AR = model.AR_list[i]
            AR = the_AR.AR
            P_value = the_AR.P_value

            str_toprint = str_toprint + 'Arellano-Bond test for AR(' + str(
                i + 1) + ') in first differences: z = ' + "{:.2f}".format(AR) + ' Pr > z =' + '{:.3f}'.format(
                P_value) + '\n'
        return (str_toprint)

    def regression_table(self, model):

        dep_name = model.variables['dep_indep'][0].name

        r_table = PrettyTable()

        r_table.field_names = [dep_name, "coef.", "Corrected Std. Err.", "z", "P>|z|", " "]

        r_table.float_format = '.7'
        regression_table = model.regression_table
        # , "z", "P>|z|", "[95% Conf. Interval]" ]
        num_indep = len(regression_table.index)

        for i in range(num_indep):
            var_name = regression_table['variable'][i]
            coeff = regression_table['coefficient'][i]
            stderr = regression_table['std_err'][i]

            z = regression_table['z_value'][i]
            p = regression_table['p_value'][i]
            sig = regression_table['sig'][i]
            r_table.add_row([var_name, coeff, stderr, z, p, sig])

        return r_table.get_string()
