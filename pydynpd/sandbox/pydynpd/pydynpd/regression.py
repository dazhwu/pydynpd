import numpy as np
from scipy import stats
import math

# import pydynpd.ar_test as ar_test
import pydynpd.gmm_module as gmm_module
#import pydynpd.specification_tests as tests
import time
import warnings
from numpy.linalg import pinv
from pandas import DataFrame

from pydynpd.dynamic_panel_model import dynamic_panel_model

#from pydynpd.model_organizer import model_oranizer
from pydynpd.model_summary import model_summary
from pydynpd.panel_data import panel_data
from sys import exit

import sys

warnings.filterwarnings("ignore", category=RuntimeWarning)


class abond:

    def __init__(self, command_str, df: DataFrame, identifiers: list):

        if len(identifiers) != 2:
            print('two variables needed')
            exit()


        pdata = panel_data(df, identifiers)

        (temp_part1_list, temp_iv_list, DGMM_list, LGMM_list, List_parts, options) =   gmm_module.process_command(pdata.T, command_str, df.columns)

        pdata.export_data(df, temp_part1_list, temp_iv_list, DGMM_list, LGMM_list, options.timedumm)

        #(self.reg_table, self.vcov, self.SS, self.hansen, self.AR_test, self.z, self.y, self.x)=    gmm_module.prepare_data(pdata.data,temp_part1_list, temp_iv_list, DGMM_list, options, List_parts, [pdata.N, pdata.T], pdata.cols)
        (reg_table,  hansen, AR_test, basic_info)=gmm_module.prepare_data(pdata.data,temp_part1_list, temp_iv_list, DGMM_list, options, List_parts, [pdata.N, pdata.T], pdata.cols, pdata.col_timedumm)
        #print(reg_table.reg_table)
        #np.savetxt("z_tab.csv", self.z, delimiter=",")
        self.model=dynamic_panel_model(identifiers, reg_table, basic_info, hansen, AR_test, options, List_parts)
        #print(self.reg_table.regression_table)
        self.form_results(self.model)

 

    def iterative_GMM(self, model, _XZ, _XZ_t, _Zy, _Zy_t):
        current_step = 1
        converge = False
        while not converge:
            previous_step = current_step
            current_step += 1
            self.GMM(model, _XZ, _XZ_t, _Zy, _Zy_t, current_step)
            beta_current = model.step_results[current_step - 1].beta
            beta_previous = model.step_results[current_step - 2].beta
            for j in range(beta_current.shape[0]):
                temp = (beta_current[j] - beta_previous[j]) ** 2
                temp2 = (beta_previous[j]) ** 2
                if j == 0:
                    nom = temp
                    denom = temp2
                else:
                    nom += temp
                    denom += temp2
            crit = np.sqrt(nom / denom)

            if crit < 0.000001:
                converge = True
                model.options.steps = current_step

 

    def form_results(self, model):
        # step = len(model.step_results)
        # the_list = model.step_results[step - 1]
        if model.name != '':
            print(' ' + model.name)

        #model.form_regression_table()
        ms = model_summary()
        ms.print_summary(model)

        #self.models.append(model)  # results = {}

    def check_model(self, model):
        tbr = False
        num_ARs = len(model.AR_list)
        last_AR = model.AR_list[num_ARs - 1]

        if last_AR.P_value > 0.05:
            if model.hansen.p_value > 0.05 and model.hansen.p_value < 0.99999:
                return True
        else:
            return False
