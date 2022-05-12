import math

import numpy as np
import pandas as pd
import scipy

import pydynpd.gmm_module as gmm_module
#from pydynpd.common_functions import get_first_diff_table, get_fod_table
from pydynpd.info import df_info, options_info

from pydynpd.panel_data import panel_data



class dynamic_panel_model(object):
    def __init__(self, identifiers, regression, basic_info, hansen, AR_test, options, list_parts):
        self.name = ''
        self.identifiers=identifiers
        self.T = basic_info.T
        self.N = basic_info.N
        self.num_obs=basic_info.num_obs
        self.num_indep=basic_info.num_indep
        self.num_instr=basic_info.num_instr
        self.dep=basic_info.dep
        self.indep=basic_info.indep
        self.form_regression_table(regression.regression_table)
        self.options = options
        self.max_obs=basic_info.max_obs
        self.min_obs=basic_info.min_obs
        self.avg_obs=basic_info.avg_obs
        self.hansen=hansen
        self.AR_list=AR_test
        self.command_str = list_parts[0] + '|' + list_parts[1]
        if list_parts[2] != '':
            self.command_str += '|' + list_parts[2]
        

        


    

    





    def calculate_MMSC_LU(self):
        self.MMSC_LU = {}
        log_n = math.log(self.num_obs)
        dif = self.num_instr - self.num_indep
        self.MMSC_LU["bic"] = self.hansen.test_value - (dif) * log_n
        self.MMSC_LU["hqic"] = self.hansen.test_value - dif * math.log(log_n) * 2.1
        self.MMSC_LU["aic"] = self.hansen.test_value - (dif) * 2

    def form_regression_table(self, the_result):

        var_names=self.indep


        coeff = the_result[:,0]
        std_err = the_result[:,1]
        z_value = the_result[:,2]
        p_value = the_result[:,3]
        sig = ['***' if p <= 0.001 else ('**' if p <= 0.01 else ('*' if p <= 0.05 else ' ')) for p in p_value]

        self.regression_table = pd.DataFrame(list(zip(var_names, coeff, std_err, z_value, p_value, sig)),
                                             columns=['variable', 'coefficient', 'std_err', 'z_value', 'p_value',
                                                      'sig'])

   