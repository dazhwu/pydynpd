import pydynpd.variable
from pydynpd.variable import gmm_var, regular_variable
from pydynpd.info import options_info
import sys
from sys import exit

import re

class temp_list:
    def __init__(self, cols):
        self.names = []
        self.lags=[]
        self.cols=cols

    def insert(self, name, lags):
        if name not in self.cols:
            return -1

        if name not in self.names:
            self.names.append(name)
            self.lags.append(lags)
        else:
            the_index=self.names.index(name)
            self.lags[the_index]+=lags

        return 0

    def purge(self):
        self.lags = [sorted(list(set(the_list))) for the_list in self.lags]
    
    def check_contiguous(self):
        for i in range(len(self.names)):
            the_list=self.lags[i]
            if the_list!=list(range(min(the_list),max(the_list)+1)):
                print( 'variable ' + self.name[i] + ' has gaps')
                exit()


class command(object):

    def __init__(self, command_str, df_col_names):
        self.command_str=command_str
        self.cols=df_col_names
        self.temp_part1_list=temp_list(df_col_names)
        self.temp_iv_list = temp_list(df_col_names)
        self.variables=None
        self.options=None

        self.list_GMM=[]

        self.parse_command()

    def parse_command(self):
        command_str=self.command_str
        parts = command_str.split('|')
        if len(parts) <= 1:
            print('There should be at least two parts in command string')
            exit()

        part_1 = parts[0]
        self.parse_dep_indep(part_1)

        part_2 = parts[1]
        self.parse_gmm_iv(part_2)

        if len(parts) == 3:
            part_3 = parts[2]
            self.options = self.parse_options(part_3)
        else:
            self.options = options_info()

        self.check_dep_indep()
        self.check_GMM()
        self.check_iv()
        self.check_three_lists()

        self.variables = {}
        self.variables['dep_indep'] = self.tbr_list(self.temp_part1_list)
        self.variables['gmm'] = self.list_GMM
        self.variables['iv'] = self.tbr_list(self.temp_iv_list)


    def parse_spaced_vars(self,list_vars, dest_list):

        prog_1 = re.compile('^L\(([0-9]{1,})[:]([0-9]{1,})\)[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
        prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')

        for var in list_vars:
            match_groups_multiple = prog_1.match(var)
            match_groups_single = prog_2.match(var)

            if match_groups_multiple:
                LB=int(match_groups_multiple.group(1))
                UB=int(match_groups_multiple.group(2))
                name=match_groups_multiple.group(3)
                ret= dest_list.insert(name, list(range(min(LB, UB), max(LB,UB)+1)))

            elif prog_2.match(var):
                lag = int(match_groups_multiple.group(1))
                name = match_groups_multiple.group(2)
                ret= dest_list.insert(name, [lag])
            else:
                name=var
                ret= dest_list.insert(name, [0])

            if ret==-1:
                return name

        return ''

    def parse_dep_indep(self, part_1):

        list_vars = part_1.split()

        ret = self.parse_spaced_vars(list_vars, self.temp_part1_list)

        if ret != '':
            print(part_1 + ':  variable ' + ret + ' does not exist')
            exit()

    def parse_gmm_iv(self,part_2):

        matching_parts = []           

        self.parse_gmmStyle(matching_parts,part_2)
        self.parse_endo_pred(matching_parts,part_2)
        self.parse_IV(matching_parts,part_2)

        part2_cpy = part_2
        for part in matching_parts:
            part2_cpy = part2_cpy.replace(part, '')

        if len(part2_cpy.strip()) > 0:
            print(part2_cpy.strip() + ': invalid GMM or IV statement')
            exit()



    def parse_gmmStyle(self, matching_parts, part_2):

        gmm_search_parts = re.findall(
            'gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{0,}[:][ ]{0,}(?:(?:[.])|(?:[0-9]{1,}))[ ]{0,}[)]', part_2)
        prog_1 = re.compile(
            '^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{0,}[:][ ]{0,}((?:[.])|(?:[0-9]{1,}))[ ]{0,}[)]$')

        for part in gmm_search_parts:

            matching_parts.append(part)
            match_groups_multiple = prog_1.match(part)
            vars = match_groups_multiple.group(1).split()

            min_lag = int(match_groups_multiple.group(2))
            if match_groups_multiple.group(3) == '.':
                max_lag = sys.maxsize
            else:
                max_lag = int(match_groups_multiple.group(3))

            self.process_GMM(vars, min_lag, max_lag,  part)


    def parse_endo_pred(self, matching_parts, part_2):

        gmm_search_parts = re.findall('endo[(][a-zA-Z_0-9 ]{1,}[)]', part_2)
        prog_1 = re.compile('^endo[(]([a-zA-Z_0-9 ]{1,})[)]$')

        for part in gmm_search_parts:
            # prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
            matching_parts.append(part)
            match_groups_multiple = prog_1.match(part)

            vars = match_groups_multiple.group(1).split()
            min_lag = 2
            max_lag = sys.maxsize

            self.process_GMM(vars, min_lag, max_lag,  part)

        gmm_search_parts = re.findall('pred[(][a-zA-Z_0-9 ]{1,}[)]', part_2)
        prog_1 = re.compile('^pred[(]([a-zA-Z_0-9 ]{1,})[)]$')

        for part in gmm_search_parts:
            matching_parts.append(part)
            match_groups_multiple = prog_1.match(part)

            vars = match_groups_multiple.group(1).split()
            min_lag = 1
            max_lag = sys.maxsize

            self.process_GMM(vars, min_lag, max_lag, part)

    def parse_IV(self, matching_parts, part_2):
        
        iv_search_parts = re.findall('iv[(].{1,}[)]', part_2)
        prog_2 = re.compile('^iv[(](.{1,})[)]$')
        for part in iv_search_parts:
            matching_parts.append(part)
            match_groups_multiple = prog_2.match(part)
            vars = match_groups_multiple.group(1).split()
            invalid_name= self.parse_spaced_vars(vars, self.temp_iv_list)
            if invalid_name !='':
                print(part + ': ' + invalid_name + ' does not exist')
                exit()

    
    def parse_options(self, part_3):
        list_options = part_3.split()
        options = options_info()


        # possible_options=[{'onestep', 'iterated'},'nolevel', 'timedumm', 'collapse']
        for option in list_options:
            if option == 'onestep':
                options.steps = 1
            elif option == 'iterated':
                options.steps = 1000
            elif option == 'nolevel':
                options.level = False
            elif option == 'timedumm':
                options.timedumm = True
            elif option == 'collapse':
                options.collapse = True
            else:
                print(option + ' is not an option allowed')
                exit()

        return (options)


    def process_GMM(self,vars, min_lag, max_lag, part):
        if min_lag > max_lag:
            print(part + ': minimum lag cannot be greater than maximum lag')
            exit()
        if min_lag < 0:
            print(part + ': lags must be non-negative')
            exit()
        if len(vars) == 0:
            print(part + ': no variable is included')
            exit()
                
        for var in vars:
            if (var not in self.cols):
                print(part + ': ' + var + ' does not exist')
                exit()
            existing_names=[v.name for v in self.list_GMM]
            if var in existing_names:
                print(part + ': ' + var + ' cannot be declared in part 2 for twice or more')
                exit()   
            temp_var = gmm_var(var, min_lag, max_lag, 0)
            self.list_GMM.append(temp_var)


    def tbr_list(self, temp_list):
        tbr = []
        for i in range(len(temp_list.names)):
            var_name=temp_list.names[i]
            lags=temp_list.lags[i]
            num_lags=len(lags)
            for j in range(lags[0], lags[0]+num_lags):
                new_var=regular_variable(var_name, j)
                tbr.append(new_var)


        return (tbr)


    def check_dep_indep(self):
     
        self.temp_part1_list.purge()

        dep = self.temp_part1_list.names[0]
        dep_lags=self.temp_part1_list.lags[0]

        if dep_lags[0]!=0:
            print('dependent variable should not be lagged on the left hand side of the model')
            exit()

        if len(dep_lags)==0:
            print('lagged dependent variable should be included')
            exit()

        
        if len(dep_lags)>0 and dep_lags[1]!=1:
            print('lag 1 of the dependent variable is not included')
            exit()

        self.temp_part1_list.check_contiguous()


    def check_GMM(self):
        dep_name=self.temp_part1_list.names[0]

        for var in self.list_GMM:
            if var.name==dep_name:
                if var.min_lag<2:
                    print('must use lag 2 or earlier of the dependent variable as instruments')
                    exit()



    def check_iv(self):
        self.temp_iv_list.purge()
        self.temp_iv_list.check_contiguous()

    def check_three_lists(self):
        gmm_names=[var.name for var in self.list_GMM]
        iv_names=self.temp_iv_list.names

        for iv_name in iv_names:
            if iv_name in gmm_names:
                print('variable '+ iv_name + ': a variable can be either in GMM style or in IV style, but not both')
                exit()

        for i in range(len(self.temp_part1_list.lags)):
            var_name = self.temp_part1_list.names[i]
            var_lags=self.temp_part1_list.lags[i]
            bool_GMM=var_name in gmm_names
            bool_IV=var_name in iv_names

            if not (bool_GMM or bool_IV):
                self.temp_iv_list.insert(var_name, var_lags)

    

            


