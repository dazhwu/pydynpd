import pydynpd.variable
from pydynpd.variable import gmm_var, regular_variable
from pydynpd.info import  options_info
import sys
from sys import exit
import re

'''  n L(1\2).n L(0\1).w k| gmm(n, 2:4) gmm(w, 1:3) iv(ys L(0\1).k)| twostep nolevel robust'
 part_1  part_2 part_3
    
    gmm
    iv
    options
'''


def parse_command(command_str, df_col_names):
    global  cols
    cols=df_col_names

    variables = {}
    parts = command_str.split('|')
    if len(parts) <= 1:
        print('There should be at least two parts in command string')
        exit()
    else:
        part_1 = parts[0]
        dep_indep_list = parse_dep_indep(part_1)

        part_2 = parts[1]
        gmm_list, iv_list = parse_gmm_iv(part_2)

        # gmm_list=gmm_iv[0]
        # iv_list=gmm_iv[1]

        if len(parts) == 3:
            part_3 = parts[2]
            options=parse_options(part_3)
        else:
            options=options_info()

    variables['dep_indep'] = dep_indep_list
    variables['gmm'] = gmm_list
    variables['iv'] = iv_list

    return((variables, options))





def parse_spaced_vars(list_vars, indep_iv):
    # for independent variables (indep_iv=0) and iv variables (indep_iv=1)
    # list_vars is a list of variables
    tbr = []
    prog_1 = re.compile('^L\(([0-9]{1,})\/([0-9]{1,})\)[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
    prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')

    for var in list_vars:
        match_groups_multiple = prog_1.match(var)
        match_groups_single = prog_2.match(var)

        if match_groups_multiple:
            new_vars=gen_list_rhs(match_groups_multiple.group(3), int(match_groups_multiple.group(1)),
                                    int(match_groups_multiple.group(2)))
        elif prog_2.match(var):
            new_vars=gen_list_rhs(match_groups_single.group(2), int(match_groups_single.group(1)),
                                    int(match_groups_single.group(1)))
        else:
            new_vars=gen_list_rhs(var, 0,0)

        if len(new_vars) == 0:
            return ([], var)

        tbr+=new_vars

    return (tbr, '')

def parse_dep_indep(part_1):
    list_vars = part_1.split()



    list_dep_indep, missing_name = parse_spaced_vars(list_vars, 0)
    if len(list_dep_indep)==0:
        print(part_1 + ': ' + missing_name + ' does not exist' )
        exit()

    return (list_dep_indep)

def parse_gmm_iv(part_2):
    # list_vars=self.part_2.split()   #cannot use space to split

    list_gmm = []
    list_iv = []
    matching_parts=[]

    #gmm_search_parts = re.findall('gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{1,}[0-9]{1,}[)]', part_2)
    #prog_1 = re.compile('^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{1,}([0-9]{1,})[)]$')
    gmm_search_parts = re.findall('gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{1,}(?:(?:[.])|(?:[0-9]{1,}))[ ]{0,}[)]', part_2)
    prog_1 = re.compile('^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{1,}((?:[.])|(?:[0-9]{1,}))[ ]{0,}[)]$')
    for part in gmm_search_parts:
        matching_parts.append(part)
        # prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
        match_groups_multiple = prog_1.match(part)

        vars = match_groups_multiple.group(1).split()

        for var_name in vars:
            if  (var_name not in cols):
                print(part + ': ' + var_name + ' does not exist')
                exit()


        min_lag = int(match_groups_multiple.group(2))
        if match_groups_multiple.group(3)=='.':
            max_lag=sys.maxsize
        else:
            max_lag = int(match_groups_multiple.group(3))


        list_gmm=process_GMM(vars, min_lag,max_lag,list_gmm, part)

    # gmm_search_parts = re.findall('gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{1,}[.][)]', part_2)
    # prog_1 = re.compile('^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{1,}([.])[)]$')
    # for part in gmm_search_parts:
    #
    #     # prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
    #     match_groups_multiple = prog_1.match(part)
    #
    #     vars = match_groups_multiple.group(1).split()
    #     min_lag = int(match_groups_multiple.group(2))
    #     max_lag = sys.maxsize
    #
    #     list_gmm=process_GMM(vars, min_lag,max_lag,list_gmm, part)

    gmm_search_parts = re.findall('endo[(][a-zA-Z_0-9 ]{1,}[)]', part_2)
    prog_1 = re.compile('^endo[(]([a-zA-Z_0-9 ]{1,})[)]$')

    for part in gmm_search_parts:
        # prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
        matching_parts.append(part)
        match_groups_multiple = prog_1.match(part)

        vars = match_groups_multiple.group(1).split()
        min_lag = 2
        max_lag = sys.maxsize

        list_gmm=process_GMM(vars, min_lag,max_lag,list_gmm, part)

    gmm_search_parts = re.findall('pred[(][a-zA-Z_0-9 ]{1,}[)]', part_2)
    prog_1 = re.compile('^pred[(]([a-zA-Z_0-9 ]{1,})[)]$')

    for part in gmm_search_parts:
        # prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
        matching_parts.append(part)
        match_groups_multiple = prog_1.match(part)

        vars = match_groups_multiple.group(1).split()
        min_lag = 1
        max_lag = sys.maxsize

        list_gmm = process_GMM(vars, min_lag, max_lag, list_gmm, part)

    iv_search_parts = re.findall('iv[(][a-zA-Z_0-9 .(/)]{1,}[)]', part_2)
    prog_2 = re.compile('^iv[(]([a-zA-Z_0-9 .(/)]{1,})[)]$')
    for part in iv_search_parts:
        matching_parts.append(part)
        match_groups_multiple = prog_2.match(part)
        vars = match_groups_multiple.group(1).split()
        temp_list, strR= parse_spaced_vars(vars, 1)
        if len(temp_list)==0:
            print(part + ': ' + strR + ' does not exist')
            exit()
        else:
            list_iv = list_iv + temp_list

    part2_cpy=part_2
    for part in matching_parts:
        part2_cpy=part2_cpy.replace(part, '')

    if len(part2_cpy.strip())>0:
        print(part2_cpy.strip() + ': invalid GMM or IV statement')
        exit()


    return ([list_gmm, list_iv])

def parse_options(part_3):
    list_options = part_3.split()
    options=options_info()
    #possible_options=[{'onestep', 'iterated'},'nolevel', 'timedumm', 'collapse']
    for option in list_options:
        if option=='onestep':
            options.steps=1
        elif option=='iterated':
            options.steps=1000
        elif option=='nolevel':
            options.level=False
        elif option=='timedumm':
            options.timedumm=True
        elif option=='collapse':
            options.collapse=True
        else:
            print(option + ' not an option allowed')
            exit()

    return(options)

def process_GMM(vars, min_lag, max_lag, list_gmm, part):
    if min_lag > max_lag:
        print(part + ': minimum lag cannot be greater than maximum lag')
        exit()
    if min_lag < 0:
        print(part + ': lags must be non-negative')
        exit()
    if len(vars)==0:
        print(part + ': no variable is included before ,')
        exit()
    for var in vars:
        temp_var = gmm_var(var, min_lag, max_lag, 0)
        list_gmm.append(temp_var)

    return list_gmm

def gen_list_rhs(name, start_lag, end_lag):

    tbr=[]
    if name in cols:
        for i in range(start_lag, (end_lag+1)):
            new_var=regular_variable(name, i)
            tbr.append(new_var)

    return(tbr)



