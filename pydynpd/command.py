import pydynpd.variable
from pydynpd.variable import gmm_var, regular_variable
from pydynpd.info import  options_info
from sys import exit
import re

'''  n L(1\2).n L(0\1).w k| gmm(n, 2:4) gmm(w, 1:3) iv(ys L(0\1).k)| twostep nolevel robust'
 part_1  part_2 part_3
    
    gmm
    iv
    options
'''


def parse_command(command_str):

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


        tbr+=new_vars

    return (tbr)

def parse_dep_indep(part_1):
    list_vars = part_1.split()

    list_dep_indep = parse_spaced_vars(list_vars, 0)

    return (list_dep_indep)

def parse_gmm_iv(part_2):
    # list_vars=self.part_2.split()   #cannot use space to split

    list_gmm = []
    list_iv = []

    gmm_search_parts = re.findall('gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{1,}[0-9]{1,}[)]', part_2)
    prog_1 = re.compile('^gmm[(]([a-zA-Z_0-9 ]{1,})[,][ ]{0,}([0-9]{1,})[ ]{1,}([0-9]{1,})[)]$')
    for part in gmm_search_parts:

        # prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
        match_groups_multiple = prog_1.match(part)
        vars = match_groups_multiple.group(1).split()
        min_lag = int(match_groups_multiple.group(2))
        max_lag = int(match_groups_multiple.group(3))
        if min_lag > max_lag:
            print( part + ': minimum lag cannot be greater than maximum lag')
            exit()
        for var in vars:
            temp_var = gmm_var(var, min_lag, max_lag, 0)
            list_gmm.append(temp_var)



    iv_search_parts = re.findall('iv[(][a-zA-Z_0-9 .(/)]{1,}[)]', part_2)
    prog_2 = re.compile('^iv[(]([a-zA-Z_0-9 .(/)]{1,})[)]$')
    for part in iv_search_parts:
        match_groups_multiple = prog_2.match(part)
        vars = match_groups_multiple.group(1).split()
        list_iv = list_iv + parse_spaced_vars(vars, 1)

    return ([list_gmm, list_iv])

def parse_options(part_3):
    list_options = part_3.split()
    options=options_info()
    possible_options=['onestep','nolevel', 'timedumm', 'collapse']
    for option in list_options:
        if option=='onestep':
            options.twosteps=False
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

def gen_list_rhs(name, start_lag, end_lag):


    tbr=[]
    for i in range(start_lag, (end_lag+1)):
        new_var=regular_variable(name, i)
        tbr.append(new_var)

    return(tbr)



