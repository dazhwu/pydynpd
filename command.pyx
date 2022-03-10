import variable
from variable import gmm_var, regular_variable
import re

'''  n L(1\2).n L(0\1).w k| gmm(n, 2:4) gmm(w, 1:3) iv(ys L(0\1).k)| twostep nolevel robust'
 part_1  part_2 part_3
    
    gmm
    iv
    options
'''


def parse_command(command_str):
    cdef:
        list dep_indep_list, gmm_list, iv_list, parts
        str part_1, part_2, part_3
        dict variables={}
    variables = {}
    parts = command_str.split('|')
    if len(parts) <= 1:
        raise Exception('two parts at least')
    else:
        part_1 = parts[0]
        dep_indep_list = parse_dep_indep(part_1)

        part_2 = parts[1]
        gmm_list, iv_list = parse_gmm_iv(part_2)

        # gmm_list=gmm_iv[0]
        # iv_list=gmm_iv[1]

        if len(parts) == 3:
            part_3 = parts[2]

    variables['dep_indep'] = dep_indep_list
    variables['gmm'] = gmm_list
    variables['iv'] = iv_list

    return(variables)





def parse_spaced_vars(list_vars, indep_iv):
    # for independent variables (indep_iv=0) and iv variables (indep_iv=1)
    # list_vars is a list of variables
    cdef list tbr = []
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
            # https://rollbar.com/blog/throwing-exceptions-in-python/

            # raise Exception(var + ' is not a valid variable')

        tbr+=new_vars

    return (tbr)

def parse_dep_indep(str part_1):
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

def gen_list_rhs(str name, int start_lag, int end_lag):
    cdef:
        list tbr
        unsigned int i

    tbr=[]
    for i in range(start_lag, (end_lag+1)):
        new_var=regular_variable(name, i)
        tbr.append(new_var)

    return(tbr)



