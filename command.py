from dataclasses import dataclass
from variable import gmm_var, iv_var, dep_var, indep_var
import re


# 'n L(1\2).n L(0\1).w k| gmm(n, 2:4) gmm(w, 1:3) iv(ys L(0\1).k)| twostep nolevel robust'
class command:
    # part_1  part_2 part_3
    #
    # gmm
    # iv
    # options
    def __init__(self, command):  # constructor method
        parts = command.split('|')
        if len(parts) <= 1:
            raise Exception('two parts at least')
        else:
            self.part_1 = parts[0]
            dep_indep = self.parse_dep_indep()
            self.dep = dep_indep[0]
            self.indep = dep_indep[1]

            self.part_2 = parts[1]
            gmm_iv = self.parse_gmm_iv()
            self.gmm = gmm_iv[0]
            self.iv = gmm_iv[1]
            if len(parts) == 3:
                self.part_3 = parts[2]

    def parse_spaced_vars(self, list_vars, indep_iv):
        # for independent variables (indep_iv=0) and iv variables (indep_iv=1)
        # list_vars is a list of variables
        tbr = []
        prog_1 = re.compile('^L\(([0-9]{1,})\/([0-9]{1,})\)[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')
        prog_2 = re.compile('^L([0-9]{1,})[.]([a-zA-Z_]{1,}[a-zA-Z_0-9]{0,})$')

        for var in list_vars:
            match_groups_multiple = prog_1.match(var)
            match_groups_single = prog_2.match(var)

            if match_groups_multiple:
                if indep_iv == 0:
                    new_var = indep_var(match_groups_multiple.group(3), int(match_groups_multiple.group(1)),
                                        int(match_groups_multiple.group(2)))
                else:
                    new_var = iv_var(match_groups_multiple.group(3), int(match_groups_multiple.group(1)),
                                     int(match_groups_multiple.group(2)), 0)
            elif prog_2.match(var):
                if indep_iv == 0:
                    new_var = indep_var(match_groups_single.group(2), int(match_groups_single.group(1)),
                                        int(match_groups_single.group(1)))
                else:
                    new_var = iv_var(match_groups_single.group(2), int(match_groups_single.group(1)),
                                     int(match_groups_single.group(1)), 0)
            else:
                # https://rollbar.com/blog/throwing-exceptions-in-python/
                if indep_iv == 0:
                    new_var = indep_var(var, 0, 0)
                else:
                    new_var = iv_var(var, 0, 0, 0)
                # raise Exception(var + ' is not a valid variable')

            tbr.append(new_var)
        return (tbr)

    def parse_dep_indep(self):
        list_vars = self.part_1.split()
        j = 0
        list_indep = []
        dep = dep_var(list_vars[0])
        list_indep = self.parse_spaced_vars(list_vars[1:len(list_vars)], 0)

        # for var in list_vars:
        #
        #     if j==0:
        #         dep=dep_var(var)
        #         j += 1
        #     else:
        #         match_groups_multiple=prog_1.match(var)
        #         match_groups_single=prog_2.match(var)
        #         print(var)
        #
        #         if match_groups_multiple:
        #
        #             indep=indep_var(match_groups_multiple.group(3), int(match_groups_multiple.group(1)), int(match_groups_multiple.group(2)))
        #         elif prog_2.match(var):
        #             indep=indep_var(match_groups_single.group(2), int(match_groups_single.group(1)), int(match_groups_single.group(1)))
        #         else:
        #             #https://rollbar.com/blog/throwing-exceptions-in-python/
        #             indep=indep_var(var, 0, 0)
        #         #raise Exception(var + ' is not a valid variable')
        #         list_indep.append(indep)
        return ([dep, list_indep])

    def parse_gmm_iv(self):
        # list_vars=self.part_2.split()   #cannot use space to split

        list_gmm = []
        list_iv = []

        gmm_search_parts = re.findall('gmm[(][a-zA-Z_0-9 ]{1,}[,][ ]{0,}[0-9]{1,}[ ]{1,}[0-9]{1,}[)]', self.part_2)
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

        iv_search_parts = re.findall('iv[(][a-zA-Z_0-9 .(/)]{1,}[)]', self.part_2)
        prog_2 = re.compile('^iv[(]([a-zA-Z_0-9 .(/)]{1,})[)]$')
        for part in iv_search_parts:
            match_groups_multiple = prog_2.match(part)
            vars = match_groups_multiple.group(1).split()
            list_iv = list_iv + self.parse_spaced_vars(vars, 1)

        return ([list_gmm, list_iv])



