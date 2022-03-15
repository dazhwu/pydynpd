from pandas import DataFrame

import numpy as np

from pydynpd.variable import regular_variable, gmm_var
from pydynpd.info import df_info, z_info, options_info
import time
from sys import exit

# https://stackoverflow.com/questions/29352511/numpy-sort-ndarray-on-multiple-columns

# https://itecnote.com/tecnote/python-efficiently-applying-a-function-to-a-grouped-pandas-dataframe-in-parallel/


def new_panel_data(df: DataFrame, identifiers, p_variables, options: options_info):
#    start=time.time()
    if len(identifiers) == 2:
        _individual = identifiers[0]
        _time = identifiers[1]
    else:
        raise Exception('two variables needed')

    variables = p_variables

    level = options.level
    timedumm=options.timedumm
    collapse=options.collapse
    method = 'fd'

    df_information = get_info(df, variables, method, _individual, _time)

    if timedumm:
        add_time_dummy(df, variables, _time, df_information.first_index, df_information.last_index)

    gmm_tables = get_gmm_table_list(df, variables, df_information, level)
    xy_tables = get_xy_table_list(df, variables, df_information)
    final_xy_tables = get_final_xy_tables(xy_tables, df_information, level)


    if level:
        z_information, z_list = build_z_level(variables, df_information, gmm_tables, collapse)
    else:
        z_information, z_list = build_z_diff(variables, df_information, gmm_tables, False, collapse)

#    print(time.time()-start)
    return ((z_list, z_information, df_information, final_xy_tables))


def get_info(df: DataFrame, variables, method, _individual, _time):
    max_lag = 0

    df['_individual'] = df[_individual].astype('category').cat.codes
    N = df['_individual'].unique().size

    df['_time'] = df[_time].astype('category').cat.codes
    T = df['_time'].unique().size

    df['_NT'] = df['_individual'] * T + df['_time']

    for var in variables['dep_indep'] + variables['iv']:
        if var.lag > max_lag:
            max_lag = var.lag

    for var in variables['gmm']:
        if var.min_lag > max_lag:
            max_lag = var.min_lag

    if method == 'fd':
        last_index = T - 1  # zero based
        first_index = max_lag + 1
    else:
        last_index = T - 2
        first_index = max_lag

    tbr = df_info(N=N, T=T, ids=['_NT'], _individual=_individual, _time=_time, max_lag=max_lag, first_index=first_index,
                  last_index=last_index)

    return (tbr)


def get_gmm_table_list(df: DataFrame, variables, df_information: df_info, level):
    iv_list = gen_ori_list(df, variables['iv'], df_information)
    gmm_list = gen_ori_list(df, variables['gmm'], df_information)
    Div_list = gen_fd_list(iv_list)
    gmm_tables = {}
    if level:  # sys-GMM
        Dgmm_list = gen_fd_list(gmm_list)
        gmm_tables['Dgmm'] = Dgmm_list

    gmm_tables['gmm'] = gmm_list
    gmm_tables['iv'] = iv_list

    gmm_tables['Div'] = Div_list

    return (gmm_tables)


def get_xy_table_list(df: DataFrame, variables: dict, df_information: df_info):
    xy_tables = {}

    num_var = len(variables['dep_indep'])
    y_list = gen_ori_list(df, variables['dep_indep'][0:1], df_information, True)
    x_list = gen_ori_list(df, variables['dep_indep'][1:num_var], df_information, True)

    Dy_list = gen_fd_list(y_list, True)
    Dx_list = gen_fd_list(x_list, True)

    xy_tables['x'] = x_list
    xy_tables['y'] = y_list
    xy_tables['Dx'] = Dx_list
    xy_tables['Dy'] = Dy_list

    return (xy_tables)


def get_final_xy_tables(xy_tables: dict, df_information: df_info, level):

    final_xy_tables = {}

    Dx_list = xy_tables['Dx']
    Dy_list = xy_tables['Dy']
    x_list = xy_tables['x']
    y_list = xy_tables['y']

    if level:  # sys-GMM
        Cx_list = []
        Cy_list = []

        height_upper = Dx_list[0].shape[0]
        height_lower = x_list[0].shape[0]
        height_total = height_upper + height_lower
        width=x_list[0].shape[1]

        for i in range(df_information.N):  # , nogil=True):  #df_information.N
            temp_y=np.empty((height_total,1), dtype='float64')
            temp_y[0:height_upper,0]=Dy_list[i][0:height_upper,0]
            temp_y[height_upper:height_total,0]=y_list[i][0:height_lower,0]
            # temp_y = np.vstack((Dy_list[i], y_list[i]))

            temp_x=np.empty((height_total, width+1))
            temp_x[0:height_upper, 0:width] = Dx_list[i][0:height_upper, 0:width]
            temp_x[height_upper:height_total, 0:width] = x_list[i][0:height_lower, 0:width]
            temp_x[0:height_upper, width]=0
            temp_x[height_upper:height_total, width] = 1
            #
            # temp_x = np.vstack((Dx_list[i], x_list[i]))
            # temp_constant = np.zeros((height_total, 1), dtype='float64')
            # temp_constant[height_upper:height_total, 0] = 1
            #
            # temp_x = np.hstack((temp_x, temp_constant))
            Cy_list.append(temp_y)
            Cx_list.append(temp_x)
    else:  # diff-GMM

        Cy_list = Dy_list
        Cx_list = Dx_list

    final_xy_tables['Cy'] = Cy_list
    final_xy_tables['Cx'] = Cx_list

    return (final_xy_tables)


def split_into_groups(arr, N, T):
    # needs to be sorted

    tbr=[]
    for i in range(0, N):
        temp_arr=np.empty((T, arr.shape[1]), dtype='float64')
        temp_arr[:]=arr[i*T:(i+1)*T,:]
        tbr.append(temp_arr)
    #tbr = np.vsplit(arr, N)
    return (tbr)


def build_z_level(variables: dict, info: df_info, gmm_tables: dict, collapse=False):

    lev_last_index = info.last_index
    lev_first_index = info.first_index - 1

    z_information, z_list = build_z_diff(variables, info, gmm_tables, True, collapse)

    width = lev_last_index - lev_first_index + 1

    gmm_vars = variables['gmm']
    iv_vars = variables['iv']
    Dgmm_list = gmm_tables['Dgmm']
    iv_list = gmm_tables['iv']

    height = len(gmm_vars) * width + len(iv_vars)

    start_row = z_information.diff_height  # z_list[0].shape[0]-height
    start_col = z_information.diff_width  # z_list[0].shape[1]-width

    for i in range(info.N):
        z = z_list[i]
        # z[start_row-1,start_col:(start_col+self.level_width)]=1
        z[z_information.diff_height + z_information.level_height - 1,
        start_col:(start_col + z_information.level_width)] = 1

        array_Dgmm = Dgmm_list[i]
        array_iv = iv_list[i]

        for var_id in range(len(gmm_vars)):
            lag = gmm_vars[var_id].min_lag - 1
            for j in range(width):
                if collapse:
                    z[start_row + var_id, start_col + j] = array_Dgmm[lev_first_index - lag + j, var_id]
                else:
                    z[start_row + var_id * width + j, start_col + j] = array_Dgmm[lev_first_index - lag + j, var_id]

        start_pos = z_information.num_gmm_instr
        for var_id in range(len(iv_vars)):
            var = iv_vars[var_id]
            z[start_pos + var_id, start_col:z_list[0].shape[1]] = array_iv[lev_first_index:(lev_last_index + 1), var_id]
    
        z[np.isnan(z)] = 0
        
    z_information.num_gmm_instr+=len(gmm_vars)
    z_information.num_instr+=z_information.level_height
    return ((z_information, z_list))


def build_z_diff(variables: dict, info: df_info, gmm_tables: dict, level, collapse=False):
    z_list = []

    gmm_vars = variables['gmm']
    iv_vars = variables['iv']

    gmm_list = gmm_tables['gmm']
    Div_list = gmm_tables['Div']

    diff_width = info.last_index - info.first_index + 1
    level_width = diff_width + 1
    if collapse:
        level_height = len(gmm_vars)  + 1  # + len(iv_vars)
    else:
        level_height = len(gmm_vars) * level_width  + 1 #+ len(iv_vars)

    num_gmm_instr, gmm_diff_info = prepare_Z_gmm_diff(variables, diff_width, info, collapse)
    iv_diff_info = prepare_Z_iv_diff(variables, diff_width, info)

    diff_height = (num_gmm_instr + iv_diff_info.shape[0])

    z_information = z_info(diff_height=diff_height, diff_width=diff_width, level_width=level_width,
                           level_height=level_height, num_gmm_instr=num_gmm_instr, num_instr=diff_height)
    
    #print(diff_height)

    if level:
        height = diff_height + level_height
        width = diff_width + level_width
    else:
        height = diff_height
        width = diff_width

    for i in range(info.N):
        z = np.zeros((height, width), dtype='float64')

        z_list.append(z)
        array_gmm = gmm_list[i]
        array_fd_iv = Div_list[i]

        var_id = 0
        for var in gmm_vars:

            for j in range(diff_width):
                row_pos = gmm_diff_info[var_id * 3 + 2, j]
                start = gmm_diff_info[var_id * 3 + 0, j]
                end = gmm_diff_info[var_id * 3 + 1, j]

                #z[row_pos:(row_pos + end - start + 1), j] = array_gmm[start:(end + 1), var_id]

                for k in range(end-start +1):
                    z[row_pos+k, j] = array_gmm[end-k, var_id]
            var_id += 1

        row_pos = num_gmm_instr

        var_id = 0
        for var_id in range(len(iv_vars)):
            var = iv_vars[var_id]
            for j in range(diff_width):
                index_to_take = iv_diff_info[var_id, j]

                z[row_pos, j] = array_fd_iv[index_to_take, var_id]

            row_pos += 1

        z[np.isnan(z)] = 0

    return ((z_information, z_list))


def prepare_Z_iv_diff(variables: dict, width, info: df_info):
    iv_vars = variables['iv']  # need to be placed at the beginning
    num_iv = len(iv_vars)

    t_info = np.empty((num_iv, width), dtype='int32')
    # num_iv_instr = 0
    var_id = 0
    for var_id in range(num_iv):
        var = iv_vars[var_id]
        t_info[var_id,] = range(info.first_index, info.last_index + 1)
        # num_iv_instr += width

    return (t_info)


def prepare_Z_gmm_diff(variables: dict, width, info: df_info, collapse=False):


    start_row = 0

    gmm_vars = variables['gmm']
    num_gmm = len(gmm_vars)
    t_info = np.empty((num_gmm * 3, width), dtype='int32')
    var_id = 0
    for var_id in range(num_gmm):
        var = gmm_vars[var_id]

        first_tend = info.first_index - var.min_lag
        last_tend = info.last_index - var.min_lag

        tend = np.arange(first_tend, last_tend + 1, dtype='int32')
        t_info[var_id * 3 + 1,] = tend


        for i in range(0, width):
            tstart = max(0, tend[i] + var.min_lag - var.max_lag)
            t_info[var_id * 3 + 0, i] = tstart
            # t_info[var_id*3+1, i]=tend[i]
            t_info[var_id * 3 + 2, i] = start_row   # physical position of the row

            num = tend[i] - tstart + 1
            # num_instru += num
            if collapse:
                if i==(width-1):
                    start_row += num
            else:
                start_row += num

    num_gmm_instr = start_row  #number of gmm instruments in diff eq
    return ((num_gmm_instr, t_info))


def gen_fd_list(ori_array_list, cut=False):
    tbr = []
    num_rows = ori_array_list[0].shape[0]
    num_cols = ori_array_list[0].shape[1]

    for ori_arr in ori_array_list:
        lag_arr = np.zeros((num_rows, num_cols), dtype='float64')
        tbr_arr = np.zeros((num_rows, num_cols), dtype='float64')
        lag_arr[range(0, 1), :] = np.NaN
        lag_arr[range(1, num_rows), :] = ori_arr[range(0, (num_rows - 1)), :]
        tbr_arr = ori_arr - lag_arr
        if cut:
            tbr_arr = tbr_arr[1:num_rows, :]
            # print('cut-----------------', start, end)
        tbr.append(tbr_arr)

    return (tbr)


def gen_ori_list(df: DataFrame, variable_list, info: df_info, cut=False):
    num_variables = len(variable_list)
    list_cols = info.ids.copy()
    # np.ndarray[np.double_t, ndim = 2] ori_arr  #pointer?
    tbr = []

    for var in variable_list:
        if var.name not in df.columns:
            print('Column ' + var.name + ' does not exist in the data set provided')
            exit()

        if var.name not in list_cols:
            list_cols.append(var.name)




    list_array = split_into_groups(make_balanced(df[list_cols].to_numpy(), info.N, info.T), info.N, info.T)

    variable_names = [var.name for var in variable_list]
    which_col = [list_cols.index(var_name) for var_name in variable_names]

    for i in range(info.N):
        ori_arr = np.empty((info.T, num_variables), dtype='float64')
        col = 0
        for j in range(num_variables):
            var = variable_list[j]
            if var.lag == 0:
                ori_arr[:, col] = list_array[i][:, which_col[j]]
            else:
                ori_arr[range(0, var.lag), col] = np.NaN
                ori_arr[range(var.lag, info.T), col] = list_array[i][range(0, info.T - var.lag), which_col[j]]

            col += 1
        if cut:
            ori_arr = ori_arr[(info.first_index - 1):(info.last_index + 1), :]
        tbr.append(ori_arr)

    return (tbr)


def make_balanced(ori, n_individual, n_time):
    # cdef:
    #     np.ndarray[np.float64_t, ndim = 2] arr_full

    arr_full = np.empty((n_individual * n_time, ori.shape[1]), dtype='float64')

    arr_full[:] = np.NaN
    # arr_full[:, 0] = np.repeat(range(0, N), T)
    arr_full[:, 0] = range(0, n_individual * n_time)
    # arr_full[:, 1] = arr_full[:, 2] % T

    mask = np.in1d(arr_full[:, 0], ori[:, 0])
    arr_full[mask, 1:arr_full.shape[1]] = ori[:, 1:ori.shape[1]]

    return (arr_full)

def add_time_dummy(df: DataFrame, variables: dict, _time: str, first_index, last_index):
    unique_time= sorted(df[_time].unique())[(first_index):(last_index+1)]

    prefix=_time + '_'
    for num in unique_time:
        name=prefix+str(num)
        df[name]=np.where(df[_time]==num, 1,0)
        new_var=regular_variable(name, 0)
        variables['dep_indep'].append(new_var)

        new_iv=regular_variable(name, 0)
        variables['iv'].append(new_var)










