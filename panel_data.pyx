from pandas import DataFrame
import numpy as np
cimport numpy


# https://stackoverflow.com/questions/29352511/numpy-sort-ndarray-on-multiple-columns

# https://itecnote.com/tecnote/python-efficiently-applying-a-function-to-a-grouped-pandas-dataframe-in-parallel/

class panel_data():

    def __init__(self, df: DataFrame, identifiers, variables, level):
        # self.ids=['_individual', '_time', '_NT']
        self.first_index = None
        self.last_index = None
        self.method = None
        self.max_lag = None
        self.level=level
        self.ids = ['_NT']

        if len(identifiers) == 2:
            self.individual = identifiers[0]
            self.time = identifiers[1]
        else:
            raise Exception('two variables needed')

        df['_individual'] = df[self.individual].astype('category').cat.codes
        self.N = df['_individual'].unique().size

        df['_time'] = df[self.time].astype('category').cat.codes
        self.T = df['_time'].unique().size
        df['_NT'] = df['_individual'] * self.T + df['_time']
        self.variables = variables
        self.get_lag_info()

        num_var=len(self.variables['dep_indep'])
        
        self.iv_list = self.gen_ori_list(df, self.variables['iv'])
        self.gmm_list = self.gen_ori_list(df, self.variables['gmm'])
        self.Div_list = self.gen_fd_list(self.iv_list)
        
        
        y_list = self.gen_ori_list(df, self.variables['dep_indep'][0:1],True)
        x_list=self.gen_ori_list(df, self.variables['dep_indep'][1:num_var], True)
    
        Dy_list = self.gen_fd_list(y_list,  True)
        Dx_list = self.gen_fd_list(x_list,  True)
        
        #if diff then do not calculate Dgmm
        #if level then 一次性生成 cy cx
        
        if self.level:  #sys-GMM
            self.Dgmm_list=self.gen_fd_list(self.gmm_list)
            
            self.Cx_list=[]
            self.Cy_list=[]
            height_upper=Dx_list[0].shape[0]
            height_lower=x_list[0].shape[0]
            height_total=height_upper+height_lower
            
            for i in range(self.N):  
                
                temp_y=np.vstack((Dy_list[i], y_list[i]))
                temp_x=np.vstack((Dx_list[i], x_list[i]))
                temp_constant=np.zeros((height_total,1), dtype='float64')
                temp_constant[height_upper:height_total,0]=1
   
                temp_x=np.hstack((temp_x, temp_constant))
                self.Cy_list.append(temp_y)
                self.Cx_list.append(temp_x)
        else:  #diff-GMM        
            
            self.Cy_list = Dy_list
            self.Cx_list = Dx_list

        
        
        self.build_z_diff()
        self.build_z_level()


        # self.df = self.make_balanced(df, self.N, self.T)

        # turn df to numpy array

    def get_lag_info(self):
        cdef int max_lag
        max_lag = 0

        for var in self.variables['dep_indep'] + self.variables['iv']:
            if var.lag > max_lag:
                max_lag = var.lag
        for var in self.variables['gmm']:
            if var.min_lag > max_lag:
                max_lag = var.min_lag

        self.max_lag = max_lag

        self.method = 'fd'
        if self.method == 'fd':
            self.last_index = self.T - 1  # zero based
            self.first_index = self.max_lag + 1
        else:
            self.last_index = self.T - 2
            self.first_index = self.max_lag


    def check_na_diff(self):
        cdef int i
        tbr = []
        for i in range(0, self.N):
            temp = self.Dxy_list[i]
            row_if_nan = np.isnan(temp).any(axis=1)
            tbr.append(row_if_nan)
        return (tbr)

    def split_into_groups(self, arr):
        # needs to be sorted

        # tbr=[]
        # for i in range(0, self.N):
        #     temp_arr=np.empty((self.T, arr.shape[1]), dtype='float64')
        #     temp_arr[:]=arr[i*self.T:(i+1)*self.T,:]
        #     tbr.append(temp_arr)
        tbr = np.vsplit(arr, self.N)
        return (tbr)

    def build_z_level(self):
        cdef int i
        lev_last_index=self.last_index
        lev_first_index=self.first_index-1

        width=lev_last_index -lev_first_index+1

        gmm_list = self.variables['gmm']
        iv_list = self.variables['iv']

        height=len(gmm_list)*width + len(iv_list)

        start_row=self.diff_height #self.z_list[0].shape[0]-height
        start_col=self.diff_width #self.z_list[0].shape[1]-width

        for i in range(self.N):
            z=self.z_list[i]

            #z[start_row-1,start_col:(start_col+self.level_width)]=1
            z[self.diff_height+self.level_height-1, start_col:(start_col+self.level_width)]=1

            array_Dgmm = self.Dgmm_list[i]

            array_iv = self.iv_list[i]


            for var_id in range(len(gmm_list)):
                lag=gmm_list[var_id].min_lag-1
                for j in range(width):
                    #print(start_row+var_id*width+j, start_col+j, lev_first_index-1+j)
                    z[start_row+var_id*width+j, start_col+j]=array_Dgmm[ lev_first_index-lag+j, var_id]


            var_id=0
            #start_pos=start_row+len(gmm_list)*width
            start_pos = self.num_gmm_instr
            for var in iv_list:
                z[start_pos + var_id, start_col:self.z_list[0].shape[1]] = array_iv[ (lev_first_index):(lev_last_index+1), var_id]
                var_id+=1


    def build_z_diff(self):
        cdef int diff_width, level_width, level_height, height, i, var_id

        diff_width = self.last_index - self.first_index + 1
        level_width=diff_width+1

        gmm_list = self.variables['gmm']
        iv_list = self.variables['iv']



        level_height=len(gmm_list)*level_width + len(iv_list)

        gmm_diff_info = self.prepare_Z_gmm_diff(diff_width)
        iv_diff_info = self.prepare_Z_iv_diff(diff_width)
        #na_list = self.check_na_diff()

        height = (self.num_gmm_instr + iv_diff_info.shape[0])


        self.diff_width=diff_width
        self.level_width=level_width
        self.diff_height=height
        self.level_height=level_height
        if self.level:
            height=height+level_height
            width=diff_width+level_width
        else:
            width=diff_width


        self.z_list = []
        for i in range(0, self.N):
            z = np.zeros((height, width), dtype='float64')

            self.z_list.append(z)
            array_gmm = self.gmm_list[i]

            array_fd_iv = self.Div_list[i]


            var_id = 0
            for var in gmm_list:
                # which_variable = self.reg_df.columns.values.tolist().index(var.name)

                for j in range(0, diff_width):
                    row_pos = gmm_diff_info[var_id * 3 + 2, j]
                    start = gmm_diff_info[var_id * 3 + 0, j]
                    end = gmm_diff_info[var_id * 3 + 1, j]
                    z[row_pos:(row_pos + end - start + 1), j] = array_gmm[start:(end + 1), var_id]

                var_id += 1

            row_pos = self.num_gmm_instr

            var_id = 0
            for var in iv_list:
                # which_variable = self.fd_iv.columns.values.tolist().index(var.name)

                for j in range(0, diff_width):
                    index_to_take = iv_diff_info[var_id, j]
                    # print (index_to_take, i, which_variable, var_id)

                    z[row_pos, j] = array_fd_iv[index_to_take, var_id]
                    # print(z[row_pos,j], '--', i*self.T+index_to_take, var_id)
                row_pos += 1

                var_id += 1
            z[np.isnan(z)] = 0



    def prepare_Z_iv_diff(self, int width):

        cdef int var_id

        iv_list = self.variables['iv']
        num_iv = len(iv_list)

        t_info = np.empty((num_iv, width), dtype='int32')
        self.num_iv_instr = 0
        var_id = 0
        for var in iv_list:
            t_info[var_id,] = range(self.first_index, self.last_index + 1)
            self.num_iv_instr += width
            var_id += 1

        return (t_info)

    def prepare_Z_gmm_diff(self, width):

        cdef int start_index, var_id, i,  num

        start_index = 0
        gmm_list = self.variables['gmm']
        num_gmm = len(gmm_list)
        t_info = np.empty((num_gmm * 3, width), dtype='int32')
        var_id = 0
        for var in gmm_list:
            # print(var.min_lag, var.max_lag)
            first_tend = self.first_index - var.min_lag
            last_tend = self.last_index - var.min_lag
            # print(first_tend, last_tend)
            # first_ti=max(0, self.first_index-var.max_lag)

            tend = np.arange(first_tend, last_tend + 1)
            t_info[var_id * 3 + 1,] = tend
            for i in range(0, width):
                tstart = max(0, tend[i] + var.min_lag - var.max_lag)
                t_info[var_id * 3 + 0, i] = tstart
                # t_info[var_id*3+1, i]=tend[i]
                t_info[var_id * 3 + 2, i] = start_index

                num = tend[i] - tstart + 1
                # num_instru += num
                start_index += num

            var_id += 1

        self.num_gmm_instr = start_index
        return (t_info)

    def gen_fd_list(self, ori_array_list,   cut=False):

        cdef int num_rows, num_cols
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
                #print('cut-----------------', start, end)
            tbr.append(tbr_arr)

        return (tbr)

    def gen_ori_list(self, df: DataFrame, variable_list,  cut=False):
        cdef int num_variables, i, j, col

        list_cols = self.ids
        num_variables = len(variable_list)



        for var in variable_list:
            if var.name not in list_cols:
                list_cols.append(var.name)

        temp_array = self.split_into_groups(self.make_balanced(df[list_cols].to_numpy(), self.N, self.T))

        variable_names = [var.name for var in variable_list]
        which_col = [list_cols.index(var_name) for var_name in variable_names]

        tbr = []

        for i in range(0, self.N):
            ori_arr = np.empty((self.T, num_variables), dtype='float64')
            col = 0
            for j in range(0, num_variables):
                var = variable_list[j]
                if var.lag == 0:
                    ori_arr[:, col] = temp_array[i][:, which_col[j]]
                else:
                    ori_arr[range(0, var.lag), col] = np.NaN
                    ori_arr[range(var.lag, self.T), col] = temp_array[i][range(0, self.T - var.lag), which_col[j]]

                col += 1
            if cut:
                ori_arr=ori_arr[(self.first_index-1):(self.last_index+1),:]
            tbr.append(ori_arr)

        return (tbr)

    # def select_columns(self, df):

    #     list_cols = self.ids
    #     for var in self.variables['gmm']:
    #         if var.name not in list_cols:
    #             list_cols.append(var.name)
    #     self.gmm = self.split_into_groups(self.make_balanced(df[list_cols].to_numpy(), self.N, self.T))
    #     self.gmm_names = list_cols

    #     list_cols = self.ids
    #     for var in self.variables['iv']:
    #         if var.name not in list_cols:
    #             list_cols.append(var.name)
    #     self.iv = self.split_into_groups(self.make_balanced(df[list_cols].to_numpy(), self.N, self.T))
    #     self.iv_names = list_cols

    def make_balanced(self, ori, int N, int T):

        arr_full = np.empty((N * T, ori.shape[1]), dtype='float64')
        arr_full[:] = np.NaN
        # arr_full[:, 0] = np.repeat(range(0, N), T)
        arr_full[:, 0] = range(0, N * T)
        # arr_full[:, 1] = arr_full[:, 2] % T

        mask = np.in1d(arr_full[:, 0], ori[:, 0])
        arr_full[mask, 1:arr_full.shape[1]] = ori[:, 1:ori.shape[1]]
        return (arr_full)
