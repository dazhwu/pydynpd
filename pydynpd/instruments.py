import numpy as np

from pydynpd.info import df_info, z_info, options_info


class instruments(object):

    def __init__(self, variables: dict, gmm_tables: dict, df_information: df_info, options: options_info):
        level = options.level
        collapse = options.collapse
        self.z_information, self.gmm_diff_info, self.iv_diff_info \
            = self.calculate_z_dimension(variables, df_information, level, collapse)

        if level:
            self.z_list = self.build_z_level(variables, gmm_tables, df_information,
                                             self.z_information, self.gmm_diff_info, self.iv_diff_info, collapse)
        else:
            self.z_list = self.build_z_diff(variables, gmm_tables, df_information,
                                            self.z_information, self.gmm_diff_info, self.iv_diff_info, False, collapse)

    def build_z_level(self, variables: dict, gmm_tables: dict, info: df_info, z_information: z_info, gmm_diff_info,
                      iv_diff_info,
                      collapse=False):
        lev_last_index = info.last_index
        lev_first_index = info.first_index - 1

        z_list = self.build_z_diff(
            variables, gmm_tables, info, z_information, gmm_diff_info, iv_diff_info, True, collapse)

        level_width = z_information.level_width
        level_height = z_information.level_height
        diff_width = z_information.diff_width
        diff_height = z_information.diff_height
        width = z_information.width
        height = z_information.height

        gmm_vars = variables['gmm']
        iv_vars = variables['iv']
        Dgmm_list = gmm_tables['Dgmm']
        iv_list = gmm_tables['iv']
        Dgmm_iv_height=int(iv_list.shape[0]/info.N)
        # height = len(gmm_vars) * width + len(iv_vars)

        start_row = diff_height  # z_list[0].shape[0]-height
        start_col = diff_width  # z_list[0].shape[1]-width

        for i in range(info.N):
            z = z_list[(i * height):(i * height + height)]
            # z[start_row-1,start_col:(start_col+self.level_width)]=1
            z[height - 1, start_col:width] = 1

            array_Dgmm = Dgmm_list[(i*Dgmm_iv_height):(i*Dgmm_iv_height+Dgmm_iv_height),:]
            array_iv = iv_list[(i*Dgmm_iv_height):(i*Dgmm_iv_height+Dgmm_iv_height),:]

            for var_id in range(len(gmm_vars)):
                lag = gmm_vars[var_id].min_lag - 1
                for j in range(level_width):
                    if collapse:
                        z[start_row + var_id, start_col +
                          j] = array_Dgmm[lev_first_index - lag + j, var_id]
                    else:
                        z[start_row + var_id * level_width + j, start_col +
                          j] = array_Dgmm[lev_first_index - lag + j, var_id]

            start_pos = z_information.num_gmm_instr
            for var_id in range(len(iv_vars)):
                var = iv_vars[var_id]
                z[start_pos + var_id,
                start_col:width] = array_iv[lev_first_index:(lev_last_index + 1), var_id]

            z[np.isnan(z)] = 0

        z_information.num_gmm_instr += len(gmm_vars)
        z_information.num_instr += z_information.level_height
        return z_list

    def calculate_z_dimension(self, variables: dict, info: df_info, level, collapse=False):
        gmm_vars = variables['gmm']

        diff_width = info.last_index - info.first_index + 1
        level_width = diff_width + 1
        if collapse:
            level_height = len(gmm_vars) + 1  # + len(iv_vars)
        else:
            level_height = len(gmm_vars) * level_width + 1  # + len(iv_vars)

        num_gmm_instr, gmm_diff_info = self.prepare_Z_gmm_diff(
            variables, diff_width, info, collapse)
        iv_diff_info = self.prepare_Z_iv_diff(variables, diff_width, info)

        diff_height = (num_gmm_instr + iv_diff_info.shape[0])

        if level:
            height = diff_height + level_height
            width = diff_width + level_width
        else:
            height = diff_height
            width = diff_width

        z_information = z_info(diff_height=diff_height, diff_width=diff_width, level_width=level_width,
                               level_height=level_height, height=height, width=width,
                               num_gmm_instr=num_gmm_instr, num_instr=diff_height)

        return (z_information, gmm_diff_info, iv_diff_info)

    def build_z_diff(self, variables: dict, gmm_tables: dict, info: df_info, z_information: z_info, gmm_diff_info,
                     iv_diff_info,
                     level, collapse=False):

        gmm_vars = variables['gmm']
        iv_vars = variables['iv']

        gmm_list = gmm_tables['gmm']
        Div_list = gmm_tables['Div']
        gmm_Div_height=int(gmm_list.shape[0]/info.N)

        diff_width = z_information.diff_width
        height = z_information.height
        width = z_information.width
        num_gmm_instr = z_information.num_gmm_instr

        z_list = np.zeros((height * info.N, width), dtype=np.float64)

        for i in range(info.N):
            z = z_list[i * height:(i + 1) * height, :]

            array_gmm = gmm_list[i * gmm_Div_height:(i + 1) * gmm_Div_height, :]
            array_fd_iv = Div_list[i * gmm_Div_height:(i + 1) * gmm_Div_height, :]

            var_id = 0
            for var in gmm_vars:

                for j in range(diff_width):
                    row_pos = gmm_diff_info[var_id * 3 + 2, j]
                    start = gmm_diff_info[var_id * 3 + 0, j]
                    end = gmm_diff_info[var_id * 3 + 1, j]

                    # z[row_pos:(row_pos + end - start + 1), j] = array_gmm[start:(end + 1), var_id]

                    for k in range(end - start + 1):
                        z[row_pos + k, j] = array_gmm[end - k, var_id]
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

        return z_list

    def prepare_Z_iv_diff(self, variables: dict, width, info: df_info):
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

    def prepare_Z_gmm_diff(self, variables: dict, width, info: df_info, collapse=False):
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

            for i in range(width):
                tstart = max(0, tend[i] + var.min_lag - var.max_lag)
                t_info[var_id * 3 + 0, i] = tstart
                # t_info[var_id*3+1, i]=tend[i]
                # physical position of the row
                t_info[var_id * 3 + 2, i] = start_row

                num = tend[i] - tstart + 1
                # num_instru += num
                if collapse:
                    if i == (width - 1):
                        start_row += num
                else:
                    start_row += num

        num_gmm_instr = start_row  # number of gmm instruments in diff eq
        return ((num_gmm_instr, t_info))
