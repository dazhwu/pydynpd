import numpy as np

from pydynpd.info import df_info, z_info, options_info


class instruments(object):

    def __init__(self, variables: dict, gmm_tables: dict, df_information: df_info, options: options_info):
        level = options.level
        collapse = options.collapse
        self.z_information, self.gmm_diff_info, self.iv_diff_info, self.gmm_level_info = self.calculate_z_dimension(
            variables, df_information, level, collapse)

        if level:
            self.z_table = self.build_z_level(variables, gmm_tables, df_information, collapse)
        else:
            self.z_table = self.build_z_diff(variables, gmm_tables, df_information, False, collapse)

    def build_z_level(self, variables: dict, gmm_tables: dict, info: df_info, collapse=False):
        last_level_index = info.last_index
        first_level_index = info.first_level_index

        z_table = self.build_z_diff(
            variables, gmm_tables, info, True, collapse)

        level_width = self.z_information.level_width
        level_height = self.z_information.level_height
        diff_width = self.z_information.diff_width
        diff_height = self.z_information.diff_height
        width = self.z_information.width
        height = self.z_information.height

        Lgmm_vars = variables['Lgmm']
        iv_vars = variables['iv']
        Lgmm = gmm_tables['Lgmm']
        iv = gmm_tables['iv']
        Lgmm_dat = Lgmm.dat
        iv_dat = iv.dat
        Lgmm_iv_height = Lgmm.unit_height

        start_row = diff_height  # z_table[0].shape[0]-height
        start_col = diff_width  # z_table[0].shape[1]-width

        for i in range(info.N):
            z = z_table[(i * height):(i * height + height),:]
            # z[start_row-1,start_col:(start_col+self.level_width)]=1
            z[height - 1, start_col:width] = 1
            array_Lgmm = Lgmm_dat[(i * Lgmm_iv_height):(i * Lgmm_iv_height + Lgmm_iv_height), :]
            array_iv = iv_dat[(i * Lgmm_iv_height):(i * Lgmm_iv_height + Lgmm_iv_height), :]

            for var_id in range(len(Lgmm_vars)):
                lag = Lgmm_vars[var_id].min_lag
                t_info = self.gmm_level_info[3 * var_id: (3 * var_id + 3), :]
                for j in range(level_width):
                    the_index = t_info[0, j]
                    if the_index >= 0:
                        the_row = start_row + t_info[2, j]

                        the_index = t_info[0, j]

                        z[the_row, start_col + j] = array_Lgmm[the_index, var_id]

            start_pos = self.z_information.num_Dgmm_instr
            for var_id in range(len(iv_vars)):
                var = iv_vars[var_id]
                z[start_pos + var_id,
                start_col:width] = array_iv[first_level_index:(last_level_index + 1), var_id]

            z[np.isnan(z)] = 0

        return z_table

    def prepare_z_gmm_level(self, variables: dict, level_width, info: df_info, collapse=False):

        gmm_vars = variables['Lgmm']
        num_gmm = len(gmm_vars)
        t_info = np.empty((num_gmm * 3, level_width),
                          dtype='int32')  # row 0: gmm_index, 2: which row, 1: empty right now

        start_row = 0
        var_id = 0
        for var_id in range(num_gmm):
            var = gmm_vars[var_id]
            for i in range(level_width):
                the_index = info.first_level_index + i
                gmm_index = the_index - var.min_lag

                t_info[var_id * 3 + 2, i] = start_row

                if gmm_index >= 1:
                    t_info[var_id * 3 + 0, i] = gmm_index
                    if collapse:
                        if i == level_width - 1:
                            start_row += 1
                    else:
                        start_row += 1

                else:
                    t_info[var_id * 3 + 0, i] = -9999

        num_Lgmm_instr = start_row  # number of gmm instruments in diff eq

        return ((num_Lgmm_instr, t_info))

    def calculate_z_dimension(self, variables: dict, info: df_info, level, collapse=False):
        Lgmm_vars = variables['Lgmm']

        diff_width = info.last_index - info.first_diff_index + 1
        level_width = info.last_index - info.first_level_index + 1

        level_height = 0
        num_Lgmm_instr = 0
        gmm_level_info = None
        if level:
            num_Lgmm_instr, gmm_level_info = self.prepare_z_gmm_level(variables, level_width, info, collapse)
            level_height = num_Lgmm_instr + 1

        num_Dgmm_instr, gmm_diff_info = self.prepare_Z_gmm_diff(
            variables, diff_width, info, collapse)
        iv_diff_info = self.prepare_Z_iv_diff(variables, diff_width, info)

        diff_height = (num_Dgmm_instr + iv_diff_info.shape[0])

        if level:
            height = diff_height + level_height
            width = diff_width + level_width
        else:
            height = diff_height
            width = diff_width

        z_information = z_info(diff_height=diff_height, diff_width=diff_width, level_width=level_width,
                               level_height=level_height, height=height, width=width,
                               num_Dgmm_instr=num_Dgmm_instr, num_Lgmm_instr=num_Lgmm_instr, num_instr=height)

        return (z_information, gmm_diff_info, iv_diff_info, gmm_level_info)

    def build_z_diff(self, variables: dict, gmm_tables: dict, info: df_info, level, collapse=False):

        gmm_vars = variables['Dgmm']
        iv_vars = variables['iv']
        Dgmm = gmm_tables['Dgmm']
        Div = gmm_tables['Div']
        Dgmm_dat = Dgmm.dat
        Div_dat = Div.dat
        gmm_Div_height = Dgmm.unit_height

        diff_width = self.z_information.diff_width
        height = self.z_information.height
        width = self.z_information.width
        num_Dgmm_instr = self.z_information.num_Dgmm_instr

        z_table = np.zeros((height * info.N, width), dtype=np.float64)

        for i in range(info.N):
            z = z_table[i * height:(i + 1) * height, :]

            array_gmm = Dgmm_dat[i * gmm_Div_height:(i + 1) * gmm_Div_height, :]
            array_fd_iv = Div_dat[i * gmm_Div_height:(i + 1) * gmm_Div_height, :]

            var_id = 0
            for var in gmm_vars:
                for j in range(diff_width):
                    row_pos = self.gmm_diff_info[var_id * 3 + 2, j]
                    start = self.gmm_diff_info[var_id * 3 + 0, j]
                    end = self.gmm_diff_info[var_id * 3 + 1, j]

                    for k in range(end - start + 1):
                        z[row_pos + k, j] = array_gmm[end - k, var_id]
                var_id += 1

            row_pos = num_Dgmm_instr

            var_id = 0
            for var_id in range(len(iv_vars)):
                var = iv_vars[var_id]
                for j in range(diff_width):
                    index_to_take = self.iv_diff_info[var_id, j]

                    z[row_pos, j] = array_fd_iv[index_to_take, var_id]

                row_pos += 1

            z[np.isnan(z)] = 0

        return z_table

    def prepare_Z_iv_diff(self, variables: dict, width, info: df_info):
        iv_vars = variables['iv']  # need to be placed at the beginning
        num_iv = len(iv_vars)

        t_info = np.empty((num_iv, width), dtype='int32')
        # num_iv_instr = 0
        var_id = 0
        for var_id in range(num_iv):
            var = iv_vars[var_id]
            t_info[var_id,] = range(info.first_diff_index, info.last_index + 1)
            # num_iv_instr += width

        return (t_info)

    def prepare_Z_gmm_diff(self, variables: dict, width, info: df_info, collapse=False):
        start_row = 0

        gmm_vars = variables['Dgmm']
        num_gmm = len(gmm_vars)
        t_info = np.empty((num_gmm * 3, width), dtype='int32')
        var_id = 0
        for var_id in range(num_gmm):
            var = gmm_vars[var_id]

            first_tend = info.first_diff_index - var.min_lag
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
