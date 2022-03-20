import numpy as np
from scipy import stats
import math
from pydynpd.info import hansen_test_info



def hansen_overid(W2_inv, N, zs,    num_instru, num_indep):

    hansen_test=np.linalg.multi_dot([zs.transpose(),W2_inv, zs])*(1.0/N)
    df=num_instru - num_indep
    crit=stats.chi2.ppf (q=0.95 , df=df)  #if hansen > crit then H0 is not supported
    p_value=1-stats.chi2.cdf(hansen_test, df)

    hansen_test=hansen_test_info(float(hansen_test), df, float(p_value), crit)


    return(hansen_test)

def AR_test(regression, step, m): #N, H, M, z_list, XZ_W,vcov, residual,residual_t, Cx_list, level, m):

    N=regression.N
    z_list=regression.z_list
    Cx_list = regression.Cx_list

    step1 = regression.result_list[0]
    step2=regression.result_list[1]
    if step==2:
        current_step=step2
    elif step==1:
        current_step=step1
    else:
        current_step=regression.result_list[step-1]

    M=current_step.M
    XZ_W=current_step._XZ_W
    M_XZ_W=current_step._M_XZ_W  #??????????????????????????????
    vcov=current_step.vcov
    residual = current_step.residual
    residual_t = current_step._residual_t

    r_height=int(residual.shape[0]/regression.N)


    x_height=int(Cx_list.shape[0]/regression.N)
    z_height=int(z_list.shape[0]/regression.N)
    x_width=Cx_list.shape[1]

    diff_width = regression.z_information.diff_width

    if regression.level:

        # r_list=[]
        # r_t_list=[]
        # x_list=[]
        r_list=np.empty((diff_width*regression.N, 1), dtype=np.float64)
        r_t_list=np.empty((1, diff_width*regression.N), dtype=np.float64)
        x_list=np.empty((diff_width*regression.N, x_width), dtype=np.float64)

        for i in range(regression.N):
            r=residual[(i*r_height):(i*r_height+r_height),:]
            r_list[(i*diff_width):(i*diff_width+diff_width),0:1]=r[0:diff_width,0:1]
            #
            x_table=Cx_list[(i*x_height):(i*x_height+x_height),:]
            x_list[(i*diff_width):(i*diff_width+diff_width),:]=x_table[0:diff_width,:]
            #
            r_t=residual_t[:, (i*r_height):(i*r_height+r_height)]            #
            r_t_list[0:1, (i*diff_width):(i*diff_width+diff_width)]=r_t[0:1, 0:diff_width]

        r_height=diff_width
        x_height=diff_width

    else:
        r_list=residual
        r_t_list=residual_t
        x_list=Cx_list


    AR_list=[]
    for j in range(1, m+1):
        lagm_list=np.ndarray((1, r_height*N), dtype='float64')
        lagm_t_list=np.ndarray((r_height*N, 1), dtype='float64')

        r0_height = int(residual.shape[0] / regression.N)

        for i in range(N):
        # calculate lag_res and lag_res_t
            r_i =r_list[(r_height*i):(r_height*i+r_height),0:1]
            lag_res_t=lagm_t_list[(i*r_height):(i*r_height+r_height),0:1]
            lag_res_t[range(0, j), 0:1] = 0 #np.NaN
            lag_res_t[range(j, r_height), 0:1] = r_i[range(0, (r_height - j)), 0:1]
            # temp=residual.shift(j)
            lag_res_t[np.isnan(lag_res_t)] = 0

            r_t_i=r_t_list[0:1, (r_height*i):(r_height*i+r_height)]
            lag_res=lagm_list[0:1, (i*r_height):(i*r_height+r_height)]
            lag_res[0, range(0, j)] = 0  # np.NaN
            lag_res[0, range(j, r_height)] = r_t_i[0,range(0, (r_height - j))]
        # temp=residual.shift(j)
            lag_res[np.isnan(lag_res)] = 0

            x = x_list[(x_height * i):(x_height * i + x_height), :]
            z=z_list[(z_height * i):(z_height * i + z_height), :]
            r_whole_i=residual[(i*r0_height):(i*r0_height+r0_height),0:1]
            r_whole_t_i=residual_t[0:1,(i*r0_height):(i*r0_height+r0_height)]

            d0_temp=lag_res @ r_i
            d1_temp=d0_temp @ r_t_i @ lag_res_t
            EX_temp=lag_res @ x

            r_t = r_t_list[0:1, (i*diff_width):(i*diff_width+diff_width)]
            temp3_temp=z @ r_whole_i @ r_t @ lag_res_t

            if i==0:
                d0=d0_temp
                d1=d1_temp
                EX=EX_temp
                temp3=temp3_temp
            else:
                d0 += d0_temp
                d1 += d1_temp
                EX += EX_temp
                temp3 += temp3_temp

        #temp3 = sum_product([z_list, H, lagm_t_list], N)
        temp2=np.linalg.multi_dot([EX, M_XZ_W, temp3])
        d2=(-2)*temp2

        d3=np.linalg.multi_dot([EX,vcov,EX.transpose()])

        AR_temp=float(d0/math.sqrt(d1+d2+d3))


        AR_list.append(AR_temp)

    return (AR_list)
        
        
