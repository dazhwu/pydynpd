import numpy as np
from scipy import stats
import math
from pydynpd.info import hansen_test_info
from pydynpd.common_functions import sum_product


def hansen_overid(ZuuZ, zs, num_instru, num_indep):

    hansen_test=np.linalg.multi_dot([zs.transpose(),np.linalg.pinv(ZuuZ), zs])
    df=num_instru - num_indep
    crit=stats.chi2.ppf (q=0.95 , df=df)  #if hansen > crit then H0 is not supported
    p_value=1-stats.chi2.cdf(hansen_test, df)

    hansen_test=hansen_test_info(float(hansen_test), df, float(p_value), crit)

    
    return(hansen_test)

def AR_test(regression, m): #N, H, M, z_list, XZ_W,vcov, residual,residual_t, Cx_list, level, m):

    N=regression.N
    z_list=regression.z_list
    if regression.twosteps:
        M=regression.M2
        XZ_W=regression.XZ_W2
        vcov=regression.vcov_step2
        residual=regression.residual2
        residual_t=regression._residual2_t
    else:
        M=regression.M1
        XZ_W=regression.XZ_W1
        vcov=regression.vcov_step1
        residual = regression.residual1
        residual_t = regression._residual1_t

    Cx_list=regression.Cx_list
    if regression.level:
        r_list=[]
        r_t_list=[]
        x_list=[]
        diff_width=regression.z_information.diff_width
        for i in range(regression.N):
            temp=residual[i][0:diff_width,0:1]
            x=Cx_list[i][0:diff_width,:]
            r_list.append(temp)
            r_t_list.append(temp.transpose())
            x_list.append(x)
    else:
        r_list=residual
        r_t_list=residual_t
        x_list=Cx_list


    AR_list=[]
    for j in range(1, m+1):
        lagm_list=[]
        lagm_t_list=[]
        for i in range(N):
            residual_i =r_list[i]
            num_rows=residual_i.shape[0]

            lag_res=np.ndarray((num_rows, 1), dtype='float64')
            lag_res[range(0, j), 0] = 0 #np.NaN
            lag_res[range(j, num_rows), 0] = residual_i[range(0, (num_rows - j)), 0]
            # temp=residual.shift(j)
            lag_res[np.isnan(lag_res)] = 0

            lagm_list.append(lag_res.transpose())
            lagm_t_list.append(lag_res)
            
        d0=sum_product([lagm_list, r_list], N)

        d1=sum_product([ lagm_list, r_list, r_t_list, lagm_t_list], N)
        #d1 = sum_product([lagm_list, H, lagm_t_list], N)
        EX=sum_product([ lagm_list, x_list], N)

        temp3=sum_product([z_list, residual, r_t_list, lagm_t_list],N)
        #temp3 = sum_product([z_list, H, lagm_t_list], N)
        temp2=np.linalg.multi_dot([EX, M, XZ_W, temp3])
        d2=(-2)*temp2

        d3=np.linalg.multi_dot([EX,vcov,EX.transpose()])

        AR_temp=float(d0/math.sqrt(d1+d2+d3))


        AR_list.append(AR_temp)

    return (AR_list)
        
        
