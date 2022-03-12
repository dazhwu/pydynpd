import numpy as np
from scipy import stats
import math
from info import hansen_test_info


def hansen_overid(N, W, zs, num_instru, num_indept):

    hansen_test=np.linalg.multi_dot([zs.transpose(),np.linalg.pinv(W*N), zs])
    df=num_instru - num_indept
    crit=stats.chi2.ppf (q=0.95 , df=df)  #if hansen > crit then H0 is not supported
    p_value=1-stats.chi2.cdf(hansen_test, df)

    hansen_test=hansen_test_info(hansen_test, df, float(p_value), crit)

    
    return(hansen_test)

def AR_test(N, H, M, z_list, XZ_W,vcov, residual,residual_t, Cx_list, m):
    #lag_residual list of list
    # array of list
    
    #区分 H1 H2 vcov_stepAR_list.append(AR_temp)
    AR_list=[]
    for j in range(1, m+1):
        lagm_list=[]
        lagm_t_list=[]
        for i in range(N):
            residual_i =residual[i]
            num_rows=residual_i.shape[0]

            lag_res=np.ndarray((num_rows, 1), dtype='float64')
            lag_res[range(0, j), 0] = 0 #np.NaN
            lag_res[range(j, num_rows), 0] = residual_i[range(0, (num_rows - j)), 0]
            # temp=residual.shift(j)
            lag_res[np.isnan(lag_res)] = 0
            # if i==0:
            #     print(residual_i)
            #     print('----')
                
            #     print(lag_res)
            # replace nan by 0
            lagm_list.append(lag_res.transpose())
            lagm_t_list.append(lag_res)
            
        d0=sum_product([lagm_list, residual], N)

        d1=sum_product([ lagm_list, residual, residual_t, lagm_t_list], N)

        EX=sum_product([ lagm_list, Cx_list], N)

        temp3=sum_product([z_list, residual, residual_t, lagm_t_list],N)

        temp2=np.linalg.multi_dot([EX, M, XZ_W, temp3])
        d2=-2*temp2

        d3=np.linalg.multi_dot([EX,vcov,EX.transpose()])

        AR_temp=float(d0/math.sqrt(d1+d2+d3))

        
        AR_list.append(AR_temp)

    return (AR_list)
        
        
def sum_product(listOflist, n_rows):
    num_elements=len(listOflist)

    for i in range(n_rows):
        list_temp=[]
        for j in range(num_elements):
            if type(listOflist[j])==list:
                var_list=listOflist[j]
                list_temp.append(var_list[i])
            elif type(listOflist[j])==np.ndarray:
                var_mat=listOflist[j]
                list_temp.append(var_mat)
            else:
                pass  #throw error
        temp=np.linalg.multi_dot(list_temp)
        if i==0:
            tbr=temp
        else:
            tbr=np.add(tbr, temp)

    return(tbr)            
    
            
        
    