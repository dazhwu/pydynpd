import numpy as np
import pandas as pd
import statsmodels.api as sm
#from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse

my_data = np.array([[5, 1, 1],
                    [3, 2, 3],
                    [1, 3, 2],
                    [3, 1, 1],
                    [4, 2, 2],
                    [7, 3, 1],
                    [7, 1, 1]])


df = pd.DataFrame(data=my_data, columns=['y', 'dummy', 'x'])
just_dummies = pd.get_dummies(df['dummy'])

step_1 = pd.concat([df, just_dummies], axis=1)
step_1.drop(['dummy'], inplace=True, axis=1)
A=step_1.to_numpy()
B=np.cov(A, rowvar=False)

for i in range(len(A)):
#    B = np.matmul( A.transpose(), A)
    B=np.cov(A, rowvar=False)
    C = np.linalg.eig(B)
    D = abs(C[0])

    if (max(D)/min(D)>10E+12):
        j = np.argmin(D[2:5]) + 2
        # A_sparse = sparse.csr_matrix(A)
        # similarities = abs(cosine_similarity(A_sparse))
        # col_sum=np.sum(similarities,axis=0)
        # j=np.argmax(col_sum[2:5])+2
        A=np.delete(A, j, 1)
        print(j)

#also can output sparse matrices
similarities_sparse = cosine_similarity(A_sparse,dense_output=False)
print('pairwise sparse output:\n {}\n'.format(similarities_sparse))