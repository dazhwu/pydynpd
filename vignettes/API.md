This document describes how to manipulate the outputs of the abond function. For instruction on what inputs are needed to run the abond function, and how to form these inputs, please refer to [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb).

## List of models returned

After regression, the abond function returns an "abond" object whose major property is "models" representing a list of models:

![Figure 1](https://raw.githubusercontent.com/dazhwu/pydynpd/main/vignettes/Images/list_models.svg)

Multiple models could be returned because users may include "?" in command string. For example, in Example 8 in [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb), the following comand string generates 5 models:

    command_str='n L(1:?).n w k | gmm(n, 2:3) pred(w k)' 
    mydpd = regression.abond(command_str, df, ['id', 'year']) 

What model(s) to be included in the list returned is based on the following rule:

-   If the user specifies THE model in the command string (i.e., she doesn't use "?" in command string), then this model is included in the list which contains only one model.

-   Otherwise, only models that statisfy all of the following conditions are included in the list returned:

    -   pass Hansen over-identification test (i.e., its P value should be greater than 5%)

    -   pass AR(2) test (i.e., its P value is greater than 5%)

    -   The P-value of Hansen over-identification test shouldn't be too high. In other words, it should be less than 99.99%. This is because a P value too close to 1 indicates a potential too-many-instrument issue.

    For example, in Example 8 in [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb), 5 models are generated but only 3 of them are returned because the other two fail to satisfy at least one of the conditions above.

    We can access each model on the list by its index which starts from 0. For example, the following code runs regression and then accesses the first model in the list returned:

        command_str='n L(1:?).n w k | gmm(n, 2:3) pred(w k)'
        mydpd = regression.abond(command_str, df, ['id', 'year'])
        m=mydpd.models[0]

In the example above, mydpd is an "abond" object, and its "models" property represents a list of models. On the other hand, m, as a dynamic panel model, is the first model on the list. Next, we discuss the properties of the dynamic panel model object.

## Dynamic penel model object

The list below shows the properties of a dynamic panel model object:

| Property         | Data Type                   | Meaning                                                                                                                                                                                 |
|------------------|:----------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| AR_list          | list                        | list of Arellano-Bond test objects                                                                                                                                                      |
| MMSC_LU          | dictionary                  | contains bic (Bayesian information criterion), aic (Akaike information criterion), and hqic (Hannan-Quinn information criterion) values suggested by Andrews, D. and Lu, B. (2001) (\*) |
| N                | int                         | number of individuals                                                                                                                                                                   |
| T                | int                         | number of time periods                                                                                                                                                                  |
| command_str      | str                         | command string of the model                                                                                                                                                             |
| final_xy_tables  | dictionary                  | each element of the dictionary is a Numpy ndarray that contains data of dependent variable or independent variables.                                                                    |
| hansen           | hansen test object          | hansen over identification test                                                                                                                                                         |
| avg_obs          | int                         | average number of observations per individual                                                                                                                                           |
| max_obs          | int                         | maximum number of observations per individual                                                                                                                                           |
| min_obs          | int                         | minimum number of observations per individual                                                                                                                                           |
| num_obs          | int                         | total number of observations                                                                                                                                                            |
| regression_table | pandas data frame           | contains all statistics in the regression table                                                                                                                                         |
| step_results     | list of step result objects |each object is related to a particular GMM step                                                                                                                                                                               |
| z_list           | Numpy ndarray               | instrument matrix                                                                                                                                                                       |

(\*) Andrews, D. and Lu, B. (2001). Consistent Model and Moment Selection Procedures for GMM Estimation with application to dynamic panel data models. Journal of Econometrics, 101(1):123--164.

### AR list property

This list contains two AR test (Arellano-Bond test) objects. For example, if m is a dynamic panel model object, then m.AR_list[0] refers to the AR test object for the first-order autocovariance of the residuals, and m.AR_list[1] for the second-order autocovariance.

#### AR test object

An AR test object has three properties:

| Property | Data Type | Meaning                                                               |
|----------|:----------|:----------------------------------------------------------------------|
| AR       | float     | AR test value                                                         |
| P_value  | float     | P value of the AR test                                                |
| lag      | int       | 1 for the first order test (AR(1)) and 2 for the second-order (AR(2)) |

For example, if you run the following code:
```{=Python}
    mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) |nolevel fod ', df, ['id', 'year'])

    m=mydpd.models[0]
    print("AR(2) test:", end=" ")
    print(m.AR_list[1].AR )
    print("P value: ", end=" ")
    print(m.AR_list[1].P_value)
```
the output will be

```{=html}
 Dynamic panel-data estimation, two-step difference GMM
 Group variable: id                               Number of obs = 611     
 Time variable: year                              Min obs per group: 4    
 Number of instruments = 36                       Max obs per group: 6    
 Number of groups = 140                           Avg obs per group: 4.36 
+------+------------+---------------------+------------+-----------+-----+
|  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
+------+------------+---------------------+------------+-----------+-----+
| L1.n | 0.0905525  |      0.1179995      | 0.7673969  | 0.4428456 |     |
| L2.n | -0.0400400 |      0.0397088      | -1.0083405 | 0.3132910 |     |
|  w   | -0.8379635 |      0.1197220      | -6.9992428 | 0.0000000 | *** |
|  k   | 0.6088287  |      0.0952578      | 6.3913796  | 0.0000000 | *** |
+------+------------+---------------------+------------+-----------+-----+
Hansen test of overid. restrictions: chi(32) = 37.921 Prob > Chi2 = 0.217
Arellano-Bond test for AR(1) in first differences: z = -1.01 Pr > z =0.314
Arellano-Bond test for AR(2) in first differences: z = -0.59 Pr > z =0.556

AR(2) test: -0.5892242871775409
P value:  0.5557108268348134
```
Note that in the general output, both AR and P values are rounded. For example in the output above, a value of -0.5892242871775409 is rounded to -0.59.

### MMSC_LU property

This property is presented when the abond function returns two or more models, so that there is a need to compare them using their bic, hqic, and aic values.

Example:

    command_str='n L(1:?).n w k | gmm(n, 2:3) pred(w k)'
    mydpd = regression.abond(command_str, df, ['id', 'year'])

    for i in range(0, len(mydpd.models)):
        print("model", end=" ")
        print(i+1, end=": bic= ")
        print(mydpd.models[i].MMSC_LU['bic'], end = "; hqic=")
        print(mydpd.models[i].MMSC_LU['hqic'], end="; aic=")
        print(mydpd.models[i].MMSC_LU['aic'])

The output is too long and therefore we only show below the section generated by the for-loop part:

    model 1: bic= -566.7507189938725; hqic=-290.96772885274686; aic=-86.12453121040255
    model 2: bic= -497.2812102841098; hqic=-256.13657126336705; aic=-73.43190220363662
    model 3: bic= -215.56036505263125; hqic=-120.02975924153786; aic=-39.93759993811322

### final_xy_tables property

This attribute is a dictionary. What are included in the dictionary depends on the transformation method used. For example, suppose m is a dynamic panel model object. If first difference method is used, the dictionary has two elements: m.final_xy_tables['Cx'] and m.final_xy_tables['Cy'] contain data of independent and dependent variables respectively. On the other hand, if forward orthogonal deviation transformation is used, then the dictionary has four elements: m.final_xy_tables['Cx'], m.final_xy_tables['Cy'], m.final_xy_tables['Diff_x'], and m.final_xy_tables['Diff_y']. The latter two elements are the first-difference matrices of independent and dependent variables respectively.

### hensen property

It is an object that has the following attributes:

| Property   | Data Type | Meaning                    | Example (suppose m is a dynamic model object)   |
|------------|:----------|:---------------------------|:--------------|
| test_value | float     | Hensen test value          |m.hensen.test_value |
| p_value    | float     | P value of the Hensen test |m.hensen.p_value |
| df         | int       | degree of freedom          |m.hensen.df |

### regression_table property

This property is a pandas dataframe with the following columns:

| column      | Meaning                                        |
|-------------|:-----------------------------------------------|
| variable    | name of the corresponding independent variable |
| coefficient | estimated coefficient                          |
| std_err     | corrected standard error                       |
| z_value     | z value                                        |
| p_value     | P value                                        |
| sig         | significance mark                              |

In the following example, we first run regression, and then use regression_table property to retrieve values in the regression table:

    mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) |nolevel fod ', df, ['id', 'year'])

    m=mydpd.models[0]
    print(m.regression_table)
    print("the coefficient of the first coefficient:")
    print(m.regression_table.iloc[0]['coefficient'])
    print("the p value of the first coefficient:")
    print(m.regression_table.iloc[0]['p_value'])

output:

     Dynamic panel-data estimation, two-step difference GMM
     Group variable: id                               Number of obs = 611     
     Time variable: year                              Min obs per group: 4    
     Number of instruments = 36                       Max obs per group: 6    
     Number of groups = 140                           Avg obs per group: 4.36 
    +------+------------+---------------------+------------+-----------+-----+
    |  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
    +------+------------+---------------------+------------+-----------+-----+
    | L1.n | 0.0905525  |      0.1179995      | 0.7673969  | 0.4428456 |     |
    | L2.n | -0.0400400 |      0.0397088      | -1.0083405 | 0.3132910 |     |
    |  w   | -0.8379635 |      0.1197220      | -6.9992428 | 0.0000000 | *** |
    |  k   | 0.6088287  |      0.0952578      | 6.3913796  | 0.0000000 | *** |
    +------+------------+---------------------+------------+-----------+-----+
    Hansen test of overid. restrictions: chi(32) = 37.921 Prob > Chi2 = 0.217
    Arellano-Bond test for AR(1) in first differences: z = -1.01 Pr > z =0.314
    Arellano-Bond test for AR(2) in first differences: z = -0.59 Pr > z =0.556

      variable  coefficient   std_err   z_value       p_value  sig
    0     L1.n     0.090552  0.118000  0.767397  4.428456e-01     
    1     L2.n    -0.040040  0.039709 -1.008341  3.132910e-01     
    2        w    -0.837964  0.119722 -6.999243  2.573495e-12  ***
    3        k     0.608829  0.095258  6.391380  1.643957e-10  ***
    the coefficient of the first coefficient:
    0.09055246263997319
    the p value of the first coefficient:
    0.442845568683905

### step_results property

This is a list of step result objects. For example, for a two-step GMM model, this list contains two objects: step_results[0] for results generated in step 1 and step_results[1] for those in step 2.

#### step_result object

Each step_result object has the following properties:

| Property | Data Type     | Meaning                                   |
|----------|:--------------|:------------------------------------------|
| W        | Numpy ndarray | Weighting matrix                          |
| beta     | Numpy ndarray | Estimated coefficients                    |
| vcov     | Numpy ndarray | Covariance matrix of coefficients         |
| std_err  | Numpy ndarray | Standard errors of coefficients           |
| residual | Numpy ndarray | one-column matrix that contains residuals |

Example:

If we run code below,

    mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) |nolevel fod ', df, ['id', 'year'])

    m=mydpd.models[0]
    print("Coefficients:")
    print(m.step_results[1].beta)  # for ste 2
    print("Standard errors of coefficients: ")
    print(m.step_results[1].std_err)

the output will be:

     Dynamic panel-data estimation, two-step difference GMM
     Group variable: id                               Number of obs = 611     
     Time variable: year                              Min obs per group: 4    
     Number of instruments = 36                       Max obs per group: 6    
     Number of groups = 140                           Avg obs per group: 4.36 
    +------+------------+---------------------+------------+-----------+-----+
    |  n   |   coef.    | Corrected Std. Err. |     z      |   P>|z|   |     |
    +------+------------+---------------------+------------+-----------+-----+
    | L1.n | 0.0905525  |      0.1179995      | 0.7673969  | 0.4428456 |     |
    | L2.n | -0.0400400 |      0.0397088      | -1.0083405 | 0.3132910 |     |
    |  w   | -0.8379635 |      0.1197220      | -6.9992428 | 0.0000000 | *** |
    |  k   | 0.6088287  |      0.0952578      | 6.3913796  | 0.0000000 | *** |
    +------+------------+---------------------+------------+-----------+-----+
    Hansen test of overid. restrictions: chi(32) = 37.921 Prob > Chi2 = 0.217
    Arellano-Bond test for AR(1) in first differences: z = -1.01 Pr > z =0.314
    Arellano-Bond test for AR(2) in first differences: z = -0.59 Pr > z =0.556

    Coefficients:
    [[ 0.09055246]
     [-0.04004001]
     [-0.83796355]
     [ 0.60882865]]
    Standard errors of coefficients: 
    [0.11799952 0.03970882 0.11972203 0.09525778]

### z_list

This is a two-dimentional numpy array which represents the GMM instrument matrix. For example, if m is a dynamic panel model object, then m.z_list is the instrument matrix used in the GMM process.
