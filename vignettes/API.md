This document describes how to manipulate the outputs of the abond function. For instruction on what inputs are needed to run the abond function, please refer to [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb).

## List of models returned

After regression, the abond function returns an "abond" object whose major property is "models" representing a list of models:

![Figure 1](https://raw.githubusercontent.com/dazhwu/pydynpd/main/vignettes/Images/list_models.svg)

Multiple models could be returned because users may include "?" in command string. For example, in Example 8 in [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb), the following comand string generates 5 models: ··· command_str='n L(1:?).n w k \| gmm(n, 2:3) pred(w k)' mydpd = regression.abond(command_str, df, ['id', 'year']) ···

What model(s) to be included in the list returned is based on the following rule:

-   If the user specifies THE model in the command string (i.e., she doesn't use "?" in command string), then this model is included in the list which contains only one model.

-   Otherwise, only models that statisfy all of the following conditions are included in the list returned:

    -   pass Hansen over-identification test (i.e., its P value should be greater than 5%)

    -   pass AR(2) test (i.e., its P value is greater than 5%)

    -   The P-value of Hansen over-identification test shouldn't be too high. In other words, it should be less than 99.99%. This is because a P value too close to 1 indicates a potential too-many-instrument issue.

    For example, in Example 8 in [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb), 5 models are generated but only 3 of them are returned because the other two fail to satisfy at least one of the conditions above.

    We can access each model on the list by its index which starts from 0. For example, the following code runs regression and then accesses the first model in the list returned:

        command_str='n L(1:?).n w k \| gmm(n, 2:3) pred(w k)'
        mydpd = regression.abond(command_str, df, ['id', 'year'])
        m=mydpd.models[0]

In the example above, mydpd is an "abond" object, and its "models" property represents a list of models. On the other hand, m, as a dynamic panel model, is the first model on the list. Next, we discuss the properties of the dynamic panel model object.

## Dynamic penel model object

The list below shows the properties of a dynamic panel model object:

| Property         | Data Type                   | Meaning                                               |
|------------------|:------------------|:----------------------------------|
| AR_list          | list                        | list of Arellano-Bond test objects                    |
| MMSC_LU          | dictionary                  | contains bic, aic, and hqic values                    |
| N                | int                         | number of individuals                                 |
| T                | int                         | number of time periods                                |
| command_str      | str                         | command string of the model                           |
| final_xy_tables  | dictionary                  | contains data of dependent and independent variables. |
| hansen           | hansen test object          | hansen over identification test                       |
| avg_obs          | int                         | average number of observations per individual         |
| max_obs          | int                         | maximum number of observations per individual         |
| min_obs          | int                         | minimum number of observations per individual         |
| num_obs          | int                         | total number of observations                          |
| regression_table | pandas data frame           | contains all statistics in the regression table       |
| step_results     | list of step result objects | discussed                                             |
| z_list           | Numpy ndarray               | instrument matrix                                     |

### AR list property

This list contains two AR test (Arellano-Bond test) objects. For example, if m is a dynamic panel model object, then m.AR_list[0] refers to the AR test object for the first-order autocovariance of the residuals, and m.AR_list[1] for the second-order autocovariance.

#### AR test object

An AR test object has three properties:

| Property | Data Type | Meaning                                                               |
|---------------|:--------------|:------------------------------------------|
| AR       | float     | AR test value                                                         |
| P_value  | float     | P value of the AR test                                                |
| lag      | int       | 1 for the first order test (AR(1)) and 2 for the second-order (AR(2)) |

For example, if you run the following code:

    mydpd = regression.abond('n L(1:2).n w k  | gmm(n, 2:4) gmm(w, 1:3)  iv(k) |nolevel fod ', df, ['id', 'year'])

    m=mydpd.models[0]
    print("AR(2) test:", end=" ")
    print(m.AR_list[1].AR )
    print("P value: ", end=" ")
    print(m.AR_list[1].P_value)

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
Note that in the general output, both AR and P values are rounded. For example, a value of -0.5892242871775409 is rounded to -0.59.

### MMSC_LU

Quarto enables you to weave together content and executable code into a finished document. To learn more about Quarto see <https://quarto.org>.

## Running Code

When you click the **Render** button a document will be generated that includes both content and the output of embedded code. You can embed code like this:

```{python}
1 + 1
```

You can add options to executable code like this

```{python}
#| echo: false
2 * 2
```

The `echo: false` option disables the printing of code (only output is displayed).
