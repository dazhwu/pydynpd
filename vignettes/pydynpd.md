---
author:
- Dazhong Wu, Jian Hua, and Feng Xu
bibliography: "ref.bib"
title: "pydynpd: A python package for dynamic panel model"
---

Department of Management, School of Business and Public Administration,
University of the District of Columbia

# Summary

We developed a python package (pydynpd) which implements all the
features in dynamic panel model with GMM (general method of moments).
These features include: (1) difference and system GMM, (2) one-step and
two-step estimators, (3) robust standard errors including the one
suggested by [@windmeijer2005finite], (4) Hansen over-identification test, (5)
Arellano-Bond test for autocorrelation, (6) time dummies, (7) allows
users to collapse instruments to reduce instrument proliferation issue,
and (8) a simple grammar for model specification. As far as we know, pydynpd is the first python package that allows researchers to estimate dynamic panel model.

What distinguishes pydynpd from any other dynamic panel model packages is its its innovative feature: the capability to search for models based on users' request, rather than just run the model specified by users as other packages do. To the best of our knowledge, there is no other econometric software/package that offers this feature, let alone dynamic panel model packages.

# Statement of need 
Over the past decade, dynamic panel model has become increasingly
popular in empirical studies. For example, researchers use dynamic panel
model to study the environmental impact of climate change ([@econometrics8030030]).
This is because many aspects of our social and natural systems are
inherently dynamic, and the GMM methods proposed by [@arellano1991some] and [@blundell1998initial] allow us to model the dynamics that
traditional static panel models are not able to capture.
Correspondingly, the growing popularity of dynamic panel model will
stimulate demand for the related packages in open source programs such
as R, Python, and Julia,

# Statement of field 
So far, there are several dynamic panel model packages in Stata and R.
Stata is a commercial software, while existing R packages have some
issues. For example, in our benchmark test R package panelvar ([@sigmund2021panel]) is more than 100 times slower than xtabond2 ([@roodman2009xtabond2]). On the other hand, R package plm ([@croissant2008panel])
takes more than twice as much time as xtabond2, and more importantly it
has calculation issue for system GMM. Due to these reasons, R
packages above are far less popular than xtabond2, according to citations they
have received.

Moreover, there is no python nor Julia package yet to estimate dynamic
panel model due to the complexity involved in implementation. Our
package contributes to the open source community because (1) it
implements all of the major features in the associated commercial packages in
Stata, (2) its innovative feature (as mentioned above) will stimulate similar or even more revolutionary features in the empirical computing community, and (3) though Python is interpreted, our package is almost as
fast as xtabond2 which was compiled as shown in figure below. This package will increase the usability of open source software in estimating dynamic panel models, because for a package to be attractive, it must be both accurate and fast. Moreover, unlike existing R
packages which rely heavily on R-specific components (that is a main
reason they are not fast), our code uses components common to any
programming language, making it easy to translate to R or Julia.

<p align = "center">
<img src = "./../Benchmark/images/Test_1.svg">
</p>
<p align = "center">
Fig.1 - Performance Comparision 
</p>

# The pydynpd package 

pydynpd is able to estimate the most complicated linear dynamic panel
models:

$$y_{it}=\sum_{j=1}^{p}\alpha_{j}y_{i,t-j}+\sum_{k=1}^{m}\sum_{j=0}^{q_{k}}\beta_{jk}r_{i,t-j}^{(k)}+\boldsymbol{\delta}\boldsymbol{d_{i,t}}+\boldsymbol{\gamma}\boldsymbol{s_{i,t}}+u_{i}+\epsilon_{it}$$

$$
\begin{equation}
\require{color}   y_{it}=\sum_{j=1}^\colorbox{yellow}{$p$}\alpha_{j}y_{i,t-j}+\sum_{k=1}^{m}\sum_{j=1}^\colorbox{yellow}{$q_{k}$}\beta_{jk}r_{i,t-j}^{(k)}+\boldsymbol{\delta}\boldsymbol{d_{i,t}}+\boldsymbol{\gamma}\boldsymbol{s_{i,t}}+u_{i}+\epsilon_{it}   
\end{equation}
$$
In the model above, $y_{i,t-j}$ ($j=1,2,\ldots,p$) denotes a group of
$p$ lagged dependent variables. $r_{i,t-j}^{(k)}$ represents a group of
$m$ endogenous variables other than lagged $y$. $\boldsymbol{d_{it}}$ is
a vector of predetermined variables which may potentially correlate with
past errors, $\boldsymbol{s_{it}}$ is a vector of exogenous variables,
and $u_{i}$ represents fixed effect. For illustration purpose, let's
consider a simpler model:




As lagged dependent variables such as $y_{i,t-1}$ are included as regressors, the
popular techniques in static panel models no longer produce consistent
results. Researchers have developed many methods to estimate dynamic
panel model. Essentially there are two types of GMM estimates,
difference GMM and system GMM. Just like other R and Stata packages, pydynpd fully implements these two methods.

Due to space limit, we focus here on general discussion of the package. A detailed statistical/technique description of exactly what pydynpd does can be found in the /vignettes folder on the github repository site. 

Consider the following equation:
$$\require{color} y_{it}=\sum_{j=1}^\colorbox{yellow}{$p$}\alpha_{j}y_{i,t-j}+\sum_{j=1}^\colorbox{yellow}{$q_{k}$}\beta_{j}r_{i,t-j}+{\delta}d_{i,t}+\gamma_{i,t}+u_{i}+\epsilon_{it}$$

The equation above is related to a group of models with different combinations of $p$ and $q_{k}$ values. Unless existing economic theory indicates what model to choose, researchers need to guess the values of $p$ and $q_{k}$ as highlighted in equation above. For example, if $p=2$ and $q_{k}=1$, then a specific model is formed:

$$\require{color} y_{it}=\alpha_{1}y_{i,t-1}+\alpha_{2}y_{i,t-2}+\beta_{j}r_{i,t-j}+{\delta}d_{i,t}+\gamma_{i,t}+u_{i}+\epsilon_{it}$$

Figure 2 shows how other packages work: a user needs to choose a specific model, then the system generates the corresponding instrument matrix and panel data with dependent/independent variables for the GMM process to produce regression results. An innovative feature of pydynpd is that it can also run in its "automatic" mode in which it doesn't require users to choose a particular model. Instead, users may let pydynpd search for the lags (e.g., $p$ and $q_{k}$) so that the corresponding models satisfy certain standards.
<p align = "center">
<img src = "Images/traditional.svg">
</p>
<p align = "center">
Fig.2 - The traditional way
</p>


Figure 3 shows how pydynpd's automatic mode works: a user indicates what values pydynpd needs to search for (e.g., the question marks in equation below), and then pydynpd tries all possible models, and returns models that pass dynamic models' specificiation tests (e.g., Hansen overidentification test and AR(2) test). Note that in Figure 3 the process named "traditional process" consistes of procedures included in the dotted box in Figure 2.

$$\require{color} y_{it}=\sum_{j=1}^\colorbox{yellow}{$?$}\alpha_{j}y_{i,t-j}+\sum_{j=1}^\colorbox{yellow}{?}\beta_{j}r_{i,t-j}+{\delta}d_{i,t}+\gamma_{i,t}+u_{i}+\epsilon_{it}$$



<p align = "center">
<img src = "Images/new_struct.svg">
</p>
<p align = "center">
Fig.3 - pydynpd's automatic mode
</p>






# References

