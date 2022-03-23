Department of Management, School of Business and Public Administration, University of the District of Columbia

# Summary {#summary .unnumbered}

Package pydynpd implements all the features needed to estimate dynamic panel models with GMM (general method of moments). These features include: (1) difference and system GMM, (2) one-step and two-step estimators, (3) robust standard errors including the one suggested by Windmeijer (2005), (4) Hansen over-identification test, (5) Arellano-Bond test for autocorrelation, (6) time dummies, and (7) collapse instruments to reduce instrument proliferation issue. To the best of our knowledge, pydynpd is the first python package to estimate dynamic panel model.

# Statement of need {#statement-of-need .unnumbered}

Over the past decade, dynamic panel model has become increasingly popular in empirical studies. For example, researchers use dynamic panel model to study the environmental impact of climate change (citation). This is because many aspects of our social and natural systems are inherently dynamic, and the GMM methods proposed by \... and \... allow us to model the dynamics that traditional static panel models are not able to capture. Correspondingly, the growing popularity of dynamic panel model will stimulate demand for the related packages in open source programs such as R, Python, and Julia,

# Statement of field  {#statement-of-field .unnumbered}

So far, there are several dynamic panel model packages in Stata and R. Stata is a commercial software, while existing R packages have many issues. For example, in our benchmark test R package panelvar is more than 100 times slower than Stata package xtabond2. On the other hand, R package plm takes more than twice as much time as xtabond2, and more importantly it has calculation issue for system GMM. On the other hand, there is no python package nor Julia package yet to estimate dynamic panel model due to the complexity involved in implementation. Our package contributes to the open source community because (1) it implements all of the major features in the commercial packages in Stata, and (2) though Python is interpretated, our package is almost as fast as xtabond2 which was compiled. Unlike existing R packages which rely heavily on R-specific components (that is a main reason they are not fast), our code uses static arraies, making it easy to translate to R or Julia.

# The pyndynpd package  {#the-pyndynpd-package .unnumbered}

Package pyndynpd is able to estimate dynamic panel models that take a form as follows:

$$y_{it}=\sum_{j=1}^{p}\alpha_{j}y_{i,t-j}+\sum_{k=1}^{m}\sum_{j=0}^{q_{k}}\beta_{jk}r_{i,t-j}^{(k)}+\boldsymbol{\delta}\boldsymbol{d_{i,t}}+\boldsymbol{\gamma}\boldsymbol{s_{i,t}}+u_{i}+\epsilon_{it}\label{eq:typical_model}$$ In the model above, $y_{i,t-j}$ ($j=1,2,\ldots,p$) denotes a group of $p$ lagged dependent variables. $r_{i,t-j}^{(k)}$ represents a group of $m$ endogeneous variables other than lagged $y$. $\boldsymbol{d_{it}}$ is a vector of predetermined variables which may potentially correlate with past errors, $\boldsymbol{s_{it}}$ is a vector of exogenous variables, and $u_{i}$ represents fixed effect. For illustration purpose, let's consider a basic form of dynamic panel model:

$$y_{it}=\alpha_{1}y_{i,t-1}+\delta d_{i,t}+u_{i}+\epsilon_{it}\label{eq:simple_model}$$

As lagged dependent variable $y_{i,t-1}$ is included as regressor, the popular techniques in static panel models, such as fixed-effect and first-difference estimators, no longer produce consistent results. Researchers have developed many methods to estimate dynamic panel model. Essentially there are two types of GMM estimates, difference GMM and system GMM.

## Difference GMM {#difference-gmm .unnumbered}

Difference GMM was developed by [@arellano1991some]. The first step in the process is to eliminate the fixed-effect term $u_{i}$. First differencing Eq [\[eq:typical_model\]](#eq:typical_model){reference-type="ref" reference="eq:typical_model"} yields:

$$\Delta y_{it}=\alpha_{1}\Delta y_{i,t-1}+\delta\Delta d_{i,t}+\Delta\epsilon_{it}\label{eq: FD}$$

In the model above, $\Delta y_{i,t-1}$correlates with $\Delta\epsilon_{i,t}$ because $\Delta y_{i,t-1}=y_{i,t-1}-y_{i,t-2}$, $\Delta\epsilon_{i,t}=\epsilon_{i,t}-\epsilon_{i,t-1}$, and $y_{i,t-1}$ is affected by $\epsilon_{i,t-1}$. As a result, estimating Eq [\[eq: FD\]](#eq: FD){reference-type="ref" reference="eq: FD"} directly produces inconsistent result. Instrumental variables are used to solve the issue. [@arellano1991some]suggest to use all lagged $y$ dated $t-2$ and earlier (i.e., $y_{i,1}$, $y_{i,2}$,\..., $y_{i,t-2}$) as instruments for $\Delta y_{i,t-1}$. Similarly, the instruments for predetermined variable $\Delta d_{it}$ include $d_{i,1}$, $d_{i,2}$,\..., $d_{i,t-1}$. Let $z_{i}$ be the instrument variable matrix for individual i:

$$z_{i}=\left[\begin{array}{ccccccccccccccccccc}
y_{i1} & 0 & 0 & 0 & 0 & 0 & \ldots & 0 & d_{i1} & d_{i2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \ldots & 0\\
0 & y_{i1} & y_{i2} & 0 & 0 & 0 & \ldots & 0 & 0 & 0 & d_{i1} & d_{i2} & d_{i3} & 0 & 0 & 0 & 0 & \ldots & 0\\
 &  &  & \vdots &  &  & \ldots &  &  &  &  &  &  & \vdots &  &  &  & \ldots & 0\\
0 & 0 & 0 & 0 & 0 & 0 & \ldots & y_{i,T-2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \ldots & d_{i,T-1}
\end{array}\right]\label{eq:Z_diff}$$

Difference GMM is based on the moment condition $E(z_{i}^{\prime}\Delta\epsilon_{i})=0$ where $\Delta\epsilon_{i}=(\Delta\textrm{\ensuremath{\epsilon_{i2}}, }\Delta\epsilon_{i3}\textrm{, }...,\Delta\epsilon_{iT})^{\prime}$ and $z_{i}$ is the instrument variable matrix. Applying this moment condition to sample data, we have $(1/N)\sum_{i=1}^{N}z{}_{i}^{\prime}(\Delta y_{i}-\theta\Delta x_{i})=0$ where $\theta=(\alpha_{1},\delta)'$ and $\Delta x_{i}=(\Delta y_{i,t-1},\Delta d_{it})$ for t=3, \... T. When the number of instruments is greater than the number of independent variables, the moment condition is overidentified and in general there is no $\theta$ available to satisfy the moment condition. Instead, we look for a $\theta$ to minimize moment condition. That is:

$$\hat{\theta}_{gmm}=\arg\min_{\theta}\left(\frac{1}{N}\sum_{i=1}^{N}(\Delta y_{i}-\theta\Delta x_{i})^{\prime}z_{i}\right)W\left(\frac{1}{N}\sum_{i=1}^{N}z^{\prime}{}_{i}(\Delta y_{i}-\theta\Delta x_{i})\right)$$

where W is the weighting matrix of the moments. There are two popularly used weighting matrixes. In a one-step GMM estimate, the weighting matrix is

$$W_{1}=\left(\frac{1}{N}Z^{\prime}H_{1}Z\right)^{-1}$$

where matrix H has twos in the main diagnols, minus ones in the first subdiagnols, and zeros elsewhere:

$$H_{1}=\left[\begin{array}{cccccc}
2 & -1 & 0 & 0 & \ldots & 0\\
-1 & 2 & -1 & 0 & \ldots & 0\\
0 & \ddots & \ddots & \ddots & \ddots & \vdots\\
\vdots & \ddots & -1 & 2 & -1 & 0\\
0 & \ddots & 0 & -1 & 2 & -1\\
0 & \ldots & 0 & 0 & -1 & 2
\end{array}\right]$$

\<mention when one-step is not good\>On the other hand, in a two-step GMM estimate, the weighting matrix is

$$W_{2}=\left(\frac{1}{N}Z^{\prime}H_{2}Z\right)^{-1}$$

where $H_{2}=\Delta\hat{\epsilon}\Delta\hat{\epsilon}^{\prime}$ and $\Delta\hat{\epsilon}$ is the residual from one-step GMM.

## System GMM {#system-gmm .unnumbered}

Compared with difference GMM, sytem GMM adds additional moment conditions, resulting in more instruments:

$$z_{i}=\left[\begin{array}{cccccccccccccccc|cccccccc}
y_{i1} & 0 & 0 & 0 & 0 & 0 & \ldots & 0 & d_{i1} & d_{i2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \ldots & 0 & 0 & 0 & 0 & \ldots & 0\\
0 & y_{i1} & y_{i2} & 0 & 0 & 0 & \ldots & 0 & 0 & 0 & d_{i1} & d_{i2} & d_{i3} & 0 & 0 & 0 & 0 & \ldots & 0 & 0 & 0 & 0 & \ldots & 0\\
 &  &  & \vdots &  &  & \ldots &  & \vdots &  &  &  &  & \vdots & \ldots & \vdots &  & \ldots & 0 &  &  &  & \ddots & 0\\
0 & 0 & 0 & 0 & 0 & 0 & \ldots & y_{i,T-2} & 0 & 0 & 0 & 0 & 0 & 0 & \ldots & d_{i,T-1} & 0 & \ldots & 0 & 0 & 0 & 0 & 0 & 0\\
\hline 0 & \ldots & 0 &  &  &  &  &  &  &  &  &  &  &  & 0 & \Delta y_{i2} & 0 & \ldots & 0 & \Delta d_{i3} & 0 &  & 0\\
\vdots &  &  &  &  &  &  &  &  &  &  &  &  &  &  & 0 & 0 & \Delta y_{i3} & \ldots & 0 & 0 & \Delta d_{i4} &  & 0\\
\vdots &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots &  & \ddots\\
0 & \ldots &  &  &  &  &  &  &  &  &  &  &  &  & \ldots & 0 & 0 & 0 & \ldots & \Delta y_{i,T-1} & 0 & 0 & \ldots & \Delta y_{i,T}
\end{array}\right]\label{eq: Z_sys}$$

$$\hat{\theta}_{gmm}=\arg\min_{\theta}\left(\frac{1}{N}\sum_{i=1}^{N}(\widetilde{y}_{i}-\theta\widetilde{x_{i}})^{\prime}z_{i}\right)W\left(\frac{1}{N}\sum_{i=1}^{N}z^{\prime}{}_{i}(\widetilde{y}-\theta\widetilde{x_{i}})\right)$$

where $$\widetilde{y}=\left(\begin{array}{c}
\Delta y_{i}\\
\hline y_{i}
\end{array}\right)\textrm{ and }\widetilde{x_{i}}=\left(\begin{array}{c|c}
\Delta x_{i} & 0\\
\hline x_{i} & 1
\end{array}\right)$$

## Robust estimation of coefficients' covariance {#robust-estimation-of-coefficients-covariance .unnumbered}

## Specification Test {#specification-test .unnumbered}

### Error serial correlation test {#error-serial-correlation-test .unnumbered}

Second-order serial correlation test: if $\epsilon_{it}$ in Eq [\[eq:simple_model\]](#eq:simple_model){reference-type="ref" reference="eq:simple_model"} is serially correlated, GMM estimates are no longer consistent. In a first-differenced model （e.g., Eq [\[eq: FD\]](#eq: FD){reference-type="ref" reference="eq: FD"}), to test whether epsilon $\epsilon_{i,t-1}$ is correlated with epsilon$\epsilon_{i,t-2}$, the second-order autocovariance of the residuals, $\textrm{AR(2)}$, is calculated as:

$$AR(2)=\frac{b_{0}}{\sqrt{b_{1}+b_{2}+b_{3}}}\textrm{ where}$$

$$b_{0}=\sum_{i=1}^{N}\Delta\hat{\hat{\epsilon}}_{i}^{\prime}L_{\Delta\hat{\hat{\epsilon}}}^{2}$$

$$b_{1}=\sum_{i=1}^{N}\text{\ensuremath{L_{\Delta\hat{\hat{\epsilon}}_{i}^{\prime}}^{2}H_{2}}\ensuremath{L_{\Delta\hat{\hat{\epsilon}}_{i}}^{2}}}$$ $$b_{2}=\textrm{-}2\left(\sum_{i=1}^{N}L_{\Delta\hat{\hat{\epsilon}}_{i}^{\prime}}^{2}x_{i}\right)\left[\left(\sum_{i=1}^{N}x_{i}^{\prime}z_{i}\right)W_{2}\left(\sum_{i=1}^{N}z_{i}^{\prime}x_{i}\right)\right]^{-1}\left(\sum_{i=1}^{N}x_{i}^{\prime}z_{i}\right)W_{2}\left(\sum_{i=1}^{N}z_{i}^{\prime}H_{2}L_{\Delta\hat{\hat{\epsilon}}_{i}}^{2}\right)$$

$$b_{3}=\left(\sum_{i=1}^{N}L_{\Delta\hat{\hat{\epsilon}}_{i}^{\prime}}^{2}x_{i}\right)\hat{V}_{\hat{\hat{\theta}}}\left(\sum_{i=1}^{N}x_{i}^{\prime}L_{\Delta\hat{\hat{\epsilon}}_{i}}^{2}\right)$$

### Hansen overidentification test {#hansen-overidentification-test .unnumbered}

Hansen overidentification test is used to check if instruments are exogeneous. Under the null hypothesis that instruments are valid, test statistic, $S$, should be close to zero:

$$S=\left(\sum_{i=1}^{N}\Delta\hat{\hat{\epsilon}}_{i}^{\prime}z_{i}\right)W_{2}\left(\sum_{i=1}^{N}z_{i}^{\prime}\Delta\hat{\hat{\epsilon}}_{i}\right)$$

# Handling instrument proliferation issue {#handling-instrument-proliferation-issue .unnumbered}

Difference GMM and system GMM may generate too many instruments, which causes several problems (citation). Package pydynpd allows users to reduce the number of instruments in two ways. First, users can control the number of instruments in command string. For example, $\textrm{gmm(w, 2 3)}$ states that only $n_{t-2}$ and $n_{t-3}$ are used as instruments, rather than all lagged $n$ dated $t-2$ and earlier. Second, users can choose to collapse the instrumental variable matrix. For example, if collapsed, matrix as in Eq [\[eq: Z_sys\]](#eq: Z_sys){reference-type="ref" reference="eq: Z_sys"} is changed to:

$$z_{i}=\left[\begin{array}{cccccccccc|cc}
y_{i1} & 0 & 0 & \ldots & 0 & d_{i1} & d_{i2} & 0 & \ldots & 0\\
y_{i1} & y_{i2} & 0 & \ldots & 0 & d_{i1} & d_{i2} & d_{i3} & \ldots & 0\\
\vdots & \vdots & \ddots & \ldots & \vdots &  &  &  & \ddots & \vdots\\
y_{i1} & y_{i2} & y_{i3} & \ldots & y_{i,T-2} & d_{i1} & d_{i2} & d_{id} & \ldots & d_{i,T-1}\\
\hline 0 &  &  & \ldots & 0 & 0 & 0 & 0 & \ldots & 0 & \Delta y_{i2} & \Delta d_{i3}\\
0 & 0 &  &  & 0 & 0 & 0 & 0 & \ldots & 0 & \Delta y_{i3} & \Delta d_{i4}\\
\vdots &  & \ddots &  & \vdots &  &  & \vdots &  & \vdots & \vdots & \vdots\\
0 & 0 & 0 & \ldots & 0 & 0 &  & 0 & \ldots & 0 & \Delta y_{i,T-1} & \Delta d_{iT}
\end{array}\right]$$

This change dramatically reduces the number of instruments. Intuitively, the number of instruments is positively associated with the width of the matrix above.

# References  {#references .unnumbered}
