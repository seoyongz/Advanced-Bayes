# Review of Advanced Bayes Class

<details>
<summary>Hierarchical Linear Model</summary>

### Model
- $y_j \sim \mathbb{R}^{n_j}$ : observation vector
- $X_j \sim \mathbb{R}^{n_j \times d}$ : design matrix 
- $\beta_j\in \mathbb{R}^d$ : subject-specific random effects
- $j=1, \ldots, m$ : subject index

$$y_j \sim \text{N}_{n_j}(X_j\beta_j,\ \sigma^2I_{n_j})\\
\beta_j \sim \text{N}_d(\mu_\beta,\ \sigma_\beta)$$

where $\sigma^2>0$, $\mu_\beta \in \mathbb{R}^d$, and $\Sigma_\beta \in \mathbb{R}^{d\times d}$ (positive definite)

### Priors
$$\mu_\beta \sim \text{N}_d(\xi,\ \Omega),\\

\sigma^2 \sim \text{Inv-}\chi^2(\nu,\ \tau^2),\\

\Sigma_\beta \sim \text{Inv-Wishart}_\rho(\Psi^{-1})$$


- $p(\beta\ |\ \sigma^2,\ \mu_\beta,\ \Sigma_\beta,\ y)$


- $p(\sigma^2\ |\ \beta,\ \mu_\beta,\ \Sigma_\beta,\ y)$


- $p(\mu_\beta\ |\ \beta,\ \sigma^2,\ \Sigma_\beta,\ y)$


- $p(\Sigma_\beta\ |\ \beta,\ \sigma^2,\ \mu_\beta,\ y)$



### Full conditional posterior distribution

</details>



<details>
<summary>Generalized Linear Model </summary>

### Model 


### Priors



### Sampling scheme



#### Independent MH algorithm


#### Random walk Metropolis algorithm


#### Data augmentation


</details>

<details>
<summary>Nonlinear Mixed Model </summary>

### Model 
$$y_{ij} = \frac{\beta_1 + u_i}{1+\exp\left\{-(\text{AGE}_{ij} - \beta_2)/\beta_3 \right\}}\\
u_i \sim \text{N}(0,\ \tau^2),\\
\epsilon_{ij} \sim \text{N}(0,\ \sigma^2)$$

### Priors
use



</details>


<details>
<summary>Nonparametric Regression Model </summary>

### Model 
$$y_i = \sin^3(2\pi x_i^3) + \epsilon_i,\\
\epsilon_i \sim \text{N}(0,\ 0.1^2)$$
Let $x_i = (2i-1)/1000,\ i=1, \ldots, n$ with $n=500$

### Priors


</details>

<details>
<summary>Gaussian Process Regression Model </summary>

### Model 
$$y_i = \mu(x_i) + \epsilon_i,\\
\epsilon_i \sim \text{N}(0,\ \sigma^2)$$
where $x_i \in \mathbb{R}^p$

### Priors
$$\mu \sim \text{GP}(0,\ k),\\

k(x,\ x') = \tau^2 \exp\left(-\frac{(x-x')^2}{l^2} \right),\\

\log(\sigma^2) \propto 1$$


### Posteior distribution
$$\begin{pmatrix} y \\ \bar{\mu} \\ \tilde{\mu} \end{pmatrix} |\ \sigma^2\sim \text{N}_{2n+m}\left(0,\ \begin{pmatrix} K(x,x)+\sigma^2I_n & K(x,x) & K(x, \tilde{x}) \\ K(x,x) & K(x,x) & K(x, \tilde{x}) \\ K(x, \tilde{x}) & K(x, \tilde{x}) & K(\tilde{x}, \tilde{x}) \end{pmatrix}\right)
$$

derive image

</details>


<details>

<summary>Finite Mixture Model </summary>

### Model 
page 10

### Priors

</details>




1. [Hierarchical Linear Models](https://www.notion.so/CH15-Hierarchical-Linear-Models-b34181ff98dd4ba085515bdcb1e80b4e)

2. [Generalized Linear Models](https://www.notion.so/CH16-Generalized-Linear-Models-d6f828054d614701acc0ba9aafbedf17)

3. [Basis Functional Models](https://www.notion.so/CH20-Basis-Functional-Model-185b7c93b7cc41d8aceecd66caf906a4)

4. [Gaussian Process Regression](https://www.notion.so/CH21-Gaussian-Process-Models-bd6f9e86ec9d4060960e138ff57fda0d)

5. [Finite Mixture Model](https://www.notion.so/CH22-Finite-Mixture-Models-e54a9682d707492f80a005d8a3084510)


