## Model Implementaion with Specific Example:

<details>
<summary>Hierarchical Linear Model</summary>

### Model
- $y_j \sim \mathbb{R}^{n_j}$ : observation vector
- $X_j \sim \mathbb{R}^{n_j \times d}$ : design matrix 
- $\beta_j\in \mathbb{R}^d$ : subject-specific random effects
- $j=1, \ldots, m$ : subject index

<div align="center">
$$y_j \sim \text{N}_{n_j}(X_j\beta_j,\ \sigma^2 I_n)$$
$$\beta_j \sim \text{N}_d(\mu_\beta,\ \sigma_\beta)$$
</div>
where $\sigma^2>0$, $\mu_\beta \in \mathbb{R}^d$, and $\Sigma_\beta \in \mathbb{R}^{d\times d}$ (positive definite)

### Priors

<div align="center">

$$\mu_\beta \sim \text{N}_d(\xi,\ \Omega),$$

$$\sigma^2 \sim \text{Inv-}\chi^2(\nu,\ \tau^2),$$

$$\Sigma_\beta \sim \text{Inv-Wishart}_\rho(\Psi^{-1})$$
</div>


</details>


<details>
<summary>Nonlinear Mixed Model </summary>

### Model 

<div align="center">

$$y_{ij} = \frac{\beta_1 + u_i}{ 1+\exp [-(\text{AGE}_{ij} - \beta_2)/\beta_3 ] },$$

$$u_i \sim \text{N}(0,\ \tau^2),$$

$$\epsilon_{ij} \sim \text{N}(0,\ \sigma^2)$$
</div>

### Priors

<div align="center">

$$p(\tau)\propto 1$$
</div>


</details>


<details>
<summary>Basis Functional Model </summary>

### Model 

<div align="center">

$$y_i = \sin^3(2\pi x_i^3) + \epsilon_i,$$

$$\epsilon_i \sim \text{N}(0,\ 0.1^2)$$
</div>

Let $x_i = (2i-1)/1000,\ i=1, \ldots, n$ with $n=500$

(a) Use **truncated power basis** with fixed $L=11$ interior uniform knots

(b) Use **polynomial radial basis** with fixed $L=11$ interior uniform knots

(c) Use **B-Spline basis** with fixed $L=11$ interior uniform knots

(d) Use **B-Spline basis** with $L\sim \text{Pois}(1)$, put the $g$-prior on the coefficents $\beta_H$ with $g=n$, 



</details>

<details>
<summary>Gaussian Process Regression Model </summary>

### Model 
<div align="center">

$$y_i = \mu(x_i) + \epsilon_i,$$

$$\epsilon_i \sim \text{N}(0,\ \sigma^2)$$
</div>

where $x_i \in \mathbb{R}^p$

### Priors
<div align="center">

$$\mu \sim \text{GP}(0,\ k),$$

$$k(x,\ x') = \tau^2 \exp\left(-\frac{(x-x')^2}{l^2} \right),$$

$$\log(\sigma^2) \propto 1$$
</div>


</details>


<details>

<summary>Finite Mixture Model </summary>

### Model 
Univariate location-scale mixture of Gaussians

<div align="center">

$$y_i\, |\, z_i \sim \text{N}(\mu_{z_i},\ \tau_{z_i}^2),$$

$$\text{P}(z_i=h) = \pi_h, \quad i=1, \ldots, n$$
</div>

### Priors
<div align="center">

$$(\pi_1,\ldots, \pi_H)\sim \text{Dirichlet}(a,\ldots, a),$$

$$\mu_h\,|\,\tau_h^2 \sim \text{N}(\mu_0,\ \kappa \tau_h^2),$$

$$\tau_h^2 \sim \text{Inv-Gamma}(a_\tau,\ b_\tau),\quad h=1, \ldots, H$$
</div>
</details>


### My summary (details about the model)
[Basis Functional Models](https://www.notion.so/CH20-Basis-Functional-Model-185b7c93b7cc41d8aceecd66caf906a4)

[Gaussian Process Regression](https://www.notion.so/CH21-Gaussian-Process-Models-bd6f9e86ec9d4060960e138ff57fda0d)


