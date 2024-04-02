# Review of Advanced Bayes Class

<details>
<summary>Hierarchical Linear Model</summary>

### Model
- $y_j \sim \mathbb{R}^{n_j}$ : observation vector
- $X_j \sim \mathbb{R}^{n_j \times d}$ : design matrix 
- $\beta_j\in \mathbb{R}^d$ : subject-specific random effects
- $j=1, \ldots, m$ : subject index
$\begin{align}
y_j &\sim \text{N}_{n_j}(X_j\beta_j,\ \sigma^2I_{n_j})\\
\beta_j&\sim \text{N}_d(\mu_\beta,\ \sigma_\beta)

\end{align}
where $\sigma^2>0$, $\mu_\beta \in \mathbb{R}^d$, and $\Sigma_\beta \in \mathbb{R}^{d\times d}$ (positive definite)

### Priors
$\begin{align}
\mu_\beta &\sim \text{N}_d(\xi,\ \Omega),\\

\sigma^2 &\sim \text{Inv-}\chi^2(\nu,\ \tau^2),\\

\Sigma_\beta &\sim \text{Inv-Wishart}_\rho(\Psi^{-1})

\end{align}$

- $p(\beta\ |\ \sigma^2,\ \mu_\beta,\ \Sigma_\beta,\ y)$


- $p(\sigma^2\ |\ \beta,\ \mu_\beta,\ \Sigma_\beta,\ y)$


- $p(\mu_\beta\ |\ \beta,\ \sigma^2,\ \Sigma_\beta,\ y)$


- $p(\Sigma_\beta\ |\ \beta,\ \sigma^2,\ \mu_\beta,\ y)$



### Full conditional posterior distribution

</details>


1. [Hierarchical Linear Models](https://www.notion.so/CH15-Hierarchical-Linear-Models-b34181ff98dd4ba085515bdcb1e80b4e)

2. [Generalized Linear Models](https://www.notion.so/CH16-Generalized-Linear-Models-d6f828054d614701acc0ba9aafbedf17)

3. [Basis Functional Models](https://www.notion.so/CH20-Basis-Functional-Model-185b7c93b7cc41d8aceecd66caf906a4)

4. [Gaussian Process Regression](https://www.notion.so/CH21-Gaussian-Process-Models-bd6f9e86ec9d4060960e138ff57fda0d)

5. [Finite Mixture Model](https://www.notion.so/CH22-Finite-Mixture-Models-e54a9682d707492f80a005d8a3084510)


