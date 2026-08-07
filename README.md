# Causal Regularization for Robust Prediction



## 📌 Abstract

Prediction is one of the most important uses of statistical methods. However, predictive systems typically struggle when the environment in which they were trained changes. The central research question is to identify and estimate stable models when only heterogeneous and potentially unknown environments are available. Existing approaches like Causal Dantzig [(Rothenh¨ausler et al.,2019)](https://arxiv.org/pdf/1706.06159) and Causal Regularization [(Kania and Wit, 2025)](https://arxiv.org/pdf/2205.01593) directly identify the causal parameter as the unique solution with equal residual moments across environments. We relax assumptions to obtain a robust best linear predictor without requiring an unperturbed environment or linear shifts. Causal interpretation is recovered under the inner-product invariance condition and when the conditional mean of the target given its causal parents remains invariant across environments.

We generalize Causal Regularization to a framework that uses two environments drawn from a common class $\mathcal{U}$ and measures the differences between them through a generalized risk difference induced by changes in second-order moments of the covariates and response. In this setting, the risk difference may be non-convex and may exhibit saddle-point geometry. We show that the parameter $\beta^\star$ is characterized as a stationary point of this risk difference. By defining a sieve of out-of-sample distributions, we derive a worst-risk decomposition over increasing classes of environments and obtain a closed-form estimator ${\beta_\gamma}$ that interpolates between OLS and the Regularization limit solution. 

Our results show that the proposed estimator can achieve prediction stability over a broad class of perturbed environments, even when trained on a limited set. The framework separates two goals: robustness as a stable prediction under weak assumptions and causality as a stronger property obtained with additional invariance structure.

---

## 🚀 Project Setup

All code is implemented in **R**.

The project is organized into two main folders:

* `functions/` — Contains all core R functions used across the project (e.g.,moments, risk computation).
* `scripts/` — Contains the main R scripts to run experiments and generate results. These scripts **source** the files from `functions/`.

## ✉️ Contact

For questions or comments about the code, please contact the authors.
