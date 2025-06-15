# CRUSH, Causal Regularization Under Shifted Heterogeneity


This repository contains the official implementation for the thesis **CRUSH – Causal Regularization Under Shifted Heterogeneity**, which explores

**CRUSH – Causal Regularization Under Shifted Heterogeneity** is built on the principles of **invariance** and **out-of-sample guarantees**. It extends the **Causal Regularization** framework by [Kania et al., 2023](https://arxiv.org/abs/2302.03012), which leverages these ideas to identify stable causal relationships from observational data. CRUSH introduces a novel extension that handles **two distinct data shifts**, improving robustness across heterogeneous environments.


## 📌 Abstract

*Causal inference seeks to uncover the cause and effect relationships between variables, often using graphs to represent dependencies and pathways. This thesis explores causal analysis followed by advanced approaches, including Causal Regularization, which integrates invariant causal prediction and out-of-sample guarantees for both observational and shifted environments. A central contribution of this work is the novel extension of the Causal Regularization method to incorporate two shifted environments, utilizing invariant prediction to develop models that maintain stable predictive relationships across diverse data distributions.*

---

## ⚙️ Project Setup

All code is implemented in **R**.

The project is organized into two main folders:

* `functions/` — Contains all core R functions used across the project (e.g.,moments, risk computation).
* `scripts/` — Contains the main R scripts to run experiments and generate results. These scripts **source** the files from `functions/`.
