# LinearMomentInequalities

This repository contains Matlab code for implementing the tests proposed in [Andrews, Roth, and Pakes (ARP)](https://jonathandroth.github.io/assets/files/arp-draft.pdf), "Inference for Linear Conditional Moment Inequalities".
The paper proposes tests for moment inequalities of the form E[Y_i(beta) - X_i delta | Z_i] <=0, where beta is the target parameter, delta is a nuisance parameter, and X_i is non-random conditional on Z_i. 

The package contains the following main functions: 

* **hybrid_test_fn** implements the hybrid test, our preferred method
* **conditional_test_fn** implements the conditional test
* **lf_critical_value_fn** calculates the least-favorable critical value
* **etahat_fn** calculates the test statistic for the LF test (the same test stat is used for the conditional/hybrid tests, but done internally)
* **lf_ci_for_linear_params** calculates an LF CI for the case where the target parameter also enters the moments linearly.
* **conditional_variance_fn** estimates the conditional variance using the method of [Abadie, Imbens, and Zheng (2014)](https://www.tandfonline.com/doi/full/10.1080/01621459.2014.928218)


The file **LinearMomentInequalities_Examples.m** contains examples for how to use all of these functions. It also illustrates some computational shortcuts that can be taken when the target parameter also enters the moments linearly. 
