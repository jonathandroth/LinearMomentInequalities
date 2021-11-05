# LinearMomentInequalities

This repository contains Matlab code for implementing the tests proposed in [Andrews, Roth, and Pakes (ARP)](https://jonathandroth.github.io/assets/files/arp-draft.pdf), "Inference for Linear Conditional Moment Inequalities".
The paper proposes tests for moment inequalities of the form <img src="http://www.sciweavers.org/tex2img.php?eq=E%5BY_i%28%5Cbeta_0%29%20-%20X_i%20%5Cdelta%20%7C%20Z_i%5D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="E[Y_i(\beta_0) - X_i \delta | Z_i]" width="146" height="19" />, where <img src="https://bit.ly/3wlcnvx" align="center" border="0" alt="\beta_0" width="22" height="19" /> is the target parameter, <img src="https://bit.ly/3COT28l" align="center" border="0" alt="\delta" width="14" height="15" /> is a nuisance parameter, and <img src="https://bit.ly/3kaxZFN" align="center" border="0" alt="X_i" width="21" height="18" /> is non-random conditional on <img src="https://bit.ly/3k9YbAl" align="center" border="0" alt="Z_i" width="21" height="18" />.

The package contains the following main functions: 

* **hybrid_test_fn** implements the hybrid test, our preferred method
* **conditional_test_fn** implements the conditional test
* **lf_critical_value_fn** calculates the least-favorable critical value
* **etahat_fn** calculates the test statistic for the LF test (the same test stat is used for the conditional/hybrid tests, but done internally)
* **lf_ci_for_linear_params** calculates an LF CI for the case where the target parameter also enters the moments linearly.
* **conditional_variance_fn** estimates the conditional variance using the method of [Abadie, Imbens, and Zheng (2014)](https://www.tandfonline.com/doi/full/10.1080/01621459.2014.928218)


The file **LinearMomentInequalities_Examples.m** contains examples for how to use all of these functions. It also illustrates some computational shortcuts that can be taken when the target parameter also enters the moments linearly. 
