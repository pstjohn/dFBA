# Dynamic Flux Balance analysis in Python

This package implements a simple dynamic flux balance algorithm in Python. The algorithm is similar to the one presented in *T. J. Hanly and M. A. Henson, “Dynamic flux balance modeling of microbial co-cultures for efficient batch fermentation of glucose and xylose mixtures,” Biotechnol. Bioeng., vol. 108, no. 2, pp. 376–385, Oct. 2010.*.

At each integration point, a LP is solved to find the optimum fluxes through the biological network. These fluxes are passed to an ODE integrator to find how external metabolite concentrations change in time. A user-provided bounds function defines the upper and lower bounds of reaction fluxes as a function of current metabolite concentrations, allowing for Michealis-Menten or other customized enzyme kinetics.

Everything is written in Cython, wrapping the C libraries GPLK for LP solutions and SUNDIALS for the ODE integrator.

## Dependencies

* GLPK v. 4.57 (Homebrew installed)
* SUNDIALS v. 2.6.2 
* opencobra/cobrapy
* Cython
