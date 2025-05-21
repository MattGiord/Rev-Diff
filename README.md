# Rev-Diff

MATLAB code for statistical analysis with continuous observations of bi-dimensional reversible diffusions.

Authors: Matteo Giordano, https://matteogiordano.weebly.com, and Kolyan Ray, https://kolyanray.wordpress.com.

This repository is associated with the article "Semiparametric Bernstein-von Mises theorems for reversible diffusions" by Matteo Giordano and Kolyan Ray. The abstract of the paper is:

"We establish a general semiparametric Bernstein--von Mises theorem for Bayesian nonparametric priors based on continuous observations in a 
periodic reversible multidimensional diffusion model. We consider a wide range of functionals satifying an approximate linearization condition, 
including several nonlinear functionals of the invariant measure. Our result is applied to Gaussian and Besov-Laplace priors, showing these 
can perform efficient semiparametric inference and thus justifying the corresponding Bayesian approach to uncertainty quantification. 
Our theoretical results are illustrated via numerical simulations."

This repository contains the MATLAB code to replicate the numerical simulation study presented in the article. It contains the following:

1. RevDiff.m, code to generate the observations and perform Bayesian nonparametric and semiparametric inference using a truncated Gaussian series prior.
2. RevDiffCoverage.m, code to perform the coverage study via repeated experiments.
3. CosCos.m, auxiliary function for the Fourier basis.
4. CosSin.m, auxiliary function for the Fourier basis.
5. SinCos.m, auxiliary function for the Fourier basis.
6. SinSin.m, auxiliary function for the Fourier basis.
7. GradCosCos.m, auxiliary function for the Fourier basis.
8. GradCosSin.m, auxiliary function for the Fourier basis.
9. GradSinCos.m, auxiliary function for the Fourier basis.
10. GradSinSin.m, auxiliary function for the Fourier basis.

For questions or for reporting bugs, please e-mail Matteo Giordano (matteo.giordano@unito.it).

Please cite the following article if you use this repository in your research: Giordano, M., and Ray, S. (2025). Semiparametric Bernstein-von Mises theorems for reversible diffusions.
