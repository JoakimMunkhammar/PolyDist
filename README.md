polydfit returns non-normalized polynomial coefficients for polynomial density distribution approximation for the data set in vector x within limits a and b.  (The function is similar to e.g. weib for Weibull parameter fit # in Scipy)

To get a proper PDF normalize polynomial function so that when integrated between a and b it equals one.

The function is based on a method of moments approach for generating  polynomial probability density distribution approximations, which should be referenced as:

[1] Munkhammar J, Mattsson L, Rydén J, (2017), Polynomial probability 
distribution estimation using the method of moments. PLoS ONE 
12(4): e0174573. doi:10.1371/journal.pone.0174573

Written by Dr. Joakim Munkhammar, PhD, and Mahmoud Shepero, Uppsala University.
Date: 2019/04/05

