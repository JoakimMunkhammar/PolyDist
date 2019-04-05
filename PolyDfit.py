#!/usr/bin/env python3
#
# polyDfit returns non-normalized polynomial coefficients for polynomial
#   density distribution approximation for the data set in vector x 
#   within limits [a,b]. (The function is similar in design to e.g. weib 
#   for Weibull parameter fit in Scipy)
#
#   Note that NumPy and SciPy has to be imported for this to work properly.
#
#   The function is based on a method of moments approach for generating 
#   polynomial probability density distribution approximations, which should
#   be referenced as:
#
#   [1] Munkhammar J, Mattsson L, Rydén J (2017). Polynomial probability 
#   distribution estimation using the method of moments. PLoS ONE 
#   12(4): e0174573. doi:10.1371/journal.pone.0174573
#
#   The output coefficients should be normalized by dividing them by 
#   the integral of the polynomial between a and b.
#   
#   N is the polynomial degree, and a,b are free parameters inherent to the
#   method. These typically represent endpoints of the distribution interval 
#   in consideration, since this distribution approximation routine will 
#   approximate a probability density distribution on the interval [a,b].
#
#   Too high polynomial order can cause the matrix inversion to become
#   close to singular, if so, use lower order.
#
#   For best results the polynomial approximation probability distributions 
#   should be tested for goodness-of-fit prior to use, e.g. by visual 
#   inspection, Kolmogorov-Smirnov tests, AIC or similar.
#
#   Written by Dr. Joakim Munkhammar, PhD, and Mahmoud Shepero, Uppsala 
#   University.
#
#   Date: 2019/04/05
#

import numpy as np

def PolyDfit(X,a,b,N):

    # Finding the values in the range [a,b]
    idx = np.logical_and(X>a,X<b)
    data = X[idx]
    dataLength = data.shape[0]

    # Calculating the moments from 0 to N    
    moments = np.zeros((dataLength,N))
    for i in range(N):
        moments[:,i] = data ** (i) 

    # Summing up the moments
    moments2 = np.mean(moments, 0) 

    # Setting up the entries of the polynomial moment matrix
    entriesRng = np.arange(1, 2*N, dtype=np.float64)
    entries = ((b**entriesRng) - (a**entriesRng)) / entriesRng

    # Setting up the polynomial moment matrix M
    M = np.zeros((N,N), dtype=np.float64)
    for col in range(N):
        M[:,col] = entries[col:col+N]

    # Solve the linear equation system from the method of moments
    # on a polynomial to get polynomial coefficients
    W = np.linalg.solve(M, moments2)

    # Return polynomial coefficients
    return(W)
