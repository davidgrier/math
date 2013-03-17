# IDL math

**IDL routines for performing calculating functions and
performing mathematical operations**

IDL is the Interactive Data Language, and is a product of
[Exelis Visual Information Solutions](http://www.exelisvis.com).

These routines are licensed under the GPLv3.

## What it does

1. **Averaging and Smoothing**

* **aziavg**: Averages a two-dimensional data set over angles.

* **azistd**: Computes the standard deviation relative to the mean of
a two-dimensional data set, averaged over angles, as a function of
radial distance from the center.

* **savgol2D**: Generates two-dimensional Savitzky-Golay smoothing
and derivative kernels.

2. **Statistics**

* **iqr**: Computes the inter-quartile range of a one-dimensional data set.

* **sigmatrim**: Computes the average and standard deviation of a data
set while ignoring contributions from outlying data points.

3. **Kernel Density Estimation**
* **kde**: Estimates the N-dimensional probability density underlying a set of
discrete samples (measurements) using the kernel density estimator method.

* **akde**: Estimates the N-dimensional probability density underlying a set of
discrete samples (measurements) using the _adaptive_ kernel density estimator method.

4. **Special Functions**

* **atanh**: Hyperbolic arctangent.

* **besselj**: Bessel function of the the first kind of order 0.

* **sphericalbesselj**: Spherical Bessel function of the first kind.

5. **Quaternions**

* **quaternion__define**: Defines an object that implements quaternion methods.

