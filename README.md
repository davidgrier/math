# IDL math

**IDL routines for computing special functions, analyzing
signals and performing mathematical operations**

IDL is the Interactive Data Language, and is a product of
[Exelis Visual Information Solutions](http://www.exelisvis.com).

These routines are licensed under the GPLv3.

## What it does

1. **Averaging and Smoothing**
    * **aziavg**: Averages a two-dimensional data set over angles.
    * **azistd**: Computes the standard deviation relative to the mean of
a two-dimensional data set, averaged over angles, as a function of
radial distance from a specified center.
    * **azimedian**: Computes the median over angles, as a function of
radial distance from a specified center.
    * **msd3**: Computes the mean-square deviation of a time series.
    * **savgol2D**: Generates two-dimensional Savitzky-Golay smoothing
and derivative kernels.
    * **zerocrossings**: Computes the coordinates of the zero crossings
of a function.
    * **extrema**: Computes the coordinates and values of the local
extrema of a function.
2. **Statistics**
    * **iqr**: Computes the inter-quartile range of a one-dimensional data set.
    * **sigmatrim**: Computes the average and standard deviation of a data
set while ignoring contributions from outlying data points.
    * **mad**: Estimate signal-to-noise ratio in a signal or image with
the median absolute deviation from the median.
3. **Kernel Density Estimation**
    * **kde**: Estimates the N-dimensional probability density underlying a set of
discrete samples (measurements) using the kernel density estimator method.
    * **akde**: Estimates the N-dimensional probability density underlying a set of
discrete samples (measurements) using the _adaptive_ kernel density estimator method.
4. **Special Functions**
    * **atanh**: Hyperbolic arctangent.
    * **bessjp**: First derivative of Bessel functions of the first kind.
    * **besselj**: Bessel function of the the first kind of order 0.
    * **sphericalbesselj**: Spherical Bessel function of the first kind.
    * **zernike**: Class of functions related to Zernike polynomials.
	* **zerojp**: Zeros of the derivatives of Bessel functions of the
      first kind.
5. **Quaternions**
    * **quaternion__define**: Defines an object that implements quaternion methods.

