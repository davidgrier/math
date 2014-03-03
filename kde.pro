;+
; NAME:
;    KDE
;
; PURPOSE:
;    Estimate the probability density underlying a set of discrete
;    samples (measurements) using the kernel density estimator method.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    rho = kde(p, x)
;
; INPUTS:
;    p: discrete measurements of the desired distribution: [ndimensions, nmeasurements]
;
;    x: values at which to estimate the probability density: [ndimensions, nsamples]
;
; KEYWORD PARAMETERS:
;    weight: weighting for sampled points.
;
; KEYWORD OUTPUTS:
;    variance: statistical variance of the result.
;
;    bias: asymptotic estimate for the statistical bias at each point.
;         The absolute value is an estimate for the error in rho.
;         NOTE: Not implemented for triangular kernel.
;
;    mse: asymptotic mean square error at each point.
;         NOTE: Not implemented for triangular kernel.
;
;    scale: smoothing factor, also called the bandwidth, used to
;        compute the density estimate
;
; KEYWORD FLAGS:
;    By default, KDE uses the Epanechnikov kernel to compute
;    the kernel density estimate, this can be overridden by
;    setting one of the following flags:
;    GAUSSIAN:   Gaussian kernel
;    TRIANGULAR: triangular kernel
;    BIWEIGHT:   biweight kernel
;
; OUTPUTS:
;    rho: probability density estimated at each value specified by x
;
; RESTRICTIONS:
;    Gaussian kernel used for multivariate systems.
;
; PROCEDURE:
;    rho_i = (1/n) \sum_{j = 1}^n K((p_j - x_i)/h) / h
;
;    where h is the estimated optimal smoothing parameter and
;    where K(z) is the selected kernel.
;
; REFERENCES:
; 1. B. W. Silverman,
;    Density Estimation for Statistics and Data Analysis
;    (CRC Press, Boca Raton, 1998)
;
; 2. B. E. Hansen, "Lecture Notes on Nonparametrics"
;    University of Wisconsin, Spring 2009
;    http://www.ssc.wisc.edu/~bhansen/718/NonParametrics1.pdf
;
; 3. Z. Ouyang, "Univariate Kernel Density Estimation"
;    Duke University, August 2005
;    http://www.stat.duke.edu/~zo2/shared/research/readings/kernelsmoothing.pdf
;
; EXAMPLE:
;    IDL> p = randomn(seed, 10000)
;    IDL> x = 2. * findgen(100)/99.
;    IDL> rho = kde(p,x)
;    IDL> plot, x, rho
;    IDL> plot, x, histogram(p, min=0, max=2, nbins=100), /noerase
;
; MODIFICATION HISTORY:
; 09/18/2010 Written by David G. Grier, New York University
; 10/08/2011 DGG Corrected for case when more than 3/4 of input data
;    have same value.  Thanks to Dan Hartung (UW Madison) 
;    for bringing this bug to light.
; 10/18/2011 DGG Added WEIGHT keyword.
; 12/09/2012 DGG Use rebin() rather than # operator for efficiency
;    and clarity in n-dimensional code.  Test scale with
;    arg_present().  Move normalization out of loops for efficiency.
;    Corrected n-dimensional normalization.  Updated usage messages.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 02/10/2014 DGG Added VARIANCE keyword.
; 02/13/2014 DGG Cast indexes to long to avoid integer overruns.
;    Cast nx to float.
; 02/25/2014 DGG Implemented MSE.
; 03/01/2014 DGG Revised MSE calculations.
; 03/03/2014 DGG Implemented BIAS.
;
; Copyright (c) 2010-2014 David G. Grier
;-

function kde_nd, x, y, $
                 weight = weight, $
                 scale = scale, $
                 variance = variance, $
                 bias = bias, $
                 mse = mse

COMPILE_OPT IDL2, HIDDEN

sx = size(x, /dimensions)
sy = size(y, /dimensions)

nd = sx[0]                      ; number of dimensions
nx = float(sx[1])               ; number of data points
ny = long(sy[1])                ; number of sampling points

if n_elements(weight) ne nx then $
   weight = replicate(1., nx)

; optimal smoothing parameter in each dimension
; Silverman Eqs. (3.30) and (3.31)
sx = stddev(x, dimension=2)
rx = fltarr(nd)
for d = 0, nd-1 do $
   rx[d] = iqr(x[d,*])
h = sx
w = where(rx gt 1e-10, ngood)
if ngood gt 0 then $
   h[w] = h[w] < rx[w]/1.34
h *= 0.9 / nx^0.2

if arg_present(scale) then scale =  h

; density estimate
; Silverman Eq. (2.15) and Table 3.1
res = fltarr(ny)
variance = fltarr(ny)
bias = fltarr(ny)
mse = fltarr(ny)
hfac = rebin(h, nd, nx, /sample)

norm = 1./(2. * !pi * total(h^2))^(nd/2.) / nx
for j = 0L, ny-1L do begin
   z = 0.5 * total(((x - rebin(y[*,j], nd, nx, /sample)) / hfac)^2, 1)
   w = where(z lt 20., ngood)
   if ngood gt 0L then begin
      ker = norm * exp(-z[w])
      val = weight[w] * ker
      res[j] = total(val)
      variance[j] = total((val - res[j])^2)/nx
      bias[j] = total(weight[w]*(1. - z[w])*ker)/2.
      mse[j] = (norm / 2.^(nd/2.)) * res[j]^2/total(ker) + bias[j]^2
   endif
endfor

return, res
end

function kde_1d, x, y, $
                 weight = weight, $
                 scale = scale, $
                 biweight = biweight, $
                 triangular = triangular, $
                 gaussian = gaussian, $
                 variance = variance, $
                 bias = bias, $
                 mse = mse

COMPILE_OPT IDL2, HIDDEN

nx = float(n_elements(x))       ; number of data points
ny = n_elements(y)              ; number of samples

if n_elements(weight) ne nx then $
   weight = replicate(1., nx)

; optimal smoothing parameter
; Silverman Eqs. (3.30) and (3.31)
sx = stddev(x)                  ; standard deviation
rx = iqr(x)                     ; interquartile range
if rx lt 1e-10 then $           ; more than 3/4 data have same value
   h = 0.9 * sx / nx^0.2 $
else $
   h = 0.9 * (sx < rx/1.34) / nx^0.2

if arg_present(scale) then scale = h

; density estimate
; Silverman Eq. (2.15) and Table 3.1
t = x/h
s = y/h
res = fltarr(ny)                ; result
variance = fltarr(ny)           ; variance in result
bias = fltarr(ny)
mse = fltarr(ny)                ; asymptotic mean-squared error

if keyword_set(biweight) then begin
   norm = (15./16.) / (h * nx)
   for j = 0L, ny-1L do begin
      z = (t - s[j])^2
      w = where(z lt 1., ngood)
      if ngood gt 0L then begin
         ker = norm * (1. - z[w])^2
         val = weight[w] * ker
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx
         bias[j] = 10./7. * norm * total(weight[w] * (3.*z[w] - 1.))
         mse[j] = res[j]^2/total(ker) / (7. * nx * h) + bias[j]^2
      endif
   endfor
endif $                     
else if keyword_set(triangular) then begin
   norm = 1./(h * nx)
   for j = 0L, ny-1L do begin
      z = abs(t - s[j])
      w = where(z lt 1., ngood)
      if ngood gt 0L then begin
         ker = norm * (1. - z[w])
         val = weight[w] * ker[w]
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx
      endif
   endfor
   res *= 1./(h*nx)
endif $                     
else if keyword_set(gaussian) then begin
   norm = 1./(sqrt(2.*!pi) * h * nx)
   for j = 0L, ny-1L do begin
      z = (t - s[j])^2/2.
      w = where(z lt 20., ngood)
      if ngood gt 0L then begin
         ker = norm*exp(-z[w])
         val = weight[w] * ker
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx
         bias[j] = total(weight[w] * (1. - z[w]) * ker)/2.
         mse[j] = (norm / sqrt(2.)) * res[j]^2 / total(ker) + bias[j]^2                  
      endif
   endfor
endif $
else begin                      ; Epanechnikov
   norm = (3./4./sqrt(5.)) / (h * nx)
   for j = 0L, ny-1L do begin
      z = (t - s[j])^2/5.
      w = where(z lt 1., ngood)
      if ngood gt 0L then begin
         ker = norm * (1. - z[w])
         val = weight[w] * ker
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx
         bias[j] = norm/5. * total(weight[w])
         mse[j] = norm * res[j]^2 / total(ker) + bias[j]^2
      endif
   endfor
endelse

return, res
end

function kde, x, y, $
              weight = weight, $
              scale = scale, $
              gaussian = gaussian, $
              biweight = biweight, $
              triangular = triangular, $
              variance = variance, $
              bias = bias, $
              mse = mse

COMPILE_OPT IDL2

umsg = 'rho = kde(p, x, [weight = w])'

if n_params() ne 2 then $
   message, umsg
  
sx = size(x)
sy = size(y)

if sx[0] gt 2 then begin
   message, umsg, /inf
   message, 'P must be organized as [ndimensions, npoints]'
endif

if sy[0] ne sx[0] then begin
   message, umsg, /inf
   message, 'P and X must have the same numbers of dimensions'
endif

if (sx[0] eq 2) and (sx[1] ne sy[1]) then begin
   message, usg, /inf
   message, 'P and X must have the same numbers of dimensions'
endif

ndims = (sx[0] eq 2) ? sx[1] : 1

if ndims gt 1 then begin
   if keyword_set(biweight) or keyword_set(triangular) then $
      message, 'Multidimensional: using Gaussian kernel', /inf
   res = kde_nd(x, y, weight = weight, scale = scale, $
                variance = variance, bias = bias, mse = mse)
endif else $
   res = kde_1d(x, y, weight = weight, scale = scale, $
                gaussian = gaussian, $
                biweight = biweight, $
                triangular = triangular, $
                variance = variance, $
                bias = bias, $
                mse = mse)

return, res
end
