;+
; NAME:
;    AKDE
;
; PURPOSE:
;    Estimate the probability density underlying a set of discrete
;    samples (measurements) using the adaptive kernel density estimator method.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    d = akde(x, s)
;
; INPUTS:
;    x: discrete samples of the desired distribution.
;
;    s: values for which the probability density is required.
;
; KEYWORD PARAMETERS:
;    alpha: adaptability exponent in range [0,1].
;        Default: alpha = 0.5
;        alpha = 0 is equivalent to non-adaptive KDE
;        alpha = 1 approaches nearest-neighbor adaptive weighting.
;
;    weight: weighting for each sample point.
;
; KEYWORD OUTPUTS:
;    variance: statistical variance of the result
;
; KEYWORD FLAGS:
;    By default, AKDE uses the Epanechnikov kernel to compute
;    the kernel density estimate for one-dimensional data, 
;    this can be overridden by setting one of the following flags:
;    GAUSSIAN: use Gaussian kernel
;    TRIANGULAR: use triangular kernel
;    BIWEIGHT: use biweight kernel
;
;    For higher-dimensional data, only the Gaussian kernel is
;    implemented and these flags are ignored.
;
; OUTPUTS:
;    d: probability density estimated at each value specified by s.
;
; RESTRICTIONS:
;    Multivariate estimates are only computed with Gaussian kernels.
;
; PROCEDURE:
;    d_i = (1/N) \sum_{j = 1}^N w_j K((x_j - s_i)/h_j)
;
;    where K(x) is the selected kernel,
;    where h_j is the estimated optimal smoothing parameter at
;    position x_j, and w_j is the weight to be assigned to the
;    measurement at position x_j.
;
; REFERENCE:
; B. W. Silverman,
; Density Estimation for Statistics and Data Analysis
; (CRC Press, Boca Raton, 1998)
;
; EXAMPLE:
;    IDL> x = randomn(seed, 1000)
;    IDL> t = 2. * findgen(100)/99.
;    IDL> d = akde(x,t)
;    IDL> plot, t, d
;    IDL> plot, t, histogram(x, min=0, max=2, nbins=100), /noerase
;
; MODIFICATION HISTORY:
; 09/18/2010 Written by David G. Grier, New York University
; 09/26/2010 DGG Added COMPILE_OPT.  First version of AKDE_ND
; 10/08/2011 DGG Corrected checks on results of where()
;    for sparsely sampled data.  Thanks to Dan Hartung (UW Madison) for
;    bringing this bug to light.
;    Added ALPHA keyword.
; 10/27/2011 DGG Major overhaul of n-dimensional code to ensure
;    consistency with Silverman Chapter 5.
; 12/09/2012 DGG replace # with rebin() in n-dimensional code for
;    efficiency and clarity.  Corrected n-dimensional normalization.
;    Upgraded usage messages.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 02/10/2014 DGG Added VARIANCE keyword.
; 02/13/2014 DGG Cast indexes to long to avoid integer overflow
;    errors. Cast nx to float.
;
; Copyright (c) 2010-2014 David G. Grier
;-

function akde_nd, x, y, $
                  weight = weight, $
                  alpha = alpha, $
                  variance = variance

COMPILE_OPT IDL2, HIDDEN

sx = size(x, /dimensions)
sy = size(y, /dimensions)

nd = sx[0]                      ; number of dimensions
nx = float(sx[1])               ; number of data points
ny = long(sy[1])                ; number of sampling points

if n_elements(weight) ne nx then $
   weight = replicate(1., nx)

; Method described by Silverman Sec. 5.3.1
; 1. pilot estimate of the density at the data points
f = kde(x, x, scale = h)

; 2. local bandwidth factor
if n_elements(alpha) ne 1 then $
   alpha = 0.5                  ; sensitivity parameter: 0 <= alpha <= 1
g = exp(mean(alog(f)))          ; geometric mean density
lambda = (g/f)^alpha            ; Eq. (5.7): factor for each data point

; 3. adaptive density estimate
res = fltarr(ny)
variance = fltarr(ny)
hfac = h # lambda               ; smoothing factor for each point

norm = 1./((2.*!pi) * total(hfac^2, 1))^(nd/2.)/nx ; normalization

for j = 0L, ny - 1L do begin
   z = 0.5 * total(((x - rebin(y[*,j], nd, nx, /sample))/hfac)^2 , 1)
   w = where(z lt 20., ngood)
   if ngood gt 0L then begin
      ker = norm[w] * exp(-z[w])
      val = weight[w] * ker
      res[j] = total(val)
      variance[j] = total((val - res[j])^2)/nx^2
   endif
endfor

return, res
end

function akde_1d, x, y, $
                  weight = weight, $
                  biweight = biweight, $
                  triangular = triangular, $
                  gaussian = gaussian, $
                  alpha = alpha, $
                  variance = variance

COMPILE_OPT IDL2, HIDDEN

nx = float(n_elements(x))       ; number of data points
ny = n_elements(y)              ; number of samples

if n_elements(weight) ne nx then $
   weight = replicate(1., nx)

; Method described by Silverman Sec. 5.3.1
; 1. pilot estimate of the density at the data points
f = kde(x, x, scale = h)        ; use Epanechnikov (p. 102)

; 2. local bandwidth factor
if n_elements(alpha) ne 1 then $
   alpha = 0.5                  ; sensitivity parameter: 0 <= alpha <= 1
g = exp(mean(alog(f)))          ; geometric mean density
lambda = (g/f)^alpha            ; Eq. (5.7)

; 3. adaptive density estimate
t = x / h
s = y / h
res = findgen(ny)
variance = fltarr(ny)

if keyword_set(biweight) then begin
   norm = (15./16.) / (h * nx)
   for j = 0L, ny-1L do begin
      z = ((t - s[j])/lambda)^2
      w = where(z lt 1., ngood)
      if ngood gt 0L then begin
         ker = norm * (1. - z[w])^2 / lambda[w]
         val = weight[w] * ker
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx^2
      endif
   endfor
endif $                     
else if keyword_set(triangular) then begin
   norm = 1./(h * nx)
   for j = 0L, ny-1L do begin
      z = abs(t - s[j])/lambda
      w = where(z lt 1., ngood)
      if ngood gt 0L then begin
         ker = norm * (1. - z[w]) / lambda[w]
         val = weight[w] * ker[w]
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx^2
      endif
   endfor
endif $                     
else if keyword_set(gaussian) then begin
   norm = 1./(sqrt(2.*!pi) * h * nx)
   for j = 0L, ny-1L do begin
      z = 0.5 * ((t - s[j])/lambda)^2
      w = where(z lt 30., ngood)
      if ngood gt 0L then begin
         ker = norm * exp(-z[w]) / lambda[w]
         val = weight[w] * ker
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx^2
      endif
   endfor
endif $
else begin                      ; Epanechnikov
   norm = 0.75 / (sqrt(5.) * h * nx)
   for j = 0L, ny-1L do begin
      z = ((t - s[j])/lambda)^2
      w = where(z lt 5., ngood)
      if ngood gt 0L then begin
         ker = norm * (1. - z[w]/5.) / lambda[w]
         val = weight[w] * ker
         res[j] = total(val)
         variance[j] = total((val - res[j])^2)/nx^2
      endif
   endfor
endelse

return, res
end


function akde, x, y, $
               weight = weight, $
               gaussian = gaussian, $
               biweight = biweight, $
               triangular = triangular, $
               alpha = alpha, $
               variance = variance

COMPILE_OPT IDL2

sx = size(x)
sy = size(y)
umsg = 'rho = akde(p, x, [weight = w], [alpha = a])'

if n_params() ne 2 then $
   message, umsg

if sx[0] gt 2 then begin
   message, umsg, /inf
   message, 'P must be organized as [ndimensions, nmeasurements]'
endif

if sy[0] ne sx[0] then begin
   message, umsg, /inf
   message, 'P and X must have the same number of dimensions'
endif

if (sx[0] eq 2) and (sx[1] ne sy[1]) then begin
   message, umsg, /inf
   message, 'P and X must have the same number of dimensions'
endif

ndims = (sx[0] eq 2) ? sx[1] : 1

if n_elements(alpha) eq 1 then begin
   if (alpha lt 0) or (alpha gt 1) then begin
      message, umsg, /inf
      message, '0 < alpha < 1', /inf
      message, "setting alpha = 0.5", /inf
      alpha = 0.5
   endif
endif

if ndims gt 1 then $
   return, akde_nd(x, y, weight = weight, alpha = alpha, variance = variance) $
else $
   return, akde_1d(x, y, $
                   weight = weight, $
                   gaussian = gaussian, $
                   biweight = biweight, $
                   triangular = triangular, $
                   alpha = alpha, $
                   variance = variance)
end
