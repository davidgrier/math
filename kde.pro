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
;    scale: smoothing factor, also called the bandwidth, used to
;        compute the density estimate
;
; KEYWORD FLAGS:
;    By default, KDE uses the Epanechnikov kernel to compute
;    the kernel density estimate, this can be overridden by
;    setting one of the following flags:
;    GAUSSIAN: use Gaussian kernel
;    TRIANGULAR: use triangular kernel
;    BIWEIGHT: use biweight kernel
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
; REFERENCE:
; B. W. Silverman,
; Density Estimation for Statistics and Data Analysis
; (CRC Press, Boca Raton, 1998)
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
;
; Copyright (c) 2010-2013 David G. Grier
;-

function kde_nd, x, y, $
                 weight = weight, $
                 scale = scale

COMPILE_OPT IDL2, HIDDEN

sx = size(x, /dimensions)
sy = size(y, /dimensions)

nd = sx[0]                      ; number of dimensions
nx = sx[1]                      ; number of data points
ny = sy[1]                      ; number of sampling points

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
hfac = rebin(h, nd, nx, /sample)

for j = 0, ny-1 do begin
   z = 0.5 * total(((x - rebin(y[*,j], nd, nx, /sample)) / hfac)^2, 1)
   w = where(z lt 20, ngood)
   if ngood gt 0 then $
      res[j] = total(weight[w]*exp(-z[w]))
endfor

res *= 1. / ((2. * !pi) * total(h^2))^(nd/2.) / nx ; normalization

return, res
end

function kde_1d, x, y, $
                 weight = weight, $
                 scale = scale, $
                 biweight = biweight, $
                 triangular = triangular, $
                 gaussian = gaussian

COMPILE_OPT IDL2, HIDDEN

nx = n_elements(x)              ; number of data points
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

if keyword_set(biweight) then begin
   for j = 0, ny-1 do begin
      z = (t - s[j])^2
      w = where(z lt 1., ngood)
      if ngood gt 0. then $
         res[j] = total(weight[w]*(1.-z[w])^2)
   endfor
   res *= 15./(16.*h*nx)
endif $                     
else if keyword_set(triangular) then begin
   for j = 0, ny-1 do begin
      z = abs(t - s[j])
      w = where(z lt 1., ngood)
      if ngood gt 0 then $
         res[j] = total(weight[w]*(1. - z[w]))
   endfor
   res *= 1./(h*nx)
endif $                     
else if keyword_set(gaussian) then begin
   for j = 0, ny-1 do begin
      z = 0.5 * (t - s[j])^2
      w = where(z lt 20, ngood)
      if ngood gt 0 then $
         res[j] = total(weight[w]*exp(-z[w]))
   endfor
   res *= 1./(sqrt(2.*!pi)*h*nx)
endif $
else begin                      ; Epanechnikov
   for j = 0, ny-1 do begin
      z = (t - s[j])^2
      w = where(z lt 5, ngood)
      if ngood gt 0 then $
         res[j] = total(weight[w]*(1.-z[w]/5.))
   endfor
   res *= 0.75/(sqrt(5.)*h*nx)
endelse

return, res
end

function kde, x, y, $
              weight = weight, $
              scale = scale, $
              gaussian = gaussian, $
              biweight = biweight, $
              triangular = triangular

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
   res = kde_nd(x, y, weight = weight, scale = scale)
endif else $
   res = kde_1d(x, y, weight = weight, scale = scale, $
               gaussian = gaussian, $
               biweight = biweight, $
               triangular = triangular)

return, res
end


