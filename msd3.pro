;+
; NAME:
;    msd3
;
; PURPOSE:
;    Computes time evolution of the mean-square displacement 
;    of a time series, including estimates of the variance.
;
; CATEGORY:
;    Time series analysis, dynamics.
;
; CALLING SEQUENCE:
;    m = msd3(x)
;
; INPUTS:
;    x: [npts] sequence of measured values assumed to be
;       obtained at equal time intervals.
;
; KEYWORD PARAMETERS:
;    nlags: Maximum lag time to compute msd.
;        Default: NPTS/2
;
; KEYWORD FLAGS:
;    dedrift: Remove linear trend from data.
;        Default: Do not remove drift.
;
; OUTPUTS:
;    m: [2,nlags] array
;        m[0,*]: mean-square displacement
;        m[1,*]: statistical error in the msd.
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 2/16/2008
; 06/10/2010 DGG Initial public release.
; 08/07/2011 DGG and Bhaskar Krishnatreya (NYU): added DEDRIFT
;    keyword to make drift removal optional.  Improved handling
;    of statistical error at the longest lag times.
; 07/26/2012 DGG Improved parameter checking.  Improved dedrift code.
;    Updated syntax.
; 07/26/2012 DGG and Michael Landy, New York University.
;    Simpler calculation of ntau.  Correct normalization of covariance
;    matrix.
; 07/27/2012 DGG fix dedrift code.
;    DGG and ML: Simplify estimate of the variance.
; 02/13/2016 DGG fix for limiting case of only one sample.
;
; Copyright (c) 2008-2016, David G. Grier, Bhaskar Krishnatreya
;    and Michael Landy
;-

function msd3, p, nlags = nlags,  $
               dedrift = dedrift, $
               nooverlap = nooverlap

  COMPILE_OPT IDL2

  umsg = 'USAGE: result = msd3(data)'

  if ~isa(p, /number, /array) then begin
     message, umsg, /inf
     return, -1
  endif

  if size(p, /n_dimensions) ne 1 then begin
     message, umsg, /inf
     message, 'DATA must be a one-dimensional array', /inf
     return, -1
  endif

  npts = n_elements(p)

  if keyword_set(dedrift) then begin
     t = findgen(npts)
     f = poly_fit(t, p, 1)
     p -= f[1] * t
  endif

  if n_elements(nlags) ne 1 then $
     nlags = (npts - 1L)/2L

  m = fltarr(2, nlags + 1)      ; storage for the answer
  
  for tau = 1., nlags do begin  ; Calculate MSD for each lag time
     d = p[tau:*] - p           ; D: all displacements at lag tau
     nsets = tau                ; Number of overlapping data sets at lag tau
     npts = floor(n_elements(d)/nsets) ; Number of samples
                                ; in each uncorrelated set
     d = d[0:npts*nsets-1]      ; Trim D into sets of equal length
     d = reform(d, nsets, npts) ; Reorganize D into uncorrelated sets
     
     if keyword_set(nooverlap) && nsets gt 1 then begin
        d = d[0,*]
        nsets = 1
     endif

     d *= d                     ; D: squared displacements at lag tau in set j
     mu = (npts gt 1) ? $
          total(d, 2)/npts : $  ; Mean-squared displacement at lag tau in set j
          d
     m[0, tau] = total(mu)/nsets ; Overall MSD at lag tau
     
     ;; Compute statistical variance of the squared displacements from
     ;; the MSD at lag time tau.
     ;; This is an estimate for the uncertainty in the MSD.
     if npts gt 1 then begin
        for j = 0, nsets-1 do $ 
           d[j, *] -= mu[j]     ; d^j_i: difference between the
                                ; i-th squared displacement in set j
                                ; and the MSD for set j
        sigma = d # transpose(d) ; Covariance matrix:
                                ; sigma_jk = sum_i d^j_i d^k_i
        sigma /= npts           ; Normalization for Gaussian statistics
        
        sd = total(d^2, 2)/npts ; Mean-squared deviation in each set
        sd = sqrt(sd # transpose(sd)) ; Normalization so that
                                ; sigma[j,k]/sd[j,k] = 1 when sets j
                                ; and k are perfectly correlated
                                ; (diagonal elements, e.g.) and 
                                ; sigma[j,k]/sd[j,k] = 0
                                ; for statistically independent sets.
        devsq = sigma^2 / sd / npts ; Products of deviations from the mean
                                ; weighted by the covariances between the
                                ; sets and normalized by the number of
                                ; data-pair products in each set.
                                ; If the sets were uncorrelated,
                                ; only the diagonal terms would survive.
        m[1, tau] = sqrt(total(devsq) / nsets^2) ; RMS variance of the MSD
     endif
  endfor

  return, m
end

;;; EXAMPLE ;;;
; Run this with .r msd3
; settings
seed = systime(/seconds)
npts = 10000L
nlags = 100L
dedrift = 1
nooverlap = 0
dx = randomn(seed, npts)        ; generate data
x = total(dx, /cumulative)

m = msd3(x, nlags = nlags, dedrift = dedrift, nooverlap = nooverlap)

t = findgen(nlags+1)
p = errorplot(t, m[0,*], m[1,*], symbol = 'o', linestyle = '', $
              xtitle = 'lag', ytitle = 'msd(lag)')
f = poly_fit(t[1:*], m[0,1:*], 2, measure_errors = m[1,1:*], sigma = s)
po = plot(t, f[0] + f[1]*t + f[2]*t^2, color = 'red', /over)
print, transpose(f)
print, s

end
