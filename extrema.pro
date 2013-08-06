;+
; NAME:
;    extrema
;
; PURPOSE:
;    Returns coordinates of extrema in a signal
;
; CATEGORY:
;    Signal processing
;
; CALLING SEQUENCE:
;    xm = extrema([x], s)
;
; INPUTS:
;    s: one-dimensional signal
;
; OPTIONAL INPUTS:
;    x: coordinates of signal samples: s = s(x)
;
; OUTPUTS:
;    xm: coordinates of extrema
;
; KEYWORD OUTPUTS:
;    values: value of signal at extrema: s(xm)
;
;    ismax: 1 if extremum is maximum, 0 otherwise
;    ismin: 0 if extremum is minimum, 1 otherwise
;
; KEYWORD FLAGS:
;    maxima: If set, only return maxima 
;    minima: If set, only return minima
;        Default: return all extrema
;
; PROCEDURE:
;    Find zero crossings of derivative of signal.
;
; EXAMPLE:
;    IDL> print, extrema(sin(findgen(100))) / !pi
;
; MODIFICATION HISTORY:
; 05/14/2013 Written by David G. Grier, New York University
; 05/15/2013 DGG Added ISMAX keyword.
; 05/18/2013 DGG Added ISMIN keyword
;
; Copyright (c) 2013 David G. Grier
;-
function extrema, arg1, arg2, $
                  maxima = maxima, $
                  minima = minima, $
                  values = values, $
                  ismax = ismax, $
                  ismin = ismin, $
                  count = count

COMPILE_OPT IDL2

s = (n_params() eq 1) ? arg1 : arg2
ds = deriv(s)

nc = zerocrossings(ds, falling = ismax, count = count)
ismin = ~ismax

if count le 0 then $
   return, []

if keyword_set(maxima) then $
   nc = nc[where(ismax, count, /NULL)] $
else if keyword_set(minima) then $
   nc = nc[where(ismin, count, /NULL)]

if count le 0 then $
   return, []

n = floor(nc)

if arg_present(values) then begin
   if keyword_set(maxima) then $
      values = s[n] > s[n+1] $
   else if keyword_set(minima) then $
      values = s[n] < s[n+1] $
   else $
      values = (s[n] > s[n+1]) * (1. - rising) + (s[n] < s[n-1]) * rising
endif

if n_params() eq 2 then begin
   fh = nc - n
   return, arg1[n] * (1. - fh) + arg1[n+1] * fh
endif

return, nc
end
