;+
; NAME:
;    zerocrossings
;
; PURPOSE:
;    Compute zero crossings of a function
;
; CATEGORY:
;    Signal analysis
;
; CALLING SEQUENCE:
;    xc = zerocrossings([x], s)
;
; INPUTS:
;    s: one-dimensional signal
;
; OPTIONAL INPUTS:
;    x: coordinates of signal: s = s(x)
;
; KEYWORD FLAGS:
;    rising: If set, return coordinates of rising zero crossings
;    falling: If set, return coordinates of falling zero crossings
;        Default: Return all zero crossings.
;
; OUTPUTS:
;    xc: coordinates of zero crossings.  If x is not specified, the
;        signal is assumed to be evenly sampled at unit intervals.
;
; EXAMPLE:
;    IDL> xc = zerocrossings(sin(findgen(100)))
;
; MODIFICATION HISTORY:
; 05/14/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;-

function zerocrossings, arg1, arg2, rising = rising, falling = falling

COMPILE_OPT IDL2

s = (n_params() eq 1) ? arg1 : arg2
a = s gt 0.

w = where(a[1:*] ne a, ncrossings)
if ncrossings le 0 then $
   return, []

if n_params() eq 1 then begin
   x1 = w
   x2 = w+1
endif else begin
   x1 = arg1[w]
   x2 = arg1[w+1]
endelse

y1 = s[w]
y2 = s[w+1]

if keyword_set(rising) then begin
   w = where(y1 lt y2, nrising)
   if nrising le 0 then $
      return, []
   return, (y2[w]*x1[w] - y1[w]*x2[w])/(y2[w] - y1[w])
endif

if keyword_set(falling) then begin
   w = where(y1 gt y2, nfalling)
   if nfalling le 0 then $
      return, []
   return, (y1[w]*x2[w] - y2[w]*x1[w])/(y1[w] - y2[w])
endif

y1 = abs(temporary(y1))
y2 = abs(temporary(y2))
return, (y1*x2 + y2*x1)/(y2 + y1)
end
