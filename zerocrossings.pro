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
; OUTPUTS:
;    xc: coordinates of zero crossings.  If x is not specified, the
;        signal is assumed to be evenly sampled at unit intervals.
;
; KEYWORD OUTPUTS:
;    rising: 1 for rising zero-crossings, 0 for falling.
;    falling: opposite of rising.
;
; EXAMPLE:
;    IDL> xc = zerocrossings(sin(findgen(100)))
;
; MODIFICATION HISTORY:
; 05/14/2013 Written by David G. Grier, New York University
;    Use RISING as a set of flags on output.
; 05/15/2013 DGG Added FALLING keyword.
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

if arg_present(rising) then $
   rising = y2 gt y1

if arg_present(falling) then $
   falling = y2 lt y1

y1 = abs(temporary(y1))
y2 = abs(temporary(y2))
return, (y1*x2 + y2*x1)/(y2 + y1)
end
