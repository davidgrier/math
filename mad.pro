;+
; NAME:
;    mad
;
; PURPOSE:
;    Computes median absolute difference (MAD) estimate for the noise
;    in an image.
;
; CATEGORY:
;    Image processing
;
; CALLING SEQUENCE:
;    noise = mad(a)
;
; INPUTS:
;    a: image data
;
; KEYWORD PARAMETERS:
;    w: width of median window.  Default: 10 [pixels]
;
; KEYWORD FLAGS:
;    fast: if set, use smooth() rather than median() to suppress noise.
;
; OUTPUTS:
;    noise: estimate for the amplitude of the Gaussian random noise
;           in the image.
;
; MODIFICATION HISTORY:
; 01/15/2013 Written by David G. Grier, New York University
; 03/22/2013 DGG Added FAST keyword flag.
;
; Copyright (c) 2013 David G. Grier
;-
function mad, a, w = w, fast = fast

COMPILE_OPT IDL2

umsg = 'noise = mad(a)'

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   return, -1.
endif

if ~isa(w, /number, /scalar) then $
   w = 10

return, keyword_set(fast) ? median(abs(float(a) - smooth(a, w))) : median(abs(float(a) - median(a, w)))
end
