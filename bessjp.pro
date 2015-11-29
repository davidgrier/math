; docformat = 'rst'

;+
; Calculates the first derivative of the first kind of Bessel function
; of order $n$ for real arguments: $J_n'(x)$
;
; :Params:
;    n : in, required, type=integer
;        order of the Bessel function
;    x : in, required, type=fltarr
;        argument of the Bessel function
;
; :Returns:
;    Double-valued array of the same dimensions as x
;    containing the value of the derivative of the Bessel function.
;
; :Author:
;    David G. Grier, New York University
;
; :History:
;    11/29/2015 Written by DGG
;
; :Copyright:
;    Copyright (c) 2015 David G. Grier
;-
function bessjp, n, x

  COMPILE_OPT IDL2

  if n eq 0 then $
     return, -beselj(x, 1, /double)

  val = 0.d * x
  w = where(abs(x) gt 1d-8, nlarge)
  if nlarge gt 0 then $
     val[w] = beselj(x[w], n-1, /double) - double(n)/x[w] * beselj(x[w], n, /double)

  return, (n_elements(val) eq 1) ? val[0] : val
end

  
