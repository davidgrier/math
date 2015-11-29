; docformat= 'rst'

;+
; Calculates the $k$-th zero of the derivative of the
; first kind of Bessel function of order $n$: $j'_{n,k}$.
;
; :Params:
;    n : in, required, type=integer
;        order of the Bessel function
;    k : in, required, type=integer
;        index of the requested zero, starting from $k = 1$
;
; :Returns:
;    Double-valued zero, $j'_{n,k}$.
;
; :Author:
;    David G. Grier, New York University
;
; :History:
;    11/29/2015 Written by DGG
;
; :Reference:
;    Abramowitz and Stegun
;    http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tzerojp_f90.txt
;
; :Copyright:
;    Copyright (c) 2015 David G. Grier
;-
function zerojp_usage, msg

  COMPILE_OPT IDL2, HIDDEN

  return, 'USAGE: zerojp(n,k): '+msg
end

function zerojp, n, k

  COMPILE_OPT IDL2

  fn = abs(double(n))
  fk = double(k)

  if k lt 1 then $
     message, zerojp_usage('k must be a positive integer')
  
  if k eq 1 then begin
     if n eq 0 then $
        return, 0.D $
     else begin                 ; Tchebychev's series for k <= n
        f1 = fn^(1.d/3.d)
        f2 = f1 * f1 * fn
        val = fn + 0.8086165D * f1 + 0.72490D/f1 - $
              0.05097D/fn + 0.0094D/f2
     endelse
  endif else begin              ; MacMahon's series for k >> n
     b0 = (fk + 0.5d*fn - 0.75d) * !dpi
     b1 = 8.D * b0
     b2 = b1 * b1
     b3 = 3.D * b1 * b2
     b5 = 5.D * b3 * b2
     b7 = 7.D * b5 * b2
     t0 = 4.D * fn * fn
     t1 = t0 + 3.D
     t3 = 4.d * ((7.D * t0 + 82.D) * t0 - 9.D)
     t5 = 32.D * (((83.D * t0 + 2075.D) * t0 - 3039.D) * t0 + 3537.D)
     t7 = 64.D * ((((6949.D * t0 + 296492.D) * t0 - 1248002.D) * t0 + $
                   7414380.D) * t0 - 5853627.D)
     val = b0 - (t1/b1) - (t3/b3) - (t5/b5) - (t7/b7)
  endelse
  
  ;;; improve solution by secant method
  v0 = 0.9D * val
  q0 = bessjp(n, v0)
  repeat begin
     q1 = bessjp(n, val)
     v1 = val - q1 * (val - v0)/(q1 - q0)
     dv = v1 - val
     v0 = val
     q0 = q1
     val = v1
  endrep until abs(dv) lt 1d-8
  
  return, val
end
