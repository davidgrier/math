;+
; NAME:
;    zernike
;
; PURPOSE:
;    Class of functions associated with normalized Zernike polynomials
;
; PUBLIC METHODS:
;    zernike.Rnm(n, m, rho): Radial Zernike polynomial
;
;    zernike.Rj(j, rho): Noll-indexed radial Zernike polynomial
;
;    zernike.Znm(n, m, rho, theta): Zernike polynomial
;
;    zernike.Zj(j, rho, theta): Noll-indexed Zernike polynomial
;
;    zernike.cZnm(n, m, x, y): Cartesian Zernike polynomial
;
;    zernike.cZj(j, x, y): Noll-indexed Cartesian Zernike polynomial
;
;    zernike.dxZnm(n, m, x, y): Derivative of Cartesian Zernike
;        polynomial along x axis
;
;    zernike.dyZnm(n, m, x, y): Derivative of Cartesian Zernike
;        polynomial along y axis
;
;    zernike.dxZj(j, x, y): Derivative of Noll-indexed Cartesian
;        Zernike polynomial along x axis
;
;    zernike.dyZj(j, x, y): Derivative of Noll-indexed Cartesian
;        Zernike polynomial along y axis
;
; ARGUMENTS:
;    n, m: scalar integer indexes of Zernike polynomials.
;        Note: not all combinations of n and m correspond to valid
;        Zernike polynomials.  Invalid combinations yield zero.
;
;    j: scalar integer Noll index.
;       Note: Values of j less than 1 yield zero return.
;
;    rho, theta: polar coordinates at which to evaluate Zernike
;       polynomials.  The radial coordinate rho should be scaled
;       for the unit circle.  Values may be scalars or arrays with
;       the same dimensions.  Returned values will have the same
;       dimensions as these inputs.
;
;    x, y: Cartesian coordinates at which to evaluate Zernike
;       polynomials.  Coordinates should be centered and scaled
;       to fit within the unit circle.  Values may be scalars or
;       arrays with the same dimensions.  Returned values will have
;       the same dimensions as these inputs.
;        
; OTHER METHODS:
;    This class includes other methods that should be considered
;    experimental, and are not intended for production use.
;
; MODIFICATION HISTORY:
; 03/15/2014 Written by David G. Grier, New York University.
;
; Copyright (c) 2014 David G. Grier
;-

;;;;; PRIVATE METHODS

;;;;;
;
; zernike::binomial()
;
; Binomial coefficients
;
function zernike::binomial, n, k

COMPILE_OPT IDL2, STATIC, HIDDEN

result = lngamma(n+1.) - lngamma(k+1.) - lngamma(n-k+1.)

m = max(result, /nan)

if m gt 88. then $
   return, exp(double(result))

if m gt 43. then $
   return, exp(result)

return, floor(0.5 + exp(result), l64 = (m gt 21.))

end

;;;;;
;
; zernike::coefficient(n, m, k)
;
function zernike::coefficient, n, m,  k

COMPILE_OPT IDL2, STATIC, HIDDEN

return, (-1)^k * $
        zernike.binomial(n-k, k) * $
        zernike.binomial(n - 2.*k, (n-m)/2. - k)
end

;;;;;
;
; zernike::noll, j, n, m
;
; Noll sequence coefficients, starting with j = 1.
;
pro zernike::noll, j, n, m

COMPILE_OPT IDL2, STATIC, HIDDEN

n1 = (sqrt(1. + 8.*j) - 1.)/2.
n = floor(n1)
if n1 eq n then $
   n--

k = (n+1.)*(n+2.)/2.
m = n - 2*floor((k-j)/2.)

if (m ne 0) and (j mod 2 eq 0) then $
   m *= -1

end

;;;;; PUBLIC METHODS

;;;;;
;
; zernike::Rnm(n, m, rho)
;
; Radial Zernike polynomial
;
function zernike::Rnm, n, _m, rho

COMPILE_OPT IDL2, STATIC

z = 0.
m = abs(_m)
for k = 0, (n-m)/2 do $
   z += zernike.coefficient(n, m, k) * rho^(n - 2*k)

return, z
end

;;;;;
;
; zernike::Rj(j, rho)
;
; Noll-indexed radial Zernike polynomial
;
function zernike::Rj(j, rho)

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.Rnm(n, m, rho)
end

;;;;;
;
; zernike::Znm(n, m, rho, theta)
;
; Zernike polynomial
;
function zernike::Znm, n, m, rho, theta

COMPILE_OPT IDL2, STATIC

z = sqrt(n+1.) * zernike.Rnm(n, m, rho)

if m ne 0 then $
   z *= sqrt(2.) * ((m gt 0) ? cos(m*theta) : -sin(m*theta))

return, z
end

;;;;;
;
; zernike::Zj(j, rho, theta)
;
; Noll-indexed Zernike polynomials
;
function zernike::Zj, j, rho, theta

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.Znm(n, m, rho, theta)
end

;;;;;
;
; zernike::cZnm(n, m, x, y)
;
; Cartesian Zernike polynomial
; J. Schwiegerling "Description of Zernike Polynomials"
;
function zernike::cZnm, n, m, x, y

COMPILE_OPT IDL2, STATIC

z = 0.
if m ge 0 then begin
   for s = 0, (n - m)/2 do begin
      for j = 0, (n - m)/2 - s do begin
         for k = 0, m/2 do begin
            z += (-1.)^k * zernike.coefficient(n, m, s) * $
                 zernike.binomial((n-m)/2 - s, j) * $
                 zernike.binomial(m, 2.*k) * $
                 x^(n - 2.*(s + j + k)) * y^(2.*(j + k))
         endfor
      endfor
   endfor
   z *= (m eq 0) ? sqrt(n + 1.) : sqrt(2.*n + 2.)
endif else begin
   for s = 0, (n + m)/2 do begin
      for j = 0, (n + m)/2 - s do begin
         for k = 0, -(m + 1)/2 do begin
            z += (-1.)^k * zernike.coefficient(n, m, s) * $
                 zernike.binomial((n+m)/2 - s,  j) * $
                 zernike.binomial(-m, 2.*k + 1) * $
                 x^(n - 2.*(s + j + k) - 1.) * y^(2.*(j + k) + 1.)
         endfor
      endfor
   endfor
   z *= sqrt(2.*n + 2.)
endelse

return, z
end

;;;;;
;
; zernike::cZj(j, x, y)
;
; Noll-indexed Cartesian Zernike polynomial
;
function zernike::cZj, j, x,  y

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.cZnm(n, m, x, y)
end

;;;;;
;
; zernike::dxZnm(n, m, x, y)
;
; Derivative of Cartesian Zernike polynomial along x
; J. Schwiegerling "Description of Zernike Polynomials"
;
function zernike::dxZnm, n, m, x, y

COMPILE_OPT IDL2, STATIC

w = where(abs(x) lt 1e-6, nsmall)
if nsmall gt 0 then x[w] = 1e-6
z = 0.
if m ge 0 then begin
   for s = 0, (n - m)/2 do begin
      for j = 0, (n - m)/2 - s do begin
         for k = 0, m/2 do begin
            z += (n - 2.*(s + j + k)) * $
                 (-1.)^k * zernike.coefficient(n, m, s) * $
                 zernike.binomial((n-m)/2 - s, j) * $
                 zernike.binomial(m, 2.*k) * $
                 x^(n - 2.*(s + j + k) - 1.) * y^(2.*(j + k))
         endfor
      endfor
   endfor
   z *= (m eq 0) ? sqrt(n + 1.) : sqrt(2.*n + 2.)
endif else begin
   for s = 0, (n + m)/2 do begin
      for j = 0, (n + m)/2 - s do begin
         for k = 0, -(m + 1)/2 do begin
            z += (n - 2.*(s + j + k) - 1.) * $
                 (-1.)^k * zernike.coefficient(n, m, s) * $
                 zernike.binomial((n+m)/2 - s,  j) * $
                 zernike.binomial(-m, 2.*k + 1) * $
                 x^(n - 2.*(s + j + k) - 2.) * y^(2.*(j + k) + 1.)
         endfor
      endfor
   endfor
   z *= sqrt(2.*n + 2.)
endelse

return, z
end

;;;;;
;
; zernike::dyZnm(n, m, x, y)
;
; Derivative of Cartesian Zernike polynomial along y
; J. Schwiegerling "Description of Zernike Polynomials"
;
function zernike::dyZnm, n, m, x, y

COMPILE_OPT IDL2, STATIC

w = where(abs(y) lt 1e-6, nsmall)
if nsmall gt 0 then y[w] = 1e-6
z = 0.
if m ge 0 then begin
   for s = 0, (n - m)/2 do begin
      for j = 0, (n - m)/2 - s do begin
         for k = 0, m/2 do begin
            z += 2.*(j + k) * $
                 (-1.)^k * zernike.coefficient(n, m, s) * $
                 zernike.binomial((n-m)/2 - s, j) * $
                 zernike.binomial(m, 2.*k) * $
                 x^(n - 2.*(s + j + k)) * y^(2.*(j + k) - 1.)
         endfor
      endfor
   endfor
   z *= (m eq 0) ? sqrt(n + 1.) : sqrt(2.*n + 2.)
endif else begin
   for s = 0, (n + m)/2 do begin
      for j = 0, (n + m)/2 - s do begin
         for k = 0, -(m + 1)/2 do begin
            z += (2.*(j + k) + 1.) * $
                 (-1.)^k * zernike.coefficient(n, m, s) * $
                 zernike.binomial((n+m)/2 - s,  j) * $
                 zernike.binomial(-m, 2.*k + 1) * $
                 x^(n - 2.*(s + j + k) - 1.) * y^(2.*(j + k))
         endfor
      endfor
   endfor
   z *= sqrt(2.*n + 2.)
endelse

return, z
end

;;;;;
;
; zernike::dxZj(j, x, y)
;
; Noll-indexed derivative of the Cartesian Zernike polynomial along x
;
function zernike::dxZj, j, x, y

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.dxZnm(n, m, x, y)
end

;;;;;
;
; zernike::dyZj(j, x, y)
;
; Noll-indexed derivative of the Cartesian Zernike polynomial along y
;
function zernike::dyZj, j, x, y

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.dyZnm(n, m, x, y)
end

;;;;;
;
; zernike::pZnm(n, m, x, y)
;
; Pseudo-Zernike polynomial with Cartesian arguments
;
function zernike::pZnm, n,  m,  x,  y

COMPILE_OPT IDL2, STATIC

rho = sqrt(x^2 + y^2)
theta = atan(y, x)
return, zernike.Znm(n, m, rho, theta)
end

;;;;; EXPERIMENTAL METHODS

;;;;;
;
; zernike::pZj(j, x, y)
;
; Noll-indexed pseudo-Zernike polynomial with Cartesian arguments
;
function zernike::pZj, j, x, y

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.pZnm(n, m, x, y)
end

;;;;;
;
; zernike::Series(coefficients, rho, theta)
;
; Sum up Noll-indexed Zernike polynomials with
; specified coefficients.
; NOTE: coefficients[0] = 0.
;
function zernike::Series, a, rho, theta

COMPILE_OPT IDL2, STATIC

nterms = n_elements(a)

if nterms le 1 then $
   return, 0.*rho

res = a[1] * zernike.Zj(1, rho, theta)
for n = 2, nterms-1 do $
   res += a[n] * zernike.Zj(n, rho, theta)

return, res
end

;;;;;
;
; zernike::pSeries(coefficients, x, y)
;
; Sum up Noll-indexed pseudo-Zernike polynomials with
; specified coefficients.
;
function zernike::pSeries, a, x, y

COMPILE_OPT IDL2, STATIC

nterms = n_elements(a)

if nterms le 1 then $
   return, 0.*x

res = a[1] * zernike.pZj(1, x, y)
for n = 2, nterms-1 do $
   res += a[n] * zernike.pZj(n, x, y)

return, res
end

;;;;;
;
; zernike::Grnm(n, m, rho, theta)
;
; Functions orthogonal to the radial derivative of the
; Zernike polynomials.
;
; Reference:
; A. Gavrielides, Opt. Lett. 7, 526-528 (1982).
;
function zernike::Grnm, n, m,  rho, theta

COMPILE_OPT IDL2, STATIC

Gxnm = 0.
if m eq 0 then begin
   for k = 0, n/2 do $
      Grnm -= zernike.coefficient(n, m, k)/(n - 2.*k + 2.) * rho^(n - 2.*k + 1.) 
   return, Grnm
endif

for k = 0, (n - m)/2 do $
   Grnm += zernike.coefficient(n, m, k) * $
           (n - 2.*k + 2.)/(n + m - 2.*k + 2.)/(n - m - 2.*k + 2.) * $
           (rho^(m-1.) - rho^(n - 2.*k + 1.))

Grnm *= sqrt(n + 1.) * ((m gt 0) ? cos(m*theta) : -sin(m*theta))
return, Grnm

end

;;;;;
;
; zernike::Gqnm(n, m, rho, theta)
;
; Functions orthogonal to the polar derivative of the
; Zernike polynomials.
;
; Reference:
; A. Gavrielides, Opt. Lett. 7, 526-528 (1982).
;
function zernike::Gqnm, n, m,  rho, theta

COMPILE_OPT IDL2, STATIC

if m eq 0 then $
   return, 0.*rho

Gqnm = 0.
for k = 0, (n - m)/2 do $
   Gqnm += zernike.coefficient(n, m, k) * $
           (n - 2.*k + 2.)/(n + m - 2.*k + 2.)/(n - m - 2.*k + 2.) * $
           ((n - 2.*k + 2)/m * rho^(m-1.) - rho^(n - 2.*k + 1.))

Gqnm *= sqrt(n + 1.) * ((m gt 0) ? cos(m*theta) : -sin(m*theta))
return, Gqnm

end

;;;;;
;
; zernike::Grj(n, m, rho, theta)
;
; Noll-indexed radial Gavrielides polynomial.
;
function zernike::Grj, j, rho, theta

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.Grnm(n, m, rho, theta)
end

;;;;;
;
; zernike::Gqj(n, m, rho, theta)
;
; Noll-indexed polar Gavrielides polynomial.
;
function zernike::Gqj, j, rho, theta

COMPILE_OPT IDL2, STATIC

zernike.noll, j, n, m
return, zernike.Gqnm(n, m, rho, theta)
end

;;;;;
;
; zernike::pGrj(j, x, y)
;
; Noll-indexed radial pseudo-Gavrielides polynomial
;
function zernike::pGrj, j, x, y

COMPILE_OPT IDL2, STATIC

rho = sqrt(x^2 + y^2)
theta = atan(y, x)
return, zernike.Grj(j, rho, theta)
end

;;;;;
;
; zernike::pGqj(j, x, y)
;
; Noll-indexed polar pseudo-Gavrielides polynomial
;
function zernike::pGqj, j, x, y

COMPILE_OPT IDL2, STATIC

rho = sqrt(x^2 + y^2)
theta = atan(y, x)
return, zernike.Gqj(j, rho, theta)
end

;;;;;
;
; zernike::cGj(j, x, y)
;
; Noll-indexed Cartesian Gavrielides polynomial
;
function zernike::cGj, j, x, y

COMPILE_OPT IDL2, STATIC

rho = sqrt(x^2 + y^2)
theta = atan(y, x)
Grj = zernike.Grj(j, rho, theta)
Gqj = zernike.Gqj(j, rho, theta)

return, list(Grj, Gqj)
end

;;;;;
;
; zernike::pSxj(j, x, y)
;
; Cartesian Zhao-Burge polynomials
;
function zernike::pSxj, j, x, y

COMPILE_OPT IDL2, STATIC

case j of
   2: return, zernike.pZj(1, x, y)
   3: return, 0.*x
   4: return, zernike.pZj(2, x, y)/sqrt(2.)
   5: return, zernike.pZj(3, x, y)/sqrt(2.)
   6: return, zernike.pZj(2, x, y)/sqrt(2.)
   7: return, zernike.pZj(5, x, y)/2.
   8: return, zernike.pZj(4, x, y)/sqrt(2.) + zernike.pZj(6, x, y)/2.
   9: return, zernike.pZj(5, x, y)/sqrt(2.)
   10: return, zernike.pZj(6, x, y)/sqrt(2.)
   11: return, zernike.pZj(8, x, y)/sqrt(2.)
   12: return, (zernike.pZj(8, x, y) + zernike.pZj(10, x, y))/2.
   13: return, (zernike.pZj(7, x, y) + zernike.pZj(9, x, y))/2.
   14: return, zernike.pZj(10, x, y)/sqrt(2.)
   15: return, zernike.pZj(9, x, y)/sqrt(2.)
   16: return, zernike.pZj(11, x, y)/sqrt(2.) + zernike.pZj(12, x, y)/2.
   17: return, zernike.pZj(3, x, y)/2.
   18: return, (zernike.pZj(12, x, y) + zernike.pZj(14, x, y))/2.
   19: return, (zernike.pZj(13, x, y) + zernike.pZj(15, x, y))/2.
   else: return, 0.*x
endcase

end

;;;;;
;
; zernike::pSyj(j, x, y)
;
; Cartesian Zhao-Burge polynomials
;
function zernike::pSyj, j, x, y

COMPILE_OPT IDL2, STATIC

case j of
   2: return, 0.*x
   3: return, zernike.pZj(1, x, y)
   4: return, zernike.pZj(3, x, y)/sqrt(2.)
   5: return, zernike.pZj(2, x, y)/sqrt(2.)
   6: return, -zernike.pZj(3, x, y)/sqrt(2.)
   7: return, zernike.pZj(4, x, y)/sqrt(2.) - zernike.pZj(6, x, y)/2.
   8: return, zernike.pZj(5, x, y)/2.
   9: return, zernike.pZj(6, x, y)/sqrt(2.)
   10: return, -zernike.pZj(5, x, y)/sqrt(2.)
   11: return, zernike.pZj(7, x, y)/sqrt(2.)
   12: return, (-zernike.pZj(7, x, y) + zernike.pZj(9, x, y))/2.
   13: return, (zernike.pZj(8, x, y) - zernike.pZj(10, x, y))/2.
   14: return, -zernike.pZj(9, x, y)/sqrt(2.)
   15: return, zernike.pZj(10, x, y)/sqrt(2.)
   16: return, zernike.pZj(3, x, y)/2.
   17: return, -zernike.pZj(11, x, y)/sqrt(2.) + zernike.pZj(12, x, y)/2.
   18: return, (-zernike.pZj(13, x, y) + zernike.pZj(15, x, y))/2.
   19: return, (zernike.pZj(12, x, y) - zernike.pZj(14, x, y))/2.
   else: return, 0.*x
endcase

end

;;;;;
;
; zernike::pZSproject(maxj, dx, dy, x, y)
;
; Project gradient field onto Zhao-Burge polynomials
;
function zernike::pZBproject, maxj, dx, dy, x, y

COMPILE_OPT IDL2, STATIC

alpha = fltarr(maxj+1)
for j = 2, maxj do $
   alpha[j] = total(dx * zernike.pSxj(j, x, y) + dy * zernike.pSyj(j, x, y))

zernike.noll, [0:maxj], n, m
return, alpha/sqrt(2.*n*(n+1))
end

;;;;; CLASS DEFINITION

;;;;;
;
; zernike__define
;
pro zernike__define

COMPILE_OPT IDL2, HIDDEN

void = {zernike, $
        inherits IDL_Object}

end
