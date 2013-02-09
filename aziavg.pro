;+
; NAME:
;    aziavg
;
; PURPOSE:
;    Average a two dimensional data set over angles about its center.
;
; CATEGORY:
;    Image Processing
;
; CALLING SEQUENCE:
;    result = aziavg( data )
;
; INPUTS:
;    data: two dimensional array of any type except string or complex
;
; KEYWORD PARAMETERS:
;    center: coordinates of center: [xc,yc].  Default is to use data's
;        geometric center.
;
;    rad: maximum radius of average [pixels]
;        Default: half the minimum dimension of the image.
;
;    weight: relative weighting of each pixel in data.
;        Default: uniform weighting.
;
;    deinterlace: If set to an even number, average only over even 
;        numbered lines.  Similarly if set to an odd number.
;        This is useful for analyzing interlaced video images.
;
; OUTPUTS:
;    result: data averaged over angles as a function of radius from
;        the center point, measured in pixels.  Result is single precision.
;
; PROCEDURE:
;    data[x,y] sits at radius r = sqrt((x-xc)^2 + (y-yc)^2) 
;        from the center, (xc,yc).  Let R be the integer part
;        of r, and dR the fractional part.  Then this point is
;        averaged into result(R) with a weight 1-dR and into
;        result(R+1) with a weight dR.
;
; RESTRICTIONS:
;    data must be two-dimensional and must not be string type
;
; MODIFICATION HISTORY:
; Written by:	David G. Grier, The University of Chicago, 07/30/1992
; 06/1994 DGG Handles complex type correctly
; 08/1999 DGG Added keyword CENTER and modernized array notation
; 01/27/2009 DGG Added DEINTERLACE keyword.  Cleaned up formatting.
;     Documentation fixes.
; 02/03/2009 DGG Return -1 if no pixels are within range of
;     the center.  Added RAD keyword.
; 06/12/2012 DGG Added WEIGHT keyword.  Added COMPILE_OPT.  Modernized
;     mathematical operations.  Improved error checking on DATA.
;     Correct complex average for double-precision.
;     Documentation upgrades.
; 07/17/2012 DGG Use rebin to calculate distance array.
; 01/24/2013 DGG do not deinterlace if deinterlace=0
;
; Copyright (c) 1992-2013 David G. Grier
;-
function aziavg, _data, $
                 center = center, $
                 rad = rad, $
                 weight = weight, $
                 deinterlace = deinterlace

COMPILE_OPT IDL2

on_error, 2			; return to calling routine on error

if ~isa(_data, /number, /array) then begin
   message, 'USAGE: result = aziavg(data)', /inf
   message, 'data: must be an array of numbers', /inf
   return, -1
endif

sz = size(_data)
if sz[0] ne 2 then begin
   message, 'Requires 2-dimensional data set', /inf
   return, -1
endif
nx = sz[1]			; width
ny = sz[2]			; height 

if n_elements(center) eq 2 then begin
   xc = double(center[0])
   yc = double(center[1])
endif else begin
   xc = 0.5D * double(nx - 1)   ; indices start at 0
   yc = 0.5D * double(ny - 1)   ; ... in y also
endelse

if isa(rad, /number, /scalar) then $ ; maximum radius
   rmax = round(rad) $
else $
   rmax = nx/2 < ny/2

if (sz[3] eq 6) or (sz[3] eq 9) then begin ; complex data
   a = dcomplex(_data)
   sum = dcomplexarr(rmax + 1)
endif else begin                ; accumlate other types into double
   a = double(_data)
   sum = dblarr(rmax + 1)
endelse

n = dblarr(rmax + 1)

; distance from center to each pixel
r = rebin((dindgen(nx) - xc)^2, nx, ny) + $
    rebin((dindgen(1, ny) - yc)^2, nx, ny)

if keyword_set(deinterlace) then begin
   n0 = deinterlace mod 2
   a = a[*,n0:*:2]
   r = r[*,n0:*:2]
endif

; limit average to maximum range
w = where(r lt rmax^2, ngood)
if ngood lt 0 then $
   return, -1                   ; no pixels in average

r = sqrt(r[w])                  ; only consider data in range
dl = a[w]                       ; data in range

; apportion data proportionally into radius bins above and below 
; target radius
ri = long(r)                    ; integer index for r
fh = r - ri                     ; fraction in higher bin
fl = 1.D - fh                   ; fraction in lower bin
if arg_present(weight) then begin
   if ~isa(weight, /number) then begin
      message, 'weight: must be a numerical data type', /inf
      return, -1
   endif
   if n_elements(weight) ne nx*ny then begin
      message, 'weight: must have same number of elements as data', /inf
      return, -1
   endif
   fh *= weight[w]
   fl *= weight[w]
endif
dh = dl * fh                    ; apportion fractions ...
dl = dl * fl

for i = 0L, ngood-1 do begin	; loop through data points in range
   ndx = ri[i]                  ; lower bin
   sum[ndx]   += dl[i] 
   n[ndx]     += fl[i]
   sum[ndx+1] += dh[i]
   n[ndx+1]   += fh[i]
endfor

return, sum/n			; normalize by number in each bin
end
