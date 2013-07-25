;+
; NAME:
;    azistd
;
; PURPOSE:
;    Calcuate the standard deviation from the mean of a two
;    dimensional data set, averaged over angles, as a function of
;    radius from the center.
;
; CATEGORY:
;    Image Processing
;
; CALLING SEQUENCE:
;    result = azistd(data, [avg], rad=rad, center=center)
;
; INPUTS:
;    data: two dimensional array of any type except string or complex
;
; KEYWORD PARAMETERS:
;    center: coordinates of center: [xc,yc].
;        Default: geometric center.
;
;    rad: maximum radius over which to compute [pixel]
;        Default: half of the minimum dimension
;
;    deinterlace: If set to an even number, compute azimuthal standard
;        deviation only over even number lines.  Similarly for odd
;        values.  Useful for analyzing interlaced video images.
;
; OUTPUTS:
;    result: standard deviation of data averaged over angles as a 
;        function of radius from the center point, measured in pixels.  
;        Result is single precision.
;
; OPTIONAL OUTPUTS:
;    avg : azimuthal average of the data
;
;    rho: distance of each pixel in DATA from center at (xc,yc).
;
; PROCEDURE:
;    data(x,y) sits at radius rho = sqrt( (x-xc)^2 + (y-yc)^2 ) 
;    from the center, (xc,yc).  Let R be the integer part
;    of rho, and dR the fractional part.  Then this point is
;    averaged into result(R) with a weight 1-dR and into
;    result(R+1) with a weight dR.
;
; RESTRICTIONS:
;	data must be two-dimensional and must not be string type
;
; MODIFICATION HISTORY:
; Written by David G. Grier, The University of Chicago, 07/30/1992
; 06/1994 DGG Handles complex type correctly
; 08/1999 DGG Added keyword CENTER and modernized array notation
; 04/2002 DGG Converted from aziavg
; 01/27/2009 DGG Added deinterlace keyword.  Documentation clean up.
; 01/27/2013 DGG Use rebin to calculate distance arrays.
;   Correctly handle deinterlace=0.  Added RAD keyword.
; 03/04/2013 DGG Added COMPILE_OPT.  Fix corner cases when bins have
;   no counts.
; 03/17/2013 DGG Small speedup in accumulation loops.  Updated usage
;   message and error testing.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 03/24/2013 DGG small efficiency improvements.
; 05/05/2013 DGG Use HISTOGRAM for calculation.  Major speed-up.
; 05/19/2013 DGG # is faster than rebin(/sample).
; 06/02/2013 DGG Added RHO keyword.
; 07/24/2013 DGG Fix critical typo
;
; Copyright (c) 1992-2013 David G. Grier
;-
function azistd, _data, avg, $
                 center = center, $
                 rad = rad, $
                 rho = rho, $
                 deinterlace = deinterlace

COMPILE_OPT IDL2

umsg = 'USAGE: result = azistd(data, [avg])'

if ~isa(_data, /number, /array) then begin
   message, umsg, /inf
   message, 'DATA must be a numeric array', /inf
   return, -1
endif
sz = size(_data)
if sz[0] ne 2 then begin
   message, umsg, /inf
   message, 'DATA must be two-dimensional', /inf
   return, -1
endif
nx = sz[1]			; width
ny = sz[2]			; height 

; center point
if n_elements(center) eq 2 then begin
   xc = double(center[0])
   yc = double(center[1])
endif else begin
   xc = 0.5D * double(nx - 1)	; indices start at 0
   yc = 0.5D * double(ny - 1)	; ... in y also
endelse

if isa(rad, /number, /scalar) then $ ; maximum radius
   rmax = round(rad) $
else $
   rmax = nx/2 < ny/2

if (sz[3] eq 6) or (sz[3] eq 9) then begin ; complex data
   a = dcomplex(_data)
   sum = dcomplexarr(rmax+1)
endif else begin                ; accumulate other types into double
   a = double(_data)
   sum = dblarr(rmax+1)
endelse
count = dblarr(rmax + 1)

rho = (dindgen(nx) - xc)^2 # replicate(1., ny) + $
      replicate(1., nx) # (dindgen(1, ny) - yc)^2

if keyword_set(deinterlace) then begin
   n0 = deinterlace mod 2
   a = a[*, n0:*:2]
   rho = rho[*, n0:*:2]
endif

rho = sqrt(rho)
fh = rho - floor(rho)
fl = 1.d - fh
ah = a * fh
al = a * fl
h = histogram(rho, min = 0, max = rmax+1, reverse_indices = n)
for i = 0L, rmax-1 do begin
   n0 = n[i]
   n1 = n[i+1]-1
   if (n1 ge n0) then begin
      ndx = n[n0:n1]
      sum[i] += total(al[ndx])
      count[i] += total(fl[ndx])
      sum[i+1] = total(ah[ndx])
      count[i+1] = total(fh[ndx])
   endif
endfor
count >= 1.d
avg = sum/count                 ; normalize by number in each bin

sum *= 0                        ; reset sum for standard deviation
for i = 0L, rmax-1 do begin     ; loop through data points in range
   n0 = n[i]
   n1 = n[i+1]-1
   if (n1 ge n0) then begin
      ndx = n[n0:n1]
      dsq = (a[ndx]-avg[i])^2
      sum[i] += total(fl[ndx]*dsq)
      sum[i+1] = total(fh[ndx]*dsq)
   endif
endfor
std = sqrt(sum/count)               ; standard deviation

return, std
end
