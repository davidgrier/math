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
;    result = aziavg(data)
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
;        the center point, measured in pixels.  Result is single
;        precision.
;
; KEYWORD OUTPUTS:
;    rho: the radial position of each pixel in DATA relative to the
;        center at (xc,yc).
;
;    values: Azimuthal average at each pixel.
;
;    deviates: difference between DATA and azimuthal average at each
;        pixel.
;
; PROCEDURE:
;    data[x,y] sits at radius rrho = sqrt((x-xc)^2 + (y-yc)^2) 
;        from the center, (xc,yc).  Let R be the integer part
;        of rho, and dR the fractional part.  Then this point is
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
; 03/04/2013 DGG Fix corner case when bins have no counts.
; 03/17/2013 DGG updated usage message.  Small speed-up in
;     accumulation loop.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 03/24/2013 DGG small efficiency improvements.
; 05/05/2013 DGG Use HISTOGRAM for computations.  Major speed-up.
; 05/19/2013 DGG # is faster than rebin(/sample)
; 06/02/2013 DGG Added RHO keyword.  Added DEVIATES keyword.  Return
;   deviates for all points, even when called with deinterlace
; 08/05/2013 DGG fix average for small r when deinterlacing
; 08/18/2013 DGG added VALUES keyword
;
; Copyright (c) 1992-2013 David G. Grier
;-
function aziavg, _data, $
                 center = center, $
                 rad = rad, $
                 weight = weight, $
                 deinterlace = deinterlace, $
                 rho = rho, $
                 deviates = deviates, $
                 values = values

COMPILE_OPT IDL2

umsg = 'USAGE: result = aziavg(data)'

if ~isa(_data, /number, /array) then begin
   message, umsg, /inf
   message, 'DATA must be a numeric array', /inf
   return, -1
endif

sz = size(_data)
if sz[0] ne 2 then begin
   message, umsg, /inf
   message, 'DATA must be a two-dimensional array', /inf
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

rmax = isa(rad, /number, /scalar) ? round(rad) : (nx/2 < ny/2)

if (sz[3] eq 6) or (sz[3] eq 9) then begin ; complex data
   a = dcomplex(_data)
   sum = dcomplexarr(rmax+1)
endif else begin                ; accumlate other types into double
   a = double(_data)
   sum = dblarr(rmax+1)
endelse
count = dblarr(rmax+1)

if isa(weight, /number, /array) then $
   a *= weight

; distance from center to each pixel
r = (dindgen(nx) - xc)^2 # replicate(1., ny) + $
      replicate(1., nx) # (dindgen(ny) - yc)^2
r = sqrt(temporary(r))

if keyword_set(deinterlace) then begin
   n0 = deinterlace mod 2
   a = a[*, n0:*:2]
   r = r[*, n0:*:2]
endif

rn = floor(r)
fh = r - rn
fl = 1.d - fh
ah = a * fh
al = a * fl
h = histogram(r, min = 0, max = rmax+1, reverse_indices = n)
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

avg = sum/(count > 1e-3)

if arg_present(values) or arg_present(deviates) then begin
   values = fl*avg[rn] + fh*avg[rn+1]
   if arg_present(deviates) then $
      deviates = a - values
endif

return, avg
end
