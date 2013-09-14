;+
; NAME:
;    azimedian
;
; PURPOSE:
;    Compute the azimuthal median of a two-dimensional data set
;    over angles about its center.
;
; CATEGORY:
;    Image Processing
;
; CALLING SEQUENCE:
;    result = azimedian(data)
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
;    values: Azimuthal median at each pixel.
;
;    deviates: difference between DATA and azimuthal median at each
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
; 08/18/13 Written by David G. Grier, New York University
; 09/14/13 DGG Compute value at largest radius.  Apportion
;    values and deviates proportionately.
;
; Copyright (c) 2013 David G. Grier
;-
function azimedian, a, $
                 center = center, $
                 rad = rad, $
                 deinterlace = deinterlace, $
                 rho = rho, $
                 deviates = deviates, $
                 values = values

COMPILE_OPT IDL2

umsg = 'USAGE: result = azimedian(data)'

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   message, 'DATA must be a numeric array', /inf
   return, -1
endif

sz = size(a)
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
   med = dcomplexarr(rmax+1)
endif else begin                ; accumulate other types into double
   med = dblarr(rmax+1)
endelse

if isa(weight, /number, /array) then $
   a *= weight

; distance from center to each pixel
r = (dindgen(nx) - xc)^2 # replicate(1., ny) + $
      replicate(1., nx) # (dindgen(ny) - yc)^2
r = sqrt(temporary(r))

if keyword_set(deinterlace) then begin
   n0 = deinterlace mod 2
   b = a[*, n0:*:2]
   r = r[*, n0:*:2]
endif else $
   b = a

h = histogram(r, min = 0, max = rmax+1, reverse_indices = n)
for i = 0L, rmax do begin
   n0 = n[i]
   n1 = n[i+1]-1
   if (n1 ge n0) then $
      med[i] = median(b[n[n0:n1]])
endfor

if arg_present(values) or arg_present(deviates) then begin
   rn = floor(r)
   fh = r - rn
   fl = 1. - fh
   values = fl*med[rn] + fh*med[rn+1]
   if arg_present(deviates) then $
      deviates = b - values
endif

return, med
end
