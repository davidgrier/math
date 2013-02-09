;+
; NAME:
;    IQR
;
; PURPOSE:
;    Compute the inter-quartile range of a one-dimensional data set
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    r = iqr(x)
;
; INPUTS:
;    x: one-dimensional numerical data set
;
; OUTPUTS:
;    r: inter-quartile range of x
;
; RESTRICTIONS:
;    Considers only finite elements of x.
;    Considers only one-dimensional data sets.
;
; PROCEDURE:
;    The data are sorted into ascending order.  The lower quartile
;    is the median of the first half of the ordered data and the
;    upper quartile is the median of the second half.  The
;    interquartile range is the difference between the upper and
;    lower quartile values.
;
; MODIFICATION HISTORY:
; 09/18/2010 Written by David G. Grier, New York University
;
; Copyright (c) 2010 David G. Grier
;
;-

function iqr, data

; only consider valid numerical data
w = where(finite(data) eq 1, ndata)
if ndata lt 1 then begin
   message, "no finite data points", /inf
   return, 0.
endif
x = data[w]

min = min(data, max=max)

if min eq max then $
   return, 0.

x = x[sort(x)]
nhi = ndata/2
; if ndata is odd, the middle data point belongs to
; both the upper and lower quartiles ...
nlo = (ndata mod 2) eq 0 ? nhi - 1 : nhi

return, median(x[nhi:*], /even) - median(x[0:nlo], /even)

end
