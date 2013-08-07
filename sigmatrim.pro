;+
; NAME:
;    sigmatrim
;
; PURPOSE:
;    Remove outlying values from a data set by requiring that
;    all values fall within SIGMA standard deviations of the mean.
;
; CATEGORY:
;    Math, statistics
;
; CALLING SEQUENCE:
;    sigmatrim, data, mean, dev, signal
;
; INPUTS:
;    data: The data to be sigmatrimmed
;
; KEYWORD PARAMETERS:
;    top: set if larger values are to be retained
;    bottom: set if smaller values are to be retained
;        Default is to trim both top and bottom.
;
;    sigma: Number of standard deviations to use for trimming.
;        Default is fixed between 1 and 3 based on a quick
;        look at the spread of the initial data.
;
; OUTPUTS:
;    mean: mean of the trimmed data
;    dev: standard deviation of the trimmed data
;    signal: indices of the retained data
;
; MODIFICATION HISTORY:
; Written by David G. Grier, AT&T Bell Laboratories, 1/92
; 03/17/2013 DGG New York University Formatting, code clean-up and
;   usage message.
; 08/07/2013 DGG test for non-convergence.
;
; Copyright (c) 1992-2013 David G. Grier
;-
pro sigmatrim, data, mean, dev, signal, $
               top = top, bottom = bottom, sigma = sigma

COMPILE_OPT IDL2

umsg = "USAGE: sigmatrim, data, mean, [dev, [signal]]"

if n_params() lt 2 then begin
   message, umsg, /inf
   return
endif

dev = stdev(data, mean)
if not keyword_set(sigma) then $
   sigma = (float(fix(10. * (max(data) - mean) / dev)) * 0.1) < 3. > 1.
nsignal = n_elements(data)
repeat begin
   osignal = nsignal
   lim = sigma * dev
   if keyword_set(top) then $
      signal = where(data gt mean - lim, nsignal) $
   else if keyword_set(bottom) then $
      signal = where(data lt mean + lim, nsignal) $
   else $
      signal = where((data gt mean - lim) and (data lt mean + lim), nsignal)
   if nsignal lt 1 then begin
      message, 'did not converge', /inf
      break
   endif
   dev = stdev(data[signal], mean)
end until(osignal eq nsignal)
end	
