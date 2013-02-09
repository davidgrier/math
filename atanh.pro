function atanh, x

if max(abs(x)) ge 1 then message, "argument out of range"
return, 0.5 * alog((1. + x)/(1. - x))
end
