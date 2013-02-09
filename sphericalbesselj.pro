function sphericalbesselj, ndx, xin

x = double(xin)

jnm1 = sin(x)/x ; j_0(x)

if ndx eq 0 then return, jnm1

jn = (jnm1 - cos(x))/x ; j_1(x)

if ndx eq 1 then return, jn

for n = 2.d, ndx do begin
    jnm2 = jnm1
    jnm1 = jn
    jn = ((2.d*n - 1.d) * jnm1 - x * jnm2)/x
endfor
return, jn
end
