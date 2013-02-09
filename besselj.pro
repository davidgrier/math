function besselj0, x

ax = abs(x)
wl = where(ax gt 8, nlarge, complement=ws, ncomplement=nsmall)

ans = x * 0.d

if (nsmall gt 0) then begin
   y = x[ws]^2
   ans1 = 57568490574.d + y * (-13362590354.d + $
                          y * (651619640.7d + $
                          y * (-11214424.18d + $
                          y * (77392.33017d + $
                          y * (-184.9052456)))))
   ans2 = 57568490411.d + y * (1029532985.d + $
                          y * (9494680.718d + $
                          y * (59272.64853 + $
                          y * (267.8532712d + $
                          y))))
   ans[ws] = ans1/ans2
endif
if (nlarge gt 0) then begin
   z = 8.0d / ax[wl]
   y = z^2
   xx = ax[wl] - 0.785398164d
   ans1 = 1.d + y * (-0.1098628627d-2 + $
                y * (0.2734510407d-4 + $
                y * (-0.2073370639d-5 + $
                y * (0.2093887211d-6))))
   ans2 = -0.1562499995d-1 + y * (0.1430488765d-3 + $
                             y * (-0.6911147651d-5 + $
                             y * (0.7621095161d-6 - $
                             y * 0.934935152d-7)))
   ans[wl] = sqrt(0.636619772/ax[wl]) * (cos(xx)*ans1 - z * sin(xx)*ans2)
endif

return, ans
end
