sound(tan(2.^[1:0.0001:25]*3))

sound(cos(log([1e5:1e1:1e6]*1e1)*1e4)./cos(([1e5:1e1:1e6])*1e-4))
sound(cos(([1e5:1e1:1e6]*1e1)*1e-2).*cos(([1e5:1e1:1e6])*1e-4))

sound(cos(.3.^cos(3.^[5:-5e-5:1]*50)*50).*log([1:5e-5:5]*30))


sound(cos(.2.^cos(2.^[4:-2e-5:0]*50)*50).*log([1:2e-5:5]*50))

sound(sinc([1:0.0001:25]*1e5)*1e18)

sound(sinc([1:0.001:25].*1e3)*1e16)

plot(sinc([1:0.001:50].*1e3).*sin([1:0.001:50])*1e16)


T = 3; t = linspace(0,T,Fs*T); x = cos(2*pi*875*t)/2.+cos(2*pi*885*t)/2.+cos(2*pi*880*t);plot(x);sound(x);
x = cos(2*pi*880*t).*exp(-500*t/(T));plot(x);sound(x);