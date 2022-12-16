Fs = 100;
f=1;
t = [1:1000]/Fs;
sig = sin(2*pi*f*t)+t*0.1;

[pks,locs] = findpeaks(sig, 'MinPeakHeight',1.6);

figure()
hold on
plot(t,sig)
scatter(t(locs),pks)