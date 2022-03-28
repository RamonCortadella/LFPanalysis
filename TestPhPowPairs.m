clearvars
close('all')
FsDC = 500;
Fs = 500;

tend= 100;
sleep_dc = sin(linspace(0,tend,tend*FsDC)*2*pi*0.1+pi/2);

sleep_ac = sin(linspace(0,tend,tend*Fs)*2*pi*0.1+pi/2);
indexes =find(sleep_ac<=0);
sleep_ac(indexes) = 0;
myrLfp = sleep_ac.*sin(linspace(0,tend,tend*Fs)*2*pi*10);
% 
% fPow = bsxfun(@plus,[10],[-1 1]')'/50; 
% fPh = bsxfun(@plus,[0.1], [-0.02 0.02]')'/50;
FNi = Fs/2;
fPow = bsxfun(@plus,[1:0.5:35],[-0.5 0.5]')'/FNi;
fPh = bsxfun(@plus,[0.02:0.005:0.4], [-0.01 0.01]')'/FNi;
% fPh = [0.07 0.13]./FNi;
% fPow = [9 11]./FNi;
display(fPh)
[out, PhMat, PowMat] = PowerPhasePairsR(sleep_dc',  fPh, myrLfp' , fPow, 1,'but',@PowerModulation);

figure(1)

plot(sleep_dc)
hold on
plot(myrLfp)
hold off

figure(2)
pcolor(out.fPh*FNi, out.fPow*FNi, out.Ramp(:,:)'); shading flat; 
colorbar

figure(3)
polarplot(out.phbins(:,1,1),sq(out.pow_dens(:,10,19))'); axis tight; box off
%%
% ph = angle(hilbert(ButFilter(xPh, 2, FrBinsPh(kPh,:)*ResampleNum, 'bandpass')));
% filtISA = ButFilter(sleep_dc, 2, fPh, 'bandpass');
% figure()
% plot(filtISA)

% figure;
% dcph = angle(hilbert(sleep_dc));
% for k=1:8
%     subplot(3,3,k);
%     hist2(([[dcph(:,k) sppow(:,2)];[dcph(:,k)+2*pi sppow(:,2)]]),20,20,'xprob','mud'); caxis([0 0.1]);
% end