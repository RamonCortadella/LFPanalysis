Par=  LoadXml('../2019-07-31T15-51-13B12784O18-T2-LongTerm-Rec1/2019-07-31T15-51-13B12784O18-T2-LongTerm-Rec1.xml');

d = dir('../*/*.dat');

%FileBase = [ComName num2str(fn)];



%% get post.chan
HVS_Chan = 41;
The_Chan = 58; rLfp = [];
for fn = 1:5
    FileName = d(fn).name;
    Lfp = LoadBinary(FileName, [The_Chan HVS_Chan], Par.nChannels,2)';
    myrLfp = resample(Lfp,1,10); % RESAMPLE TO 100 HZ !!!
    rLfp = [rLfp; myrLfp];
    dur(fn) = length(myrLfp);
end

%% load dc
%dc_Chan = [65 63 6 40 32];
dc_Chan = [65 63 6 8 30 40 32 38];
rdc = [];
for fn = 1:5
    FileName = d(fn).name;
    dc = LoadBinary(FileName, dc_Chan, Par.nChannels,2)';
    myrdc= resample(dc,1,100); % resample to 10 hz
    rdc = [rdc; myrdc];
end
rdcf = ButFilter(rdc,4,0.1/5,'low');

%% load movement
for fn = 1:5
    FileName = d(fn).name;
    mymv = load([FileName(1:end-3) 'MovVar.mat']);
    mv(fn)=mymv.MovVar;
    durmv(fn) = size(mv(fn).HeadSpeed,1);
end

% now concattenate with taking into account that movement file is shorter
% than lfp
cmv = CatStruct(mv);

%% load states
% sleep

cdur = cumsum(dur); cdur = [0; cdur(1:end-1)']; % this is in 100hz samples
clear States
States.Labels = {'SWS','REM','HVS','AWKTH','AWKNOTH','REAR'};
%lls = dir('*1.sts.*'); lls = lls(2:end);
% cnt=1;
% for k=1:length(lls) 
%     lab = (lls(k).name(max(findstr(lls(k).name,'.'))+1:end));
%     if ~( strcmp(lab,'SLEEP')  |   strcmp(lab,'RUN') |  strcmp(lab,'QUIET'))
%         States.Labels{cnt} = lab;
%         cnt=cnt+1;
%     end
% end
    
for fn = 1:5
    FileName = d(fn).name;
    for s=1:length(States.Labels)
        if fn==1 States.(States.Labels{s}) = []; end
        StateFileName = [FileName(1:end-3) 'sts.' States.Labels{s}];
        if FileExists(StateFileName)
            myPer = load(StateFileName); % is in seconds ..  
            States.(States.Labels{s}) = [States.(States.Labels{s}); cdur(fn) + round(myPer/1000*100) ]; % 100hz s.rate
        end
    end
end
States.nLabels = length(States.Labels);

%%

myrLfp = SelectPeriods(rLfp, States.SWS, 'c',1);
fLfp = ButFilter(myrLfp,4,[8 16]/50,'bandpass');
sppow = abs(hilbert(fLfp));
% sppow = resample(sppow,1,10);
%% load dc
%dc_Chan = [65 63 6 40 32];
dc_Chan = [65 63 6 8 30 40 32 38];
rdc100 = [];
for fn = 1:5
    FileName = d(fn).name;
    dc = LoadBinary(FileName, dc_Chan, Par.nChannels,2)';
    myrdc= resample(dc,1,10); % resample to 100 hz
    rdc100 = [rdc100; myrdc];
end
rdc100f = ButFilter(rdc100,4,0.1/50,'low');

sleep_dc = SelectPeriods(rdc100f,round(States.SWS),'c',1);
%sleep_dc = resample(sleep_dc,1,10);


%% 
fPow = [bsxfun(@plus,[1:1:15],[-0.5 0.5]')'; bsxfun(@plus, [15:5:40],[-5 5]')']/50;
fPh = bsxfun(@plus,[0.03:0.02:0.5], [-0.02 0.02]')'/50;
[out, PhMat, PowMat] = PowerPhasePairs(sleep_dc,  fPh, myrLfp , fPow,1,'but',@PowerModulation);
%%
figure
for k=1:2
    for l=1:8
        subplot2(2,8, k,l);
        pcolor(out.fPh*50, out.fPow*50, out.Ramp(:,:,k,l)'); shading flat; 
    end;
end


%%
figure;
dcph = angle(hilbert(sleep_dc));
for k=1:8
    subplot(3,3,k);
    hist2(([[dcph(:,k) sppow(:,2)];[dcph(:,k)+2*pi sppow(:,2)]]),20,20,'xprob','mud'); caxis([0 0.1]);
end
%% same for theta and REM
myrLfp1 = SelectPeriods(rLfp, States.REM, 'c',1);
sleep_dc1 = SelectPeriods(rdc100f,round(States.REM),'c',1);
fPow1 = [bsxfun(@plus,[1:0.5:15],[-0.3 0.3]')']/50;
[out1, PhMat1, PowMat1] = PowerPhasePairs(sleep_dc1,  fPh, myrLfp1 , fPow1,1,'but',@PowerModulation);

%%
figure
for k=1:2
    for l=1:8
        subplot2(2,8, k,l);
        pcolor(out1.fPh*50, out1.fPow*50, out1.Ramp(:,:,k,l)'); shading flat; 
    end;
end

%%
fLfp = ButFilter(myrLfp,4,[5 12]/50,'bandpass');
thpow = abs(hilbert(fLfp));
thpow = resample(thpow,1,10);

sleep_dc = resample(sleep_dc,1,10);
dcph = angle(hilbert(sleep_dc));
figure;
%dcph = angle(hilbert(sleep_dc));
for k=1:8
    subplot(3,3,k);
    hist2(([[dcph(:,k) thpow(:,2)];[dcph(:,k)+2*pi thpow(:,2)]]),20,20,'xprob','mud'); caxis([0 0.1]);
end

%% now load all channels
rLfpAll = [];
for fn = 1:5
    FileName = d(fn).name;
    Lfp = LoadBinary(FileName, RepCh, Par.nChannels,2)';
    myrLfp = resample(Lfp,1,10); % RESAMPLE TO 100 HZ !!!
    rLfpAll= [rLfpAll; myrLfp];
    
end

%% REM sleep
nT = size(rLfpAll,1);
In = WithinRanges([1:nT]',States.REM) & ~WithinRanges([1:nT]',States.HVS);
myrLfpAll = rLfpAll(In,:);
sleep_dc1 = rdc100f(In,:);


fPow1 = [bsxfun(@plus,[1:0.5:15],[-0.3 0.3]')']/50;
fPh1 = bsxfun(@plus,[0.03:0.02:0.2], [-0.02 0.02]')'/50;
%[out2] = PowerPhasePairs(sleep_dc1,  fPh1, myrLfpAll , fPow1,1,'but',@PowerModulation);

Shuffle.MaxShift= 100*100/10; % after 10 times resampling will have 50sec shift.
Shuffle.Type = 'shift'; Shuffle.nShuffle = 1000;
[out2sh] = PowerPhasePairs(sleep_dc1(:,7),  fPh1, myrLfpAll(:,13) , fPow1,10,@PowerModulation, Shuffle);

fPow = [bsxfun(@plus,[1:1:15],[-0.5 0.5]')'; bsxfun(@plus, [15:5:40],[-5 5]')']/50;
[out2gam] = PowerPhasePairs(sleep_dc1(:,7),  fPh1, myrLfpAll(:,13) , fPow,10,@PowerModulation);

%%compute for gamma 


%% 

%% 
figure
subplot(221);
pcolor(out2.fPow*50, [1:64], sq(out2.Ramp(4,:,:,4))'); shading flat; 

subplot(222);
pcolor(out2.fPow*50, [1:64], sq(out2.Rth(4,:,:,4))'); shading flat; CircColormap

subplot(223);
pcolor(out2.fPow*50, [1:64], sq(out2.Ramp(4,:,:,6))'); shading flat; 

subplot(224);
pcolor(out2.fPow*50, [1:64], sq(out2.Rth(4,:,:,6))'); shading flat; CircColormap


%% loop through freqs and plot over phase channels
for k=1:size(fPow1,1);
    figure(3323);clf
    subplot(211);
    pcolor([1:64], [1:7], sq(out2.Ramp(4,k,:,[1:4 6:8]))'); shading flat; 

    subplot(212);
    pcolor([1:64], [1:7], sq(out2.Rth(4,k,:,[1:4 6:8]))'); shading flat; CircColormap
    title(num2str(out2.fPow(k)*50));
    waitforbuttonpress
end


%%
figure
for k=1:8
   subplot2(2,8,1,k);
    pcolor( out2.fPh*50, out2.fPow*50, out2.Ramp(:,:,11,k)');
    
    subplot2(2,8,2,k);
    pcolor( out2.fPh*50, out2.fPow*50, out2.Rth(:,:,11,k)'); CircColormap; shading flat;
end

%%
nT = size(rLfpAll,1);
In = WithinRanges([1:nT]',States.AWKTH);
myrLfpAll = rLfpAll(In,:);
sleep_dc1 = rdc100f(In,:);


fPow1 = [bsxfun(@plus,[1:0.5:15],[-0.3 0.3]')']/50;
fPh1 = bsxfun(@plus,[0.03:0.02:0.2], [-0.02 0.02]')'/50;
[out3] = PowerPhasePairs(sleep_dc1,  fPh1, myrLfpAll , fPow1,1,'but',@PowerModulation);

figure
clf
for k=1:8
    subplot2(3,8,1,k);
    pcolor( out3.fPh*50, out3.fPow*50, out3.MI(:,:,11,k)');
   
   subplot2(3,8,2,k);
    pcolor( out3.fPh*50, out3.fPow*50, out3.Ramp(:,:,11,k)');
    
    subplot2(3,8,3,k);
    pcolor( out3.fPh*50, out3.fPow*50, out3.Rth(:,:,11,k)'); CircColormap; shading flat;
end

for k=1:size(fPow1,1);
    figure(3323);clf
    subplot(211);
    imagesc([1:64], [1:7], sq(mean(out3.Ramp(2:5,k,:,[1:4 6:8])))'); colorbar

    subplot(212);
    imagesc([1:64], [1:7], sq(circmean(out3.Rth(2:5,k,:,[1:4 6:8])))');  CircColormap; colorbar
    title(num2str(out3.fPow(k)*50));
    waitforbuttonpress
end

%theta_power = abs(hilbert(ButFilter(
%% compute the same for SWS
nT = size(rLfpAll,1);
In = WithinRanges([1:nT]',States.SWS);
myrLfpAll = rLfpAll(In,:);
sleep_dc1 = rdc100f(In,:);

fPow1 = [bsxfun(@plus,[1:0.5:15],[-0.3 0.3]')']/50;
fPh1 = bsxfun(@plus,[0.03:0.02:0.2], [-0.02 0.02]')'/50;
%[out4] = PowerPhasePairs(sleep_dc1,  fPh1, myrLfpAll , fPow1,1,'but',@PowerModulation);
% reduced to 1 ch AC and 1 DC 
Shuffle.MaxShift= 100*100/10; % after 10 times resampling will have 100sec shift.
Shuffle.Type = 'shift'; Shuffle.nShuffle = 5000;
[out4sh] = PowerPhasePairs(sleep_dc1(:,7),  fPh1, myrLfpAll(:,13) , fPow1,10,@PowerModulation, Shuffle);
[out4mud] = PowerPhasePairs(sleep_dc1(:,7),  fPh1, myrLfpAll(:,13) , fPow1,10,@PowerModulation);

%% 
%% plot maps
gch = setdiff([1:66],[1 56]); % all channels of 6x11 array map
bch = [ 1 56 57 35  46 47 48 66]; % bad channels in the array to gray out in the displa
for k=1:size(fPow1,1);
    figure(3324);clf
    dcy = rem([1:8],4); dcy(dcy==0)=4;
    dcx = double([1:8]>4)+1;
    for l=1:8
        if l==5 continue; end
        subplot2(4,12,dcy(l),dcx(l));% Ramp REM
        map(gch) = sq(mean(out2.Ramp(2:5,k,:,l)));
        map(bch)=nan;
        cscale = [0 max(reshape(mean(out2.Ramp(2:5,k,:,:)),[],1))];
        imagescnan(reshape(map,11,6),cscale,0,l==8);
        if l==1; title('REM Ramp'); end
       
        subplot2(4,12,dcy(l),dcx(l)+2); %Rrh REM
        map(gch) = sq(circmean(out2.Rth(2:5,k,:,l)));
        map(bch)=nan;
        imagescnan(reshape(map,11,6),[-pi pi],1,l==8);
        if l==1; title('REM Phase'); end
       
        subplot2(4,12,dcy(l),dcx(l)+4); %Ramp RUN
        map(gch) = sq(mean(out3.Ramp(2:5,k,:,l)));
        map(bch)=nan;
        cscale = [0 max(reshape(mean(out3.Ramp(2:5,k,:,:)),[],1))];
        imagescnan(reshape(map,11,6),cscale,0,l==8);
       if l==1; title('RUN Ramp'); end
       
         subplot2(4,12,dcy(l),dcx(l)+6); %Rth RUN
        map(gch) = sq(circmean(out3.Rth(2:5,k,:,l)));
        map(bch)=nan;
        imagescnan(reshape(map,11,6),[-pi pi],1,l==8);
        if l==1; title('RUN Phase'); end
       
        subplot2(4,12,dcy(l),dcx(l)+8); %Ramp SWS
        map(gch) = sq(mean(out4.Ramp(2:5,k,:,l)));
        map(bch)=nan;
        cscale = [0 max(reshape(mean(out4.Ramp(2:5,k,:,:)),[],1))];
        imagescnan(reshape(map,11,6),cscale,0,l==8);
       if l==1; title('SWS Ramp'); end
       
        
         subplot2(4,12,dcy(l),dcx(l)+10); %Rth SWS
        map(gch) = sq(circmean(out4.Rth(2:5,k,:,l)));
        map(bch)=nan;
        imagescnan(reshape(map,11,6),[-pi pi],1,l==8);
        if l==1; title('SWS Phase'); end
       
        
    end
    title(num2str(out2.fPow(k)*50));
    pause
end

%% display fPh x fPow 
myout = out2;
refch = 13;
figure
clf
for k=1:8
    subplot2(3,8,1,k);
    pcolor( myout.fPh(2:end)*50, myout.fPow*50, myout.MI(2:end,:,refch,k)');
   
   subplot2(3,8,2,k);
    pcolor( myout.fPh(2:end)*50, myout.fPow*50, myout.Ramp(2:end,:,refch,k)');
    
    subplot2(3,8,3,k);
    pcolor( myout.fPh(2:end)*50, myout.fPow*50, myout.Rth(2:end,:,refch,k)'); CircColormap; shading flat;
end

%% do this in time domain for specific bands 
%% 1. SWS, spindle 
In = WithinRanges([1:nT]',States.SWS);
myrLfpAll = rLfpAll(In,:);
sleep_dc1 = rdc100f(In,:);

% filter DC signal in ISA band
isadc = ButFilter(sleep_dc1,2,[0.07 0.2]/50, 'bandpass');
risadc = resample(isadc,1,10);
RefCh = 3;
lm = LocalMinima(isadc(:,RefCh), 100*2, 0)/10; % bring to 10hz sample rate
%lm(lm<=10*WinSec | lm>= size(lfp_power,1)-10*WinSec)=[];
% filter LFP in spindle band
fmyLfp = ButFilter(myrLfpAll,2,[10 13]/50,'bandpass');
lfp_power = resample(abs(hilbert(fmyLfp)),1,10); % at 10 hz sr.
% smooth to the ISA band
slfp_power = ButFilter(lfp_power,4,[0.015 0.2]/5,'bandpass');

WinSec = 20; %sec

DCSegs = GetSegs(risadc, lm - 10*WinSec, 2*WinSec*10+1);
tSegs = [-10*WinSec:10*WinSec]/10;
ISA_Amp = isadc(lm*10,RefCh);
OutlISA = ISA_Amp<prctile(ISA_Amp, 1);
% remove outlier triggers
lm = lm(~OutlISA);
ISA_Amp = ISA_Amp(~OutlISA);

PowSegs = GetSegs(slfp_power, lm - 10*WinSec, 2*WinSec*10+1);

[s si]= sort(ISA_Amp);
figure;
imagesc(tSegs,[], sq(zscore(PowSegs(:,si,4)))');

%mPowerSegs  = TriggeredAv(lfp_power, 10*WinSec, 10*WinSec, lm);
mPowSegs = sq(median(PowSegs,2));
wmPowSegs = sq(median(bsxfun(@times,PowSegs,ISA_Amp'),2))/mean(ISA_Amp);

zmPowSegs = zscore(mPowSegs); 
mDCSegs = sq(median(DCSegs,2));

%[~, PeakInd] = ; 

%% extract lags of trig PowSegs
for k=1:64
    for l=1:64
        [~,lag(k,l)] = max(xcorr(mPowSegs(:,k),mPowSegs(:,l)));
    end
end
lag = lag-400;
lag(55,:) = 0; lag(:,55)=0;

%%
figure(29); clf
subplot(221);
imagesc(tSegs, [], zmPowSegs'); 
subplot(223);
PlotTraces(mDCSegs(:,[1:4 6:8]),[-10*WinSec:10*WinSec]/10, 10, 3);

subplot(224);
imagesc( wmPowSegs'); 
% plot map of peak power
map(gch) = sq(zmPowSegs(218, :));
map(bch)=nan;
subplot(222);
imagescnan(reshape(map,11,6),[],0,1);


%% try coherence analysis on spindle power and SWS ISA DC
figure;
mtchd([risadc(:,3) lfp_power(:,4)], 2^12, 10, 2^11,[],2, 'linear',[],[0.015 0.5]);

%%

gT = [0.1e6:2e6];
mbch = [1 23 48 55 58];
mgch = setdiff([1:64], mbch);
fmyLfp1 = ButFilter(myrLfpAll(gT,mgch),2,[2 20]/50,'bandpass');
[ic gA gW] = fastica (fmyLfp1', 'lastEig', 20, 'numOfIC', 9);
A = nan(64,9);
A(mgch,:) = gA;

ic=ic';
figure; 
clf
for k=1:9
    subplot(3,3,k);
    sgnA = sign(nanmedian(A(:,k)));
    MapBrainCom64(sgnA*A(:,k));
end
[y,f,phi] = mtchd(ic, 2^8,100,2^7,[],1.5);
figure;
PlotMatrix(f,y);
ForAllSubplots('xlim([1 20])');

[out5] = PowerPhasePairs(sleep_dc1(gT,3),  fPh1, ic , fPow1,1,'but',@PowerModulation);

figure;
for k=1:9
    subplot(3,3,k);
    imagesc(out5.fPh(2:end)*50, out5.fPow*50, sq(out5.Ramp(2:end,:,k))');
end
figure
lfp_power(lfp_power(:)<0)=0;
[W,H ] = nmfsc(lfp_power', 5,0.2,0.01,'temp',0);
figure
for k=1:5
    subplot(3,2,k);
    MapBrainCom64(W(:,k));
end


%% REM and theta
nT = size(rLfpAll,1);
In = WithinRanges([1:nT]',States.REM) & ~WithinRanges([1:nT]',States.HVS);
myrLfpAll = rLfpAll(In,:);
sleep_dc1 = rdc100f(In,:);

% filter DC signal in ISA band
isadc = ButFilter(sleep_dc1,2,[0.07 0.2]/50, 'bandpass');
risadc = resample(isadc,1,10);
RefCh = 3;
lm = LocalMinima(isadc(:,RefCh), 100*2, 0)/10; % bring to 10hz sample rate
%lm(lm<=10*WinSec | lm>= size(lfp_power,1)-10*WinSec)=[];
% filter LFP in spindle band
fmyLfp = ButFilter(myrLfpAll,2,[7 10]/50,'bandpass');
lfp_power = resample(abs(hilbert(fmyLfp)),1,10); % at 10 hz sr.
% smooth to the ISA band
slfp_power = ButFilter(lfp_power,4,[0.015 0.2]/5,'bandpass');

WinSec = 20; %sec

DCSegs = GetSegs(risadc, lm - 10*WinSec, 2*WinSec*10+1);
tSegs = [-10*WinSec:10*WinSec]/10;
ISA_Amp = isadc(lm*10,RefCh);
OutlISA = ISA_Amp<prctile(ISA_Amp, 1);
% remove outlier triggers
lm = lm(~OutlISA);
ISA_Amp = ISA_Amp(~OutlISA);

PowSegs = GetSegs(slfp_power, lm - 10*WinSec, 2*WinSec*10+1);

[s si]= sort(ISA_Amp);
figure;
imagesc(tSegs,[], sq(zscore(PowSegs(:,si,4)))');

%mPowerSegs  = TriggeredAv(lfp_power, 10*WinSec, 10*WinSec, lm);
mPowSegs = sq(nanmedian(PowSegs,2));
wmPowSegs = sq(nanmedian(bsxfun(@times,PowSegs,ISA_Amp'),2))/mean(ISA_Amp);

zmPowSegs = zscore(mPowSegs); 
mDCSegs = sq(nanmedian(DCSegs,2));

%[~, PeakInd] = ; 

%% extract lags of trig PowSegs
for k=1:64
    for l=1:64
        [~,lag(k,l)] = max(xcorr(mPowSegs(:,k),mPowSegs(:,l)));
    end
end
lag = lag-400;
lag(55,:) = 0; lag(:,55)=0;

%%
figure(29); clf
subplot(221);
imagesc(tSegs, [], zmPowSegs'); 
subplot(223);
PlotTraces(mDCSegs(:,[1:4 6:8]),[-10*WinSec:10*WinSec]/10, 10, 3);

subplot(224);
imagesc( wmPowSegs'); 
% plot map of peak power
map(gch) = sq(zmPowSegs(250, :));
map(bch)=nan;
subplot(222);
imagescnan(reshape(map,11,6),[],0,1);



%% PARAFAC on Ramp 

mat = sq(mean(out2.Ramp(2:5,:,:,[1:4 6:8])));
rmat = reshape(mat, size(mat,1), []);
figure
imagesc(out2.fPow*50, [], rmat');
[LU, LR, FSr, VT] = erpPCA( rmat', 4);

[Factors,it,err,corcondia] = parafac(mat,3,[],2);
figure; 
for k=1:3
        subplot2(3,3,k,1);
        plot(out4.fPow*50, Factors{1}(:,k));
        
        subplot2(3,3,k,2);
        MapBrainCom64( Factors{2}(:,k));
        
        subplot2(3,3,k,3);
        bar(Factors{3}(:,k));
end


%% plot phase resolved power for select channels/ freqs.
figure
imagesc(out2.phbins(:,1,1,1,1),out2.fPow*50,sq(mean(out3.pow_dens(:,2:5,:, 13,7),2))'); colormap jet; axis xy

thfi = find(out2.fPow*50>7 & out2.fPow*50<10);
thfi1 = find(out2.fPow*50>5 & out2.fPow*50<7);
spfi = find(out2.fPow*50>10 & out2.fPow*50<13);
%%
figure(18);clf
nc=4; nr=6;

subplot2(nr, nc, 1, 1);
imagesc(out2.phbins(:,1,1,1,1),out2.fPow*50,sq(mean(out2.pow_dens(:,3,:, 13,7),2))'); colormap jet; axis xy
ylabel('LFP power frequency (Hz)');
title('REM'); 

subplot2(nr, nc, 1, 2);
imagesc(out2.phbins(:,1,1,1,1),out2.fPow*50,sq(mean(out3.pow_dens(:,4,:, 13,7),2))'); colormap jet; axis xy
title('RUN'); 

subplot2(nr, nc, 1, 3);
imagesc(out2.phbins(:,1,1,1,1),out2.fPow*50,sq(mean(out4.pow_dens(:,2:5,:, 13,2),2))'); colormap jet; axis xy
title('SWS'); 

subplot2(nr, nc, 1, 4);
imagesc(out2.phbins(:,1,1,1,1),out2.fPow*50,sq(mean(out2.pow_dens(:,3,:, 13,4),2))'); colormap jet; axis xy
ylabel('LFP power frequency (Hz)');
title('REM'); 


subplot2(nr, nc, 2, 1);
mm = sq(mean(out2.pow_dens(:,2:end,thfi,13,7),3));
imagesc(out2.phbins(:,1,1,1,1),out2.fPh*50,mm'); colormap jet; axis xy
ylabel('ISA frequency (Hz)');

subplot2(nr, nc, 2, 2);
mm = sq(mean(out3.pow_dens(:,2:end,thfi,13,7),3));
imagesc(out2.phbins(:,1,1,1,1),out2.fPh*50,mm'); colormap jet; axis xy

subplot2(nr, nc, 2, 3);
mm = sq(mean(out4.pow_dens(:,2:end,spfi,13,2),3));
imagesc(out2.phbins(:,1,1,1,1),out2.fPh*50,mm'); colormap jet; axis xy

subplot2(nr, nc, 2, 4);
mm = sq(mean(out2.pow_dens(:,2:end,thfi1,13,7),3));
imagesc(out2.phbins(:,1,1,1,1),out2.fPh*50,mm'); colormap jet; axis xy
ylabel('ISA frequency (Hz)');


subplot2(nr, nc, 3, 1);
imagesc(out2.phbins(:,1,1,1,1),[],sq(mean(mean(out2.pow_dens(:,2:5,thfi, 13,[1:4 6:8]),2),3))'); colormap jet; axis xy
xlabel('ISA phase, rad '); 
ylabel('DC channels');
subplot2(nr, nc, 3, 2);
imagesc(out2.phbins(:,1,1,1,1),[],sq(mean(mean(out3.pow_dens(:,2:5,thfi, 13,[1:4 6:8]),2),3))'); colormap jet; axis xy
xlabel('ISA phase, rad'); 
subplot2(nr, nc, 3, 3);
imagesc(out2.phbins(:,1,1,1,1),[],sq(mean(mean(out4.pow_dens(:,2:5,spfi, 13,[1:4 6:8]),2),3))'); colormap jet; axis xy
xlabel('ISA phase, rad'); 
subplot2(nr, nc, 3, 4);
imagesc(out2.phbins(:,1,1,1,1),[],sq(mean(mean(out2.pow_dens(:,2:5,thfi1, 13,[1:4 6:8]),2),3))'); colormap jet; axis xy
xlabel('ISA phase, rad '); 

subplot2(nr, nc, 4, 1);
map = sq(mean(mean(out2.Ramp(2:5,thfi, :,7),1),2));
MapBrainCom64(map);

subplot2(nr, nc, 4, 2);
map = sq(mean(mean(out3.Ramp(2:5,thfi, :,7),1),2));
MapBrainCom64(map);

subplot2(nr, nc, 4, 3);
map = sq(mean(mean(out4.Ramp(2:5,spfi, :,2),1),2));
MapBrainCom64(map);

subplot2(nr, nc, 4, 4);
map = sq(mean(mean(out2.Ramp(2:5,thfi1, :,7),1),2));
MapBrainCom64(map);


subplot2(nr, nc, 5, 1);
map = sq(mean(mean(out2.Ramp(2:5,thfi, 13,:),1),2)); map(5)=nan;
imagescnan(reshape(map,4,2),[],0,1);title('ISA modulation');

subplot2(nr, nc, 5, 2);
map = sq(mean(mean(out3.Ramp(2:5,thfi, 13,:),1),2)); map(5)=nan;
imagescnan(reshape(map,4,2),[],0,1);title('ISA modulation');

subplot2(nr, nc, 5, 3);
map = sq(mean(mean(out4.Ramp(2:5,spfi, 13,:),1),2)); map(5)=nan;
imagescnan(reshape(map,4,2),[],0,1);title('ISA modulation');

subplot2(nr, nc, 5, 4);
map = sq(mean(mean(out2.Ramp(2:5,thfi1, 13,:),1),2)); map(5)=nan;
imagescnan(reshape(map,4,2),[],0,1);title('ISA modulation');


subplot2(nr, nc, 6, 1);
map = sq(circmean(sq(circmean(out2.Rth(2:5,thfi, 13,:))))); map(5)=nan;
imagescnan(reshape(map,4,2),[-pi pi], 1, 1); title('ISA phase');

subplot2(nr, nc, 6, 2);
map = sq(circmean(sq(circmean(out3.Rth(2:5,thfi, 13,:))))); map(5)=nan;
imagescnan(reshape(map,4,2),[-pi pi], 1, 1); title('ISA phase');

subplot2(nr, nc, 6, 3);
map = sq(circmean(sq(circmean(out4.Rth(2:5,thfi, 13,:))))); map(5)=nan;
imagescnan(reshape(map,4,2),[-pi pi], 1, 1); title('ISA phase');

subplot2(nr, nc, 6, 4);
map = sq(circmean(sq(circmean(out2.Rth(2:5,thfi1, 13,:))))); map(5)=nan;
imagescnan(reshape(map,4,2),[-pi pi], 1, 1); title('ISA phase');


%% reduced version of the figure
% prepare data
spec(2).State = 'REM';
spec(2).In = WithinRanges([1:nT]',States.REM) & ~WithinRanges([1:nT]',States.HVS);
spec(1).State = 'RUN';
spec(1).In = WithinRanges([1:nT]',States.AWKTH);
spec(3).State = 'SWS';
spec(3).In = WithinRanges([1:nT]',States.SWS);

for k=1:3
spec(k).mylfp = rLfpAll(spec(k).In,13);
spec(k).mysleep_dc = rdc100f(spec(k).In,7);
spec(k).myisadc = ButFilter(spec(k).mysleep_dc,2,[0.07 0.2]/50, 'bandpass');
[spec(k).y, spec(k).f, spec(k).t] = mtchglong(WhitenSignal(spec(k).mylfp),2^8,100, 2^7,2^7-2^5,2.5,'linear',[],[1 20]);
spec(k).t = spec(k).t+2^6/100;
end


%% plot figure
k=3;
figure(333);clf
subplot(211);
imagesc(spec(k).t, spec(k).f, log10(spec(k).y)'); axis xy; colormap jet
subplot(212);
plot([1:size(spec(k).myisadc,1)]/100, spec(k).myisadc); axis tight
linkx

%%
spec(k).per = xlim;
spec(3).ca = [0 4.5]; spec(1).ca=[]; spec(2).ca=[];

%%collect out's 
outs = struct([]);
outs =CopyStruct(out3, outs,1); 
outs =CopyStruct(out2, outs,1); 
outs =CopyStruct(out4, outs,1); 
outs_sh = struct([]);
outs_sh = CopyStruct(out2sh, outs_sh,1); 
outs_sh = CopyStruct(out4sh, outs_sh,1); 

spec(1).fi = thfi; 
spec(2).fi = thfi; 
spec(3).fi = spfi; 



%% final new figure on ISA-LFP 
% only for REM and SWS
figure(888); clf
set(gcf, 'PaperPositionMode', 'auto');
nc = 2; nr = 3;
for l=1:2
    k=l+1;
    %subplot2(nr, nc, 1, k);

    [h, crit_p, adj_p]=fdr_bh(outs_sh(l).Rpval,0.0001);
    plot_mat = outs_sh(l).Ramp;

   % plot_mat = outs_sh(l).Rzscore;
    plot_mat(adj_p>0.0001)=nan;
    
    
%     imagesc(spec(k).t, spec(k).f, log10(spec(k).y)'); axis xy; colormap jet
%     xlim(spec(k).per); if ~isempty(spec(k).ca), caxis(spec(k).ca); end
%     title(spec(k).State); if k==1 ylabel('Frequency'); end
%     
%     subplot2(nr, nc, 2, k);
%     plot([1:size(spec(k).myisadc,1)]/100, spec(k).myisadc); axis tight
%     xlim(spec(k).per); box off ; axis off; legend('hide');
%     
    subplot2(nr, nc, 1, l);
    imagescnan({outs(k).fPh(2:end)*50, outs(k).fPow*50, plot_mat'},[],0,1); axis xy
    xlabel('ISA phase frequency'); if l==1 ylabel('LFP power frequency'); end
    title(spec(k).State);
    
    subplot2(nr, nc, 3, l);
 %   imagesc(outs(k).phbins(:,1,1,1,1),outs(k).fPow*50,sq(mean(outs(k).pow_dens(:,2:5,:, 13,7),2))'); colormap jet; axis xy; colorbar
    polarplot(outs(k).phbins(:,1,1,1,1),sq(mean(mean(outs(k).pow_dens(:,2:5,spec(k).fi, 13,7),2),3))'); axis tight; box off
    %hold on;
    %title(spec(k).State);
    
 
    %xlabel('ISA phase, rad'); if k==1 ylabel('LFP power , Hz'); end
        
    subplot2(nr, nc, 2, l);
    map = sq(mean(mean(outs(k).Ramp(2:5,spec(k).fi, :,2),1),2));
    MapBrainCom64(map);
    xlabel('ML'); ylabel('AP');
end

%AxisPropForIllustrator(14);
%print -depsc2 -r300 fig_isa_lfp.eps


%% examples fig
figure(3433);clf
nc = 2; nr = 2;
for l=1:2
    k=l+1;
    subplot2(nr, nc, 1, l);
    imagesc(spec(k).t, spec(k).f, log10(spec(k).y)'); axis xy; colormap jet
    xlim(spec(k).per); if ~isempty(spec(k).ca), caxis(spec(k).ca); end
    title(spec(k).State); if k==1 ylabel('Frequency'); end
    
    subplot2(nr, nc, 2, l);
    plot([1:size(spec(k).myisadc,1)]/100, spec(k).myisadc); axis tight
    xlim(spec(k).per); box off ; axis off; legend('hide');
end

 
%% try high frequency
myLfpSel = [];
for fn = 1:5
    FileName = d(fn).name;
    Lfp = LoadBinary(FileName, RepCh(13), Par.nChannels,2)';
    myLfpSel= [myLfpSel; Lfp];
    
end
dcsel = [];
for fn = 1:5
    FileName = d(fn).name;
    dc = LoadBinary(FileName, dc_Chan(7), Par.nChannels,2)';
    dcsel = [dcsel; dc];
end
fdcsel = ButFilter(dcsel,4,0.1/500,'low');


nT = size(myLfpSel,1);
In = WithinRanges([1:nT]',States.REM*10) & ~WithinRanges([1:nT]',States.HVS*10); % NOW IN 1KHZ
myLfp = myLfpSel(In,:);
mydc = fdcsel(In,:);
fPowGam = bsxfun(@plus, [15:5:200],[-10 10]')'/500; % NB!!! this has to be adjusted to the signal sample rate - 1kHz. 
fPhGam = bsxfun(@plus,[0.04:0.01:0.5], [-0.03 0.03]')'/500;
Shuffle.MaxShift= 1000*100/100; % after 10 times resampling will have 50sec shift.
Shuffle.Type = 'shift'; Shuffle.nShuffle = 1000;
[outgam] = PowerPhasePairs(mydc,  fPhGam, myLfp , fPowGam,100,@PowerModulation,Shuffle);
outgam.fPow = mean(fPowGam,2)*500;
outgam.fPh = mean(fPhGam,2)*500;

%%

%%% THETA modulation of GAMMA
%%same for theta 
fPhTh = [bsxfun(@plus,[1:0.5:15],[-0.3 0.3]')']/500;
fPowGam = bsxfun(@plus, [30:5:200],[-10 10]')'/500; % NB!!! this has to be adjusted to the signal sample rate - 1kHz. 
Shuffle.MaxShift= 1000*5/10; % after 10 times resampling will have 50sec shift.
Shuffle.Type = 'shift'; Shuffle.nShuffle = 1000;
[outthgam] = PowerPhasePairs(myLfp,  fPhTh, myLfp , fPowGam,10,@PowerModulation, Shuffle);
outthgam.fPow = mean(fPowGam,2)*500;
outthgam.fPh = mean(fPhTh,2)*500;


for k=1:64
    myLfpSel = [];
    for fn = 1:5
        FileName = d(fn).name;
        Lfp = LoadBinary(FileName, RepCh(k), Par.nChannels)';
        myLfpSel= [myLfpSel; Lfp];
        
    end
    myLfp = myLfpSel(In,:);
    [outthgamall{k}] = PowerPhasePairs(myLfp,  fPhTh, myLfp , fPowGam,10,@PowerModulation);
    [outgamall{k}] = PowerPhasePairs(mydc,  fPhGam, myLfp , fPowGam,100,@PowerModulation);
end

for k=1:64
    tmp(k)=outgamall{k};
   
end
outgamallnew= CatStruct(tmp); 
    
%% plot thgam modulation maps

figure(7);clf
subplot2(2,4,1,1);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgamallnew.Ramp(:,4:end,13)'); axis xy; colorbar; xlim([1 12]);
subplot2(2,4,2,1);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgamallnew.Rth(:,4:end,13)'); axis xy; CircColormap; colorbar; xlim([1 12]);
subplot2(2,4,1,2);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgamallnew.Ramp(:,4:end,32)'); axis xy; colorbar;xlim([1 12]);
subplot2(2,4,2,2);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgamallnew.Rth(:,4:end,32)'); axis xy; CircColormap; colorbar;xlim([1 12]);

powfi = find(outthgam.fPow>120 & outthgam.fPow<160);
phfi = find(outthgam.fPh>6.5 & outthgam.fPh<7.5);
mR = sq(mean(outthgamallnew.Ramp(phfi,powfi,:)));
subplot2(2,4,1,3);
MapBrainCom64(mR);

mR = sq(circmean(outthgamallnew.Rth(phfi,powfi,:)));
subplot2(2,4,2,3);
MapBrainCom64(mR(:),[-pi pi],1,1);

powfi = find(outthgam.fPow>70 & outthgam.fPow<110);
phfi = find(outthgam.fPh>3.5 & outthgam.fPh<5.5);
mR = sq(mean(mean(outthgamallnew.Ramp(phfi,powfi,:)),2));
subplot2(2,4,1,4);
MapBrainCom64(mR);

mR = sq(circmean(sq(circmean(outthgamallnew.Rth(phfi,powfi,:)))));
subplot2(2,4,2,4);
MapBrainCom64(mR(:),[-pi pi],1,1);


% %% do parafac on thgam
% %addpath /Users/antsiro/Work/Code/matlab//matlab/toolboxes/nway300
% nfac = 2;
% fac = parafac(outthgamallnew.Ramp.*exp(sqrt(-1)*outthgamallnew.Rth),nfac);
% phfac = cellfun(@angle,fac,'UniformOutput',false);
% fac = cellfun(@abs,fac,'UniformOutput',false);

% 
%  figure(1); clf
% for k=1:nfac
%     subplot2(3,nfac,1,k);
%     plot(outthgam.fPh, fac{1}(:,k)); xlabel('fPh'); axis tight
%     title(['factor # ' num2str(k)]);
%     subplot2(3,nfac,2,k);
%     plot(outthgam.fPow, fac{2}(:,k)); xlabel('fPow');axis tight
%     subplot2(3,nfac,3,k);
%     MapBrainCom64(fac{3}(:,k));
% end
    
%%
figure(6);clf
subplot(221);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgam.Ramp(:,4:end)'); axis xy; colorbar
subplot(222);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgam.Rth(:,4:end)'); axis xy; CircColormap; colorbar
subplot(223);
imagesc(outthgam.fPh, outthgam.fPow(4:end), outthgam.MI(:,4:end)'); axis xy;  colorbar



%%
figure(4);clf
subplot(221);
[h, crit_p, adj_p]=fdr_bh(outgam.Rpval,0.0001);
plot_mat = outgam.Ramp;
plot_mat(adj_p>0.0001)=nan;
imagescnan({outgam.fPh(3:end), outgam.fPow(4:end), plot_mat(3:end,4:end)'}); axis xy; xlim([0.06 0.2])
%imagesc(outgam.fPh(3:end), outgam.fPow(4:end), outgam.Ramp(3:end,4:end)'); axis xy
subplot(222);
imagesc(outgam.fPh(3:end), outgam.fPow(4:end), outgam.Rth(3:end,4:end)'); axis xy; CircColormap; xlim([0.06 0.2])
subplot(223);
[h, crit_p, adj_p]=fdr_bh(outgam.MIpval,0.0001);
plot_mat = outgam.MI;
plot_mat(adj_p>0.0001)=nan;
imagescnan({outgam.fPh(3:end), outgam.fPow(4:end), plot_mat(3:end,4:end)'}); axis xy; xlim([0.06 0.2])
%imagesc(outgam.fPh(3:end), outgam.fPow(4:end), outgam.MI(3:end,4:end)'); axis xy
subplot(224);
gamfi = find(outgam.fPow>70 & outgam.fPow<120);
polarplot(outgam.phbins(:,1,1),sq(mean(mean(outs(k).pow_dens(:,1:7,gamfi),2),3))'); axis tight; box off
   
figure; plot(outgam.phbins(:,1,1),sq(mean(mean(outs(k).pow_dens(:,1:7,gamfi),2),3))')
phdens = sq(mean(mean(outs(k).pow_dens(:,1:7,gamfi),2),3))';
r = sum(phdens'.*exp(sqrt(-1)*outgam.phbins(:,1,1)))./sum(phdens);


%%
% 
% myLfpSel = [];
% for fn = 1:5
%     FileName = d(fn).name;
%     Lfp = LoadBinary(FileName, RepCh(13), Par.nChannels,2)';
%     myLfpSel= [myLfpSel; Lfp];
%     
% end
    %%%% Figure 5h
    figure(10);clf
    subplot(221);
    [h, crit_p, adj_p]=fdr_bh(outgam.Rpval,0.0001);
    plot_mat = outgam.Ramp;
    plot_mat(adj_p>0.0001)=nan;
    imagescnan({outgam.fPh(3:end), outgam.fPow(4:end), plot_mat(3:end,4:end)'},[],0,1); axis xy; xlim([0.06 0.2])
    xlabel('ISA phase frequency'); ylabel('LFP power frequency');

    %imagesc(outgam.fPh(3:end), outgam.fPow(4:end), outgam.Ramp(3:end,4:end)'); axis xy
    subplot(223);
    [h, crit_p, adj_p]=fdr_bh(outthgam.Rpval,0.001);
    plot_mat = outthgam.Ramp;
    plot_mat(adj_p>0.001)=nan;
    imagescnan({outthgam.fPh(3:end), outthgam.fPow(4:end), plot_mat(3:end,4:end)'},[],0,1); axis xy; xlim([2 12]);
    %imagesc(outthgam.fPh(3:end), outthgam.fPow(4:end), outthgam.Ramp(3:end,4:end)'); axis xy
    xlabel('LFP phase frequency'); ylabel('LFP power frequency');
    subplot(222)
    polarplot(outgam.phbins(:,1,1),sq(mean(mean(outgam.pow_dens(:,1:7,gamfi),2),3))'); axis tight; box off

    subplot(224)
    powfi = find(outthgam.fPow>120 & outthgam.fPow<160);
    phfi = find(outthgam.fPh>6.5 & outthgam.fPh<7.5);
    polarplot(outthgam.phbins(:,1,1),sq(mean(mean(outthgam.pow_dens(:,phfi,powfi),2),3))'); axis tight; box off
