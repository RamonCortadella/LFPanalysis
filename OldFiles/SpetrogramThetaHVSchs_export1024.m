close('all')
%% Define settings

directory = '../../../data/ASIC1024/B14062W18-T1-rat01601/DatData/';
Par =  LoadXml(strcat(directory,'Rec_map_SIG_B14062W18-T1_54_S-ACinterAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
d = dir(strcat(directory,'Rec_map_SIG_B14062W18-T1_54_S-ACinterAC.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');

LFPfs = 250;
DownSampDC = 100;

timeRes = 1;
timeResDC = 100;

Fs=LFPfs;
FsDC = Fs/DownSampDC;

% Parameters PSD
SampsRes = timeRes * Fs;
SampsResDC = timeResDC * FsDC;
window = floor(log2(SampsRes));
windowDC = floor(log2(SampsResDC));
nperov = 2^(window-4);
nperovDC = 2^(windowDC-4);
nfft = 2^window;
nfftDC = 2^windowDC;
nCh = 1024;
nChCol = 32;

% Parameters spectrogram
FactTheta = 6;
FactHVS = 2;
nFFTtheta = 512*FactTheta;
nFFTHVS = 512*FactHVS;
Fmax = 50;
FmaxDC = 50;

ChIndHVS= 229; % centre of the 9 channels averaged
ChIndTheta = 93;

StartStop = false;
Tinit =100; %initial and final times of the spectrogram 
Tend = 6500;

%% get ephys data
ACLfp = [];
DCLfp = [];
for fn = 1:length(d)
    FileName = [d(fn).folder ,'/',d(fn).name];
    Lfp = LoadBinaryDAT(FileName, [0:nCh-1], Par.nChannels,1)';
    
    SplitName = split(d(fn).name,'-');
    
%     if SplitName{3}(end-5:end-4)=='AC'
%     myrLfp = resample(Lfp,1,10); % RESAMPLE TO 100 HZ !!!
    LfpGeom = Lfp;
%     elseif  SplitName{3}(end-5:end-4)=='DC' 
%         LfpGeomDC = Lfp;
%     end
end

%% Compute PSD DC/AC channels 

[wx,armodel] = WhitenSignal(LfpGeom,[],[],[],1);
% [wxDC,armodel] = WhitenSignal(LfpGeomDC,[],[],[],1);

figure(1)
subplot(2,2,1)
[ppxAC, fAC] = pwelch(LfpGeom,nfft,nperov,nfft,Fs);
loglog(fAC,ppxAC)
title('AC-channels')
subplot(2,2,2)
[ppx, fAC] = pwelch(wx,nfft,nperov,nfft,Fs);
loglog(fAC,ppx)
title('whitened AC-channels')
% subplot(2,2,3)
% [ppxDC, fDC] = pwelch(LfpGeomDC(100*FsDC:end,:),nfftDC,nperovDC,nfftDC,FsDC);
% loglog(fDC,ppxDC)
% title('DC-channels')
% subplot(2,2,4)
% [ppx, fDC] = pwelch(wxDC(100*FsDC:end,:),nfftDC,nperovDC,nfftDC,FsDC);
% loglog(fDC,ppx)
% title('whitened DC-channels')

%% plot power map at specific freq

ChIndThetaList= [ChIndTheta-1-nChCol:ChIndTheta+1-nChCol, ChIndTheta-1:ChIndTheta+1, ChIndTheta-1+nChCol:ChIndTheta+1+nChCol];
ChIndHVSList = [ChIndHVS-1-nChCol:ChIndHVS+1-nChCol, ChIndHVS-1:ChIndHVS+1, ChIndHVS-1+nChCol:ChIndHVS+1+nChCol];

FreqInterest = 8;

[bX, bY]= meshgrid([1:nChCol]',[1:nChCol]');
Coord(:,1) = reshape(repmat([1:nChCol],nChCol,1),[],1);
Coord(:,2)= repmat([1:nChCol]',nChCol,1);

figure()
[ppxAC, fAC] = pwelch(LfpGeom,nfft,nperov,nfft,Fs);
Findex = find(abs(fAC-FreqInterest)==min(abs(fAC-FreqInterest)));
% 
ppxAC(Findex,ChIndThetaList)=ppxAC(Findex,ChIndThetaList)*100;
map = ppxAC(Findex,:);
F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
bF = F(bX,bY);
imagesc(bF);
display(fAC(Findex))
% % caxis([2e4 5e5])% 7Hz caxis([1e4 1e6]) 4Hz caxis([1e3 1.5e5]) 1Hz caxis([1e5 2e6])

title(strcat('Power @',string(FreqInterest),'Hz AC-channels'))

% subplot(1,2,2)
% map = ppxDC(Findex,:);
% F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
% bF = F(bX,bY);
% imagesc(bF);
% display(fDC(Findex))
% % caxis([2e2 1e4])%7Hz caxis([1e2 4e3]), 4Hz caxis([1e1 4e2]), 1Hz caxis([1e3 2e4]), 0.5Hz caxis([1e2 4e3]), 0.3Hz caxis([2e2 1e4])
% title(strcat('Power @',string(FreqInterest),'Hz DC-channels'))


%% config spectrogram (mtcsglong)

if StartStop == true
    Sinit=floor(Tinit*Fs);
    Sfin = floor(Tend*Fs);
%     SinitDC = floor(Tinit*FsDC);
%     SfinDC = floor(Tend*FDCs);
else
    Sinit = 1;
%     SinitDC = 1;
    Sfin = length(LfpGeom(:,1));
%     SfinDC = length(LfpGeomDC(:,1));
end

NFmaxTheta = floor(Fmax*nFFTtheta/Fs);
NFmaxHVS = floor(FmaxDC*nFFTHVS/Fs);


%% compute spectrogram whitened and mean-of-sveral-channels for AC theta and HVS channels
[ws,wf,wt] = mtcsglong(mean(wx(:,ChIndHVSList),2),nFFTHVS,Fs);%,1024);
[s,f,t] = mtcsglong(mean(LfpGeom(:,ChIndHVSList),2),nFFTHVS,Fs);%,1024);

[wst,wft,wtt] = mtcsglong(mean(wx(:,ChIndThetaList),2),nFFTtheta,Fs);%,1024);
[st,ft,tt] = mtcsglong(mean(LfpGeom(:,ChIndThetaList),2),nFFTtheta,Fs);%,1024);

%% save spectrograms and raw signals
WhitenedSpecHVS = struct('ws',ws,'wf',wf,'wt',wt);
RawSpecHVS = struct('s',s,'f',f,'t',t);
SpecHVS = struct('WhitenedSpecHVS',WhitenedSpecHVS,'RawSpecHVS',RawSpecHVS);
% save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-SpecHVS.mat'),'SpecHVS')


WhitenedSpecStates = struct('ws',wst,'wf',wft,'wt',wtt);
RawSpecStates = struct('s',st,'f',ft,'t',tt);
SpecStates = struct('WhitenedSpecStates',WhitenedSpecStates,'RawSpecStates',RawSpecStates);
% save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-SpecStates.mat'),'SpecStates')

%save DC channel Theta and HVS detection
% LfpDC_Chs = LfpGeomDC(:,ChIndThetaList);
% LfpDC_HVSChs = LfpGeomDC(:,ChIndHVSList);
% LfpAC_Chs = LfpGeom(:,ChIndThetaList);
% LfpAC_HVSChs = LfpGeom(:,ChIndHVSList);
% LFP_traces = struct('LfpDC_Chs',LfpDC_Chs,'LfpAC_Chs',LfpAC_Chs);
% save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-Traces.mat'),'LFP_traces')

%% plot spectrograms and AC/DC signals
figure()
%plot AC channel spectrogram
ax1 = subplot(3,1,1);
h= pcolor(t,f,log10(abs(s')));

ylim([Fs*FactTheta/nFFTtheta Fmax])
caxis([min(min(log10(abs(s'))))+1.9 max(max(log10(abs(s'))))-0.3])
set(gca,'YScale','log')
set(gca,'YDir','normal')
colormap('jet')
set(h, 'EdgeColor', 'none')

ax2 = subplot(3,1,2); %plot sig AC
plot((Sinit+(1:1:length(LfpGeom(Sinit:Sfin,1))))/LFPfs, LfpGeom(Sinit:Sfin,ChIndHVS))

% ax3 = subplot(3,1,3); %plot sigDC
% plot((SinitDC+(1:1:length(LfpGeomDC(SinitDC:SfinDC,1))))/FsDC, LfpGeomDC(SinitDC:SfinDC,ChIndHVS))

set(h, 'EdgeColor', 'none')
linkaxes([ax1,ax2,ax3],'x')
xlim([Sinit/Fs Sfin/Fs])


