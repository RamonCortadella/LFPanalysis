close('all')
directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec1interAC.xml'));
d = dir(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec1inter*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');

LFPfs = 651.04166667;
DownSampDC = 100;
%% get ephys data
ACLfp = [];
DCLfp = [];
for fn = 1:length(d)
    FileName = [d(fn).folder ,'/',d(fn).name];
    Lfp = LoadBinaryDAT(FileName, [0:255], Par.nChannels,1)';
    
    SplitName = split(d(fn).name,'-');
    
    if SplitName{3}(end-5:end-4)=='AC'
%     myrLfp = resample(Lfp,1,10); % RESAMPLE TO 100 HZ !!!
        LfpGeom = Lfp;
    elseif  SplitName{3}(end-5:end-4)=='DC' 
        LfpGeomDC = Lfp;
    end
end

%%

% for i = 1:length(Par.AnatGrps)
%     for ii = 1:length(Par.AnatGrps(i).Channels)
%         LfpGeom(:,ii+(i-1)*16) = ACLfp(:,Par.AnatGrps(i).Channels(ii)+1);
%         LfpGeomDC(:,ii+(i-1)*16) = DCLfp(:,Par.AnatGrps(i).Channels(ii)+1);
%     end
% end
% 
% figure()
% loglog(pwelch(LfpGeom(:,1)',2^14,2^10,2^14,Fs))


%%
close('all')
timeRes = 10;
timeResDC = 100;

Fs=LFPfs;
FsDC = Fs/DownSampDC;
SampsRes = timeRes * Fs;
SampsResDC = timeResDC * FsDC;

window = floor(log2(SampsRes));
windowDC = floor(log2(SampsResDC));
nperov = 2^(window-4);
nperovDC = 2^(windowDC-4);
nfft = 2^window;
nfftDC = 2^windowDC;

% 

[wx,armodel] = WhitenSignal(LfpGeom,[],[],[],1);
[wxDC,armodel] = WhitenSignal(LfpGeomDC,[],[],[],1);

figure(1)
subplot(2,2,1)
[ppxAC, fAC] = pwelch(LfpGeom,nfft,nperov,nfft,Fs);
loglog(fAC,ppxAC)
title('AC-channels')
subplot(2,2,2)
[ppx, fAC] = pwelch(wx,nfft,nperov,nfft,Fs);
loglog(fAC,ppx)
title('whitened AC-channels')
subplot(2,2,3)
[ppxDC, fDC] = pwelch(LfpGeomDC(100*FsDC:end,:),nfftDC,nperovDC,nfftDC,FsDC);
loglog(fDC,ppxDC)
title('DC-channels')
subplot(2,2,4)
[ppx, fDC] = pwelch(wxDC(100*FsDC:end,:),nfftDC,nperovDC,nfftDC,FsDC);
loglog(fDC,ppx)
title('whitened DC-channels')
% figure(2)
% [ppx, f] = pwelch(zscore(LfpGeom')',2^14,2^10,2^14,Fs);
% loglog(f,ppx)
% figure(4)
% wx2 = whitening(LfpGeom,Fs,'freq',[1 48]);
% [ppx, f] = pwelch(wx2,2^14,2^10,2^14,Fs);
% loglog(f,ppx)
% 
% 
%% plot power map at specific freq

ChIndHVS= 229;
ChIndTheta = 93;
% ChIndHVS = 18;

ChIndThetaList= [ChIndTheta-1-16:ChIndTheta+1-16, ChIndTheta-1:ChIndTheta+1, ChIndTheta-1+16:ChIndTheta+1+16];
ChIndHVSList = [ChIndHVS-1-16:ChIndHVS+1-16, ChIndHVS-1:ChIndHVS+1, ChIndHVS-1+16:ChIndHVS+1+16];

FreqInterest = 8;

[bX, bY]= meshgrid([1:16]',[1:16]');
Coord(:,1) = reshape(repmat([1:16],16,1),[],1);
Coord(:,2)= repmat([1:16]',16,1);

figure()
% title('Power of ~1Hz oscillation')
% subplot(1,2,1)
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
% 
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
FactTheta = 6;
FactHVS = 2;
nFFTtheta = 512*FactTheta;
nFFTHVS = 512*FactHVS;
Fmax = 50;
FmaxDC = 50;

StartStop = false;
Tinit =100;
Tend = 6500;

if StartStop == true
    Sinit=floor(Tinit*Fs);
    Sfin = floor(Tend*Fs);
    SinitDC = floor(Tinit*FsDC);
    SfinDC = floor(Tend*FDCs);
else
    Sinit = 1;
    SinitDC = 1;
    Sfin = length(LfpGeom(:,1));
    SfinDC = length(LfpGeomDC(:,1));
end

NFmaxTheta = floor(Fmax*nFFTtheta/Fs);
NFmaxHVS = floor(FmaxDC*nFFTHVS/Fs);

%% save spectrograms and raw signals
% compute spectrogram whitened and mean-of-sveral-channels for AC theta ch
[ws,wf,wt] = mtcsglong(mean(wx(:,ChIndHVSList),2),nFFTHVS,Fs);%,1024);
[s,f,t] = mtcsglong(mean(LfpGeom(:,ChIndHVSList),2),nFFTHVS,Fs);%,1024);

WhitenedSpecHVS = struct('ws',ws,'wf',wf,'wt',wt);
RawSpecHVS = struct('s',s,'f',f,'t',t);
SpecHVS = struct('WhitenedSpecHVS',WhitenedSpecHVS,'RawSpecHVS',RawSpecHVS);
save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-SpecHVS.mat'),'SpecHVS')

[ws,wf,wt] = mtcsglong(mean(wx(:,ChIndThetaList),2),nFFTtheta,Fs);%,1024);
[s,f,t] = mtcsglong(mean(LfpGeom(:,ChIndThetaList),2),nFFTtheta,Fs);%,1024);

WhitenedSpecStates = struct('ws',ws,'wf',wf,'wt',wt);
RawSpecStates = struct('s',s,'f',f,'t',t);
SpecStates = struct('WhitenedSpecStates',WhitenedSpecStates,'RawSpecStates',RawSpecStates);
save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-SpecStates.mat'),'SpecStates')

% compute spectrogram whitened and mean-of-sveral-channels for AC HVS ch
% [ws,wf,wt] = mtcsglong(wx(:,ChIndHVS),nFFT,Fs);%,1024);
% [s,f,t] = mtcsglong(mean(LfpGeom(:,ChIndHVSList),2),nFFT,Fs);%,1024);
% 
% WhitenedSpecAC = struct('ws',ws,'wf',wf,'wt',wt);
% RawSpecAC = struct('s',s,'f',f,'t',t);
% SpecAC_HVS = struct('WhitenedSpecAC',WhitenedSpecAC,'RawSpecAC',RawSpecAC);
% save(strcat('../MatlabData/',SplitName{3}(1:end-6),'-SpecAC_ThetaCh.mat'),'SpecAC_HVS')

%save DC channel Theta and HVS detection
LfpDC_Chs = LfpGeomDC(:,ChIndThetaList);
% LfpDC_HVSChs = LfpGeomDC(:,ChIndHVSList);
LfpAC_Chs = LfpGeom(:,ChIndThetaList);
% LfpAC_HVSChs = LfpGeom(:,ChIndHVSList);
LFP_traces = struct('LfpDC_Chs',LfpDC_Chs,'LfpAC_Chs',LfpAC_Chs);
save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-Traces.mat'),'LFP_traces')

%% plot spectrograms and AC/DC signals
figure()
%plot AC channel spectrogram
ax1 = subplot(3,1,1);
h= pcolor(t,f(1:NFmaxTheta),log10(abs(s(:,1:NFmaxTheta)')));

ylim([Fs*FactTheta/nFFTtheta Fmax])
caxis([min(min(log10(abs(s(:,1:NFmaxTheta)'))))+1.9 max(max(log10(abs(s(:,1:NFmaxTheta)'))))-0.3])
set(gca,'YScale','log')
set(gca,'YDir','normal')
colormap('jet')
set(h, 'EdgeColor', 'none')

ax2 = subplot(3,1,2); %plot sig AC
plot((Sinit+(1:1:length(LfpGeom(Sinit:Sfin,1))))/LFPfs, LfpGeom(Sinit:Sfin,ChIndHVS))


%plot DC channel spectrogram
% wxDChpf = ButFilter(wxDC,2,[0.01 50]/Fs*2,'bandpass');

% [s,f,t] = mtcsglong(wxDChpf(:,ChInd),nFFTDC,Fs);
% save(strcat('../MatlabData/',SplitName{3},'-WhitenedSpectDC.mat'),'s','f','t'

% [ws,wf,wt] = mtcsglong(mean(wxDChpf(:,[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]),2),nFFT,Fs);%,1024);
% [ppxDC, fDC] = pwelch(mean(wx(:,[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]),2),nFFTDC*2,nFFTDC,nFFTDC*2,Fs);%[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]
% figure()
% subplot(2,1,1)
% loglog(fDC,ppxDC)
% title('Whitened DC-channel exported')

% [s,f,t] = mtcsglong(mean(LfpGeomDC(:,[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]),2),nFFTDC,Fs);%,1024);
% [ppxDC, fDC] = pwelch(mean(LfpGeomDC(:,[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]),2),2^14,2^10,2^14,Fs);%[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]
% subplot(2,1,2)
% loglog(fDC,ppxDC)
% title('DC-channel exported')

% WhitenedSpecDC = struct('ws',ws,'wf',wf,'wt',wt);
% RawSpecDC = struct('s',s,'f',f,'t',t);
% SpecDC = struct('WhitenedSpecDC',WhitenedSpecDC,'RawSpecDC',RawSpecDC);
% save(strcat('../MatlabData/',SplitName{3}(1:end-6),'-SpecDC.mat'),'SpecDC')

% figure(3)
% ax3 = subplot(4,1,3); %plot specDC
% h= pcolor(t,f(1:NFmaxDC),log10(abs(s(:,1:NFmaxDC)')));

% ylim([Fs/nFFTDC FmaxDC])
% caxis([min(min(log10(abs(s(:,1:NFmaxDC)')))) max(max(log10(abs(s(:,1:NFmax)'))))])
% set(gca,'YScale','log')
% set(gca,'YDir','normal')
% colormap('jet')

ax3 = subplot(3,1,3); %plot sigDC
plot((SinitDC+(1:1:length(LfpGeomDC(SinitDC:SfinDC,1))))/FsDC, LfpGeomDC(SinitDC:SfinDC,ChIndHVS))

set(h, 'EdgeColor', 'none')
linkaxes([ax1,ax2,ax3],'x')
xlim([Sinit/Fs Sfin/Fs])

%% plot trigger Mocap
% figure()
% hold on
% Tstab = 10; %Time removed from the beggining to remove artifactual peaks not related to trigger
% TmaxInit = 17; %End of time where trigger is searched
% TmaxEnd = 2500;
% 
% time = [1:1:length(LfpGeomDC(:,1))]/Fs;
% plot(time(Tstab*floor(LFPfs):TmaxInit*floor(LFPfs)), LfpGeomDC(Tstab*floor(LFPfs):TmaxInit*floor(LFPfs),81))
% plot(time(end-TmaxEnd*floor(LFPfs):end), LfpGeomDC(end-TmaxEnd*floor(LFPfs):end,81))
% hold off
% 
% SigStart = LfpGeomDC(Tstab*floor(LFPfs):TmaxInit*floor(LFPfs),81);
% 
% 
% indexTrigger = find(abs(SigStart)==max(abs(SigStart)));
% indexTrigger = indexTrigger(1);
% 
% MocapStart = time(Tstab*floor(LFPfs)+indexTrigger);
% 
% display(LfpGeomDC(Tstab*floor(LFPfs)+indexTrigger-1,81))
% 
% indexTrigger = find(abs(LfpGeomDC(end-TmaxEnd*floor(LFPfs):end,81))==max(abs(LfpGeomDC(end-TmaxEnd*floor(LFPfs):end,81))));
% indexTrigger = indexTrigger(1);
% 
% MocapEnd = time(indexTrigger-1+length(LfpGeomDC(:,81))-TmaxEnd*floor(LFPfs));
% 
% display(LfpGeomDC(indexTrigger-1+length(LfpGeomDC(:,81))-TmaxEnd*floor(LFPfs),81))
% 
% %% save Mocap trigger
% MocapTrigger = [MocapStart MocapEnd];
% save(strcat('../MatlabData/',SplitName{3}(1:end-6),'-MocapTrigger.mat'),'MocapTrigger')

%% signal map 
% figure()
% for i = 1:length(Par.AnatGrps)
%     for ii = 1:length(Par.AnatGrps(i).Channels)
%         
%         subplot(16,16,i+(ii-1)*16)
%         plot(ACLfp(Sinit:Sfin,Par.AnatGrps(i).Channels(ii)+1))
%         ylim([-7000 7000])
% 
%     end
% end

%% spectrogram
% [s,f,t] = spectrogram(ACLfp(:,196),2^window, nperseg , nfft, LFPfs);
% 
% ax3 = subplot(6,1,3);
% h= pcolor(t,f(3:98),flip(log10(abs(s(3:98,:))),2));
% 
% caxis([3 6])
% % set(gca,'YScale','log')
% set(gca,'YDir','normal')
% colormap('jet')
% set(h, 'EdgeColor', 'none')
% 
% ax4 = subplot(6,1,4);
% plot((1:1:length(ACLfp(:,196)))/LFPfs, ACLfp(:,1))
% linkaxes([ax1,ax2],'x')