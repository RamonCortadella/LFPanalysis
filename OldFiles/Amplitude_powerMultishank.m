% Here the modulation of LFP power by ISA phase is computed for a selected
% brain state
close('all')
clearvars
%% Define settings

LFPfs = 12500;
SpeedDCFact = 10;
SpeedACFact = 10;

nCh =256;

MinDurationState = 1; %Min duration HVS 1s 

Fs=LFPfs;
DownSamp2Pow = 1;
OrderFiltAC = 2;
OrderFiltDC = 2;
GaussianFact = 0.5;

timeResSpect = 0.05;
SampsResSpect = timeResSpect * Fs;
window = floor(log2(SampsResSpect));
nperov = 2^(window-4);
nfft = 2^window;
fMaxPow = 5;
fMinPow = 20;

FminAC = 1;
FmaxAC = 8;
FminDC = 8;
FmaxDC = 20;
SlowDown = 30;

% Par =  LoadXml('../../../data/ASICintraMUX/B13907W21-T1-rat01601/DatData/SIG_040522_B13907W23-M4_InVitro-5-interAC.xml'); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
% videoOutPath = strcat("../../../data/ASICintraMUX/B13907W21-T1-rat01601/DatData/",'Videos/deltaTop-SpindleBot-Gauss0p5-CnstScale-MMF-512.avi');
Par =  LoadXml('../../../data/ASICintraMUX/B13907W21-T1-rat01601/DatData/2022.06.14-22.46.45-Rec-M6-Rec1.dat.xml'); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat("../../../data/ASICintraMUX/B13907W21-T1-rat01601/DatData/",'Videos/deltaTop-SpindleBot-CnstScale-MSh.avi');

%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};


d = dir(strcat("../../../data/ASICintraMUX/B13907W21-T1-rat01601/DatData/",'2022.06.14-22.46.45-Rec-M6-Rec1.dat.dat'));%DC-LP30Hz-Notch50-100Hz.dat');

nCh = 256;
nChCol = 32;
nChRow = 8;

ACLfp = [];
DCLfp = [];
for fn = 1:length(d)
    FileName = [d(fn).folder ,'/',d(fn).name];
    Lfp = LoadBinaryDAT(FileName, [0:nCh-1], Par.nChannels,1)';
%     LfpBP = ButFilter(Lfp,2,[1 50]/Fs*2,'bandpass');
    SplitName = split(d(fn).name,'-');

    LfpGeom = zeros([size(Lfp)]);
    
    for i = 1:length(Par.AnatGrps)
        for ii = 1:length(Par.AnatGrps(i).Channels)
            LfpGeom(:,ii+(i-1)*nChCol) = Lfp(:,Par.AnatGrps(i).Channels(ii)+1);
            display(Par.AnatGrps(i).Skip(ii))
        end
    end
    
end
 

ECoG2Video(LfpGeom, LfpGeom, videoOutPath, Fs, DownSamp2Pow, SpeedDCFact,SpeedACFact, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,false,GaussianFact,nCh,nChCol,nChRow,false)

