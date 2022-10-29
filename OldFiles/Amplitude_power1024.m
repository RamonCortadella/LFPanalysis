% Here the modulation of LFP power by ISA phase is computed for a selected
% brain state
close('all')
clearvars
%% Define settings

LFPfs = 250;
SpeedDCFact = 1;
nCh =1024;

MinDurationState = 1; %Min duration HVS 1s 

Fs=LFPfs;
DownSamp2Pow = 1;
OrderFiltAC = 2;
OrderFiltDC = 2;
GaussianFact = 1;

timeResSpect = 0.05;
SampsResSpect = timeResSpect * Fs;
window = floor(log2(SampsResSpect));
nperov = 2^(window-4);
nfft = 2^window;
fMaxPow = 40;
fMinPow = 20;

FminAC = 1;
FmaxAC = 8;
FminDC = 0.2;
FmaxDC = 8;
SlowDown = 30;

Par =  LoadXml('../../../data/ASIC1024/B14062W18-T1-rat01601/DatData/Rec_map_SIG_B14062W18-T1_54_S-ACinterAC.xml'); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat("../../../data/ASIC1024/B14062W18-T1-rat01601/DatData/",'Videos/delta-beta1024.avi');

%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};


d = dir(strcat("../../../data/ASIC1024/B14062W18-T1-rat01601/DatData/",'Rec_map_SIG_B14062W18-T1_54_S-AC','interAC.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');


ACLfp = [];
DCLfp = [];
for fn = 1:length(d)
    FileName = [d(fn).folder ,'/',d(fn).name];
    Lfp = LoadBinaryDAT(FileName, [0:nCh-1], Par.nChannels,1)';

    SplitName = split(d(fn).name,'-');

    if SplitName{3}(end-5:end-4)=='AC'
        LfpGeom = Lfp;
    end
end
 
nCh = 1024;
nChCol = 32;

ECoG2Video(LfpGeom(floor(length(LfpGeom(:,1))*95.6/100):floor(length(LfpGeom(:,1))*98/100),:), LfpGeom(floor(length(LfpGeom(:,1))*21.6/100):floor(length(LfpGeom(:,1))*22/100),:), videoOutPath, Fs, DownSamp2Pow, 10,10, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,false,GaussianFact,nCh,nChCol)

