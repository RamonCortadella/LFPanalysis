% Here the modulation of LFP power by ISA phase is computed for a selected
% brain state
close('all')
clearvars
%% Define settings

LFPfs = 1000;
SpeedDCFact = 1;
nCh =512;
nChCol = 32;
nChRow = 16;

MinDurationState = 1; %Min duration HVS 1s 

Fs=LFPfs;
DownSamp2Pow = 1;
OrderFiltAC = 2;
OrderFiltDC = 2;
GaussianFact = 0.5;
double = true;

timeResSpect = 0.05;
SampsResSpect = timeResSpect * Fs;
window = floor(log2(SampsResSpect));
nperov = 2^(window-4);
nfft = 2^window;
fMaxPow = 5;
fMinPow = 20;

FminAC = 5;
FmaxAC = 20;
FminDC = 0.01;
FmaxDC = 20;
SlowDown = 30;

Par =  LoadXml('../../../data/ASIC512/B13907W21-T1-rat01601/DatData/SIG_B13907W21-T1_InVivo-9_F0-ACinterAC.xml'); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat("../../../data/ASIC512/B13907W21-T1-rat01601/DatData/",'Videos/bipolar10Hz-512-LimitedScale.avi');

%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};


d = dir(strcat("../../../data/ASIC512/B13907W21-T1-rat01601/DatData/",'SIG_B13907W21-T1_InVivo-9_F0-AC','interAC.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');


ACLfp = [];
DCLfp = [];
for fn = 1:length(d)
    FileName = [d(fn).folder ,'/',d(fn).name];
    Lfp = LoadBinaryDAT(FileName, [0:nCh-1], Par.nChannels,1)';

    SplitName = split(d(fn).name,'-');

    if SplitName{4}(end-5:end-4)=='AC'
        LfpGeom = Lfp;
    end
end
 

ECoG2Video(LfpGeom(floor(length(LfpGeom(:,1))*(60)/100):floor(length(LfpGeom(:,1))*(60.3)/100),:), LfpGeom(floor(length(LfpGeom(:,1))*(60)/100):floor(length(LfpGeom(:,1))*(60.3)/100),:), videoOutPath, Fs, DownSamp2Pow, 1,1, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,true,GaussianFact,nCh,nChCol,nChRow,double)

