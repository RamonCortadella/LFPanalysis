% Here the modulation of LFP power by ISA phase is computed for a selected
% brain state
close('all')
clearvars
%% Define settings

LFPfs = 651.04166667;
DownSampDC = 20;
SpeedDCFact = 10;
SpeedACFact = 10;

MinDurationState = 1; %Min duration HVS 1s 

Fs=LFPfs;
FsDC = Fs/DownSampDC;
ty = 'DC';

DownSamp2Pow = 1;
OrderFiltAC = 2;
OrderFiltDC = 2;
GaussianFact = 0.5;

timeResSpect = 0.05;
SampsResSpect = timeResSpect * Fs;
window = floor(log2(SampsResSpect));
nperov = 2^(window-4);
nfft = 2^window;
fMaxPow = 16;
fMinPow = 8;

ChIndHVS= 229; % centre of the 9 channels averaged

ChNum = 256;

M = [-16 0 16] + [-1 0 1]';
ChIndHVSList= ChIndHVS  + M;


% ChIndTheta = 90;
FminAC = 0.02;
FmaxAC = 0.08;
FminDC = 0.02;
FmaxDC = 0.08;
SlowDown = 30;

nCh = 256;
nChCol = 16;
nChRow = 16;

directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/ClippedMapped/B13289O14-DH1-Rec9interAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat(directory,'MatlabData/Videos/ISAPh-SpindlePow.avi');

%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};
LfpGeomEventsAC = {};

for i = [9]
    d = dir(strcat(directory,'DatData/ClippedMapped/B13289O14-DH1-Rec',int2str(i),'inter*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');


    ACLfp = [];
    DCLfp = [];
    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:255], Par.nChannels,1)';

        SplitName = split(d(fn).name,'-');

        if SplitName{3}(end-5:end-4)== ty
            LfpGeom = Lfp;
        end
    end
    
 
    Periods = load(strcat(directory,'MatlabData/Rec',int2str(i),'-PerStates.mat'));
    
    
    for per = Periods.PerStates.PerSWS' %THE' %SWS
        if per(2)-per(1) >= MinDurationState
            counter = counter +1;
            orderVec(counter) = per(2)-per(1);
                
            if per(1) == 0 
                continue
            end
            if ty == 'DC'
                Fs =FsDC;
            end
            perIndAC = floor(per*Fs);
            LfpGeomEventsAC{counter} = LfpGeom(perIndAC(1):perIndAC(2),:);
            TimeAC{counter} = linspace(perIndAC(1)/Fs,perIndAC(2)/Fs,perIndAC(2)-perIndAC(1)+1);
        end
    end
    
    clearvars LfpGeomDC
    clearvars LfpGeom
end
[x, sortInd] = sort(orderVec);

%% plot average HVS trace and power video
sig = LfpGeomEventsAC{sortInd(counter)};

LenLongestEvent = length(sig(:,1));
t=linspace(0,LenLongestEvent/(Fs),LenLongestEvent+1);
t = t(1:end-1);
%compute power time series

figure
for i =1:ChNum
    display('Compute rms power  ->');display(i);
    [s,f,tspec] = mtcsglong(sig(:,i),nfft,Fs);%,nperov,nfft);
    display(max(max(s)))
    if i == 1
        SmapInt = zeros(length(t),ChNum);
    end
%     EventDC = zeros([LenLongestHVS,ChNum]);
    % resample DC signal to have the same sampling rate  as power
    if length(tspec) >= length(sig(:,1))
        EventDC2 = Interpolate([t', sig(:,i)], tspec,'trim','off');
        indexNaN = find(isnan(EventDC2(:,2)));
        szDC = size(EventDC2);
%         EventDC(1:szDC(1),:) = EventDC2(1:indexNaN-1,2);
        EventDC = EventDC2(1:indexNaN-1,2);
        sig(:,i) = EventDC; 
        
    else
        ind = find((f<=fMaxPow) & (f>= fMinPow));
        s = sum(s(:,ind),2);
        s2 = Interpolate([tspec, s], t,'trim','off');
        indexNaN = find(isnan(s2(:,2)));
        szDC = LenLongestEvent;
%         EventDC(1:szDC(1),:) = EventDC2(1:indexNaN-1,2);
        s = s2(1:indexNaN-1,2);
    end
      
    SmapInt(1:length(s),i) = s;
end
% h= pcolor(t,f(1:NFmaxTheta),log10(abs(s(:,1:NFmaxTheta)')));


ECoG2Video(sig, SmapInt, videoOutPath, Fs, DownSamp2Pow, SpeedACFact,SpeedDCFact, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,false,GaussianFact,nCh,nChCol,nChRow)
% (floor(length(sig(:,1))/2):end,:)
