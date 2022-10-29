% Here the HVS events are charcterized
close('all')
clearvars
%% Define settings
stats = false;

LFPfs = 651.04166667;
DownSampDC = 100;
SpeedDCFact = 1;

MinDurationState = 1; %Min duration HVS 1s 

Fs=LFPfs;
FsDC = Fs/DownSampDC;
DownSamp2Pow = 1;
OrderFiltAC = 2;
OrderFiltDC = 2;
GaussianFact = 0.5;

timeResSpect = 0.6;
SampsResSpect = timeResSpect * Fs;
window = floor(log2(SampsResSpect));
nperov = 2^(window-4);
nfft = 2^window;
fMaxHVS = 8;
fMinHVS = 5;

ChIndHVSAnt= 229; % centre of the 9 channels averaged
ChIndHVSPost= 235; % centre of the 9 channels averaged

ChNum = 256;

M = [-16 0 16] + [-1 0 1]';
ChIndHVSListAnt= ChIndHVSAnt  + M;
ChIndHVSListPost= ChIndHVSPost  + M;


% ChIndTheta = 90;
FminAC = 0.01;
FmaxAC = 0.5;
FminDC = 0.005;
FmaxDC = 0.5;
SlowDown = 1;

directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec9interAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat(directory,'MatlabData/Videos/HVSamp-power-nogaussian.avi');

%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};

for i = [1 2 3 4 5 6 9]
    d = dir(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec',int2str(i),'inter*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');
    
    display('Rec number');display(i);
    
    ACLfp = [];
    DCLfp = [];
    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:ChNum-1], Par.nChannels,1)';

        SplitName = split(d(fn).name,'-');

        if SplitName{3}(end-5:end-4)=='AC'
            LfpGeom = Lfp;
        elseif  SplitName{3}(end-5:end-4)=='DC' 
            LfpGeomDC = Lfp;
        end
    end
    
 
    Periods = load(strcat(directory,'MatlabData/Rec',int2str(i),'-PerStates.mat'));
    
    
    for per = Periods.PerStates.PerHVS' %THE' %SWS
        if per(2)-per(1) >= MinDurationState
            counter = counter +1;
            orderVec(counter) = per(2)-per(1);
%                 
            if per(1) == 0 
                continue
            end
            perIndAC = floor(per*Fs);
            perIndDC = floor(per*FsDC);
            LfpGeomEventsAC{counter} = LfpGeom(perIndAC(1):perIndAC(2),:);
            TimeAC{counter} = linspace(perIndAC(1)/Fs,perIndAC(2)/Fs,perIndAC(2)-perIndAC(1)+1);
            LfpGeomEventsDC{counter} = LfpGeomDC(perIndDC(1):perIndDC(2),:);
            TimeDC{counter} = linspace(perIndDC(1)/FsDC,perIndDC(2)/FsDC,perIndDC(2)-perIndDC(1)+1);
        end
    end
    
    clearvars LfpGeomDC
    clearvars LfpGeom
end
[x, sortInd] = sort(orderVec);

figure
hist(orderVec,20)
%% plot average HVS trace Anterior/Posterior
LenLongestHVS = length(mean(LfpGeomEventsDC{sortInd(counter)}(:,ChIndHVSListAnt),2));

EventDCAnt = zeros([LenLongestHVS,1]);
EventDCPost = zeros([LenLongestHVS,1]);

sz = length(EventDCAnt);
HVSStackAnt = zeros([counter,sz]);

sz = length(EventDCPost);
HVSStackPost = zeros([counter,sz]);


figure
for i = 1:counter
    display('Stacking Anterior/Posterior, event number');display(i);
    
    EventDC2Ant = mean(LfpGeomEventsDC{sortInd(i)}(:,ChIndHVSListAnt),2)-mean(LfpGeomEventsDC{sortInd(i)}(1,ChIndHVSListAnt),2);
    EventDC2Post = mean(LfpGeomEventsDC{sortInd(i)}(:,ChIndHVSListPost),2)-mean(LfpGeomEventsDC{sortInd(i)}(1,ChIndHVSListPost),2);
    if any(abs(EventDC2Ant)>=1500)
        continue
    end
    if any(abs(EventDC2Post)>=1500)
        continue
    end
    EventDCAnt(1:length(EventDC2Ant)) = EventDC2Ant;
    EventDCPost(1:length(EventDC2Post)) = EventDC2Post;

%     if i == 1
%         HVSStackAnt = EventDCAnt;
%         HVSStackPost = EventDCPost;
%     else
%         HVSStackAnt = cat(2,HVSStackAnt,EventDCAnt);
%         HVSStackPost = cat(2,HVSStackPost,EventDCPost);
%     end
    HVSStackAnt(i,:) = EventDCAnt;
    HVSStackPost(i,:) = EventDCPost;
    
    EventDCAnt = zeros([LenLongestHVS,1]); %reset before next loop
    EventDCPost = zeros([LenLongestHVS,1]);
end
MeanHVSAnt = mean(HVSStackAnt,1);
MeanHVSPost = mean(HVSStackPost,1);

figure
hold on
t=linspace(0,LenLongestHVS/(FsDC),LenLongestHVS+1);
t = t(1:end-1);
plot(t,MeanHVSAnt,'b')
plot(t,MeanHVSPost,'r')
legend('Anterior','Posterior')
xlabel('time(s)')
hold off
figure
hold on
plot(t,HVSStackAnt,'b')
plot(t,MeanHVSAnt,'r')
ylim([-300 300])
hold off
%% plot average HVS trace and power video

LenLongestHVS = length(mean(LfpGeomEventsDC{sortInd(counter)}(:,ChIndHVSListAnt),2));
EventDC = zeros([LenLongestHVS,ChNum]);

LenLongestHVSAC = length(mean(LfpGeomEventsAC{sortInd(counter)}(:,ChIndHVSListAnt),2));
EventAC = zeros([LenLongestHVSAC,ChNum]);

sz = size(EventDC);
HVSStack = zeros([counter,sz(1),sz(2)]);

sz = size(EventAC);
HVSStackAC = zeros([counter,sz(1),sz(2)]);

for i = 1:counter
    display('Compute HVS mean ->');display(i);

    EventDC2 = LfpGeomEventsDC{sortInd(i)}-LfpGeomEventsDC{sortInd(i)}(1,:);
    sz = size(EventDC2);
    EventDC(1:sz(1),:) = EventDC2;

    EventAC2 = LfpGeomEventsAC{sortInd(i)};
    szAC = size(EventAC2);
    EventAC(1:szAC(1),:) = EventAC2;

    if any(abs(EventDC)>=1500)
        continue
    end

%     if i == 1
%         HVSStack = EventDC;
%         HVSStackAC = EventAC;
%     else
%         HVSStack = cat(3,HVSStack,EventDC);
%         HVSStackAC = cat(3,HVSStackAC,EventAC);
%     end
    HVSStack(i,:,:) = EventDC;
    HVSStackAC(i,:,:) = EventAC;
    
    EventDC = zeros([LenLongestHVS,ChNum]); %reset before next loop
    EventAC = zeros([LenLongestHVSAC,ChNum]); %reset before next loop

end
MeanHVS = mean(HVSStack,1);

MeanHVSAC = mean(HVSStackAC,1);

%compute power HVS time series

figure
for i =1:ChNum
    display('Compute rms power HVS ->');display(i);
    [s,f,tspec] = mtcsglong(MeanHVSAC(1,:,i),nfft,Fs);%,nperov,nfft);
    if i == 1
        SmapInt = zeros(length(t),ChNum);
    end
%     EventDC = zeros([LenLongestHVS,ChNum]);
    % resample DC signal to have the same sampling rate  as power
    if length(tspec) >= length(HVSStack(1,:,1))
        EventDC2 = Interpolate([t', MeanHVS(:,i)], tspec,'trim','off');
        indexNaN = find(isnan(EventDC2(:,2)));
        szDC = size(EventDC2);
%         EventDC(1:szDC(1),:) = EventDC2(1:indexNaN-1,2);
        EventDC = EventDC2(1:indexNaN-1,2);
        MeanHVS(:,i) = EventDC; 
        
    else
        ind = find((f<=fMaxHVS) & (f>= fMinHVS));
        s = sum(s(:,ind),2);
        s2 = Interpolate([tspec, s], t,'trim','off');
        indexNaN = find(isnan(s2(:,2)));
        szDC = size(EventDC2);
%         EventDC(1:szDC(1),:) = EventDC2(1:indexNaN-1,2);
        s = s2(1:indexNaN-1,2);
    end
      
    SmapInt(1:length(s),i) = s;
end
% h= pcolor(t,f(1:NFmaxTheta),log10(abs(s(:,1:NFmaxTheta)')));

MeanHVS = squeeze(MeanHVS(1,:,:));
ECoG2Video(MeanHVS, SmapInt, videoOutPath, FsDC, DownSamp2Pow, SpeedDCFact, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,true,GaussianFact)

