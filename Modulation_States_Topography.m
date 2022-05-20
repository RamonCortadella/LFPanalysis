% Here the modulation of LFP power by ISA phase is computed for a selected
% brain state
close('all')
clearvars
%% Define settings
stats = true;

LFPfs = 651.04166667;
DownSampDC = 100;
NumCh = 256;

MinDurationState = 1; %Min duration SWS 200s %Min duration REM 55s1
MinDurationStateREM = 55;

timeRes = 10;
timeResDC = 100;

Fs=LFPfs;
FsDC = Fs/DownSampDC;
SpeedDCFact=10;
OrderFiltAC = 4;
OrderFiltDC = 2;

ChIndHVS= 229; % centre of the 9 channels averaged
ChIndTheta = 229;%93 for theta
% ChIndTheta = 90;
FminAC = 2;
FmaxAC = 10;
FminDC = 0.05;
FmaxDC = 0.2;
SlowDown = 30;
CnstScale = false;
% GaussianFact= 0.5; % define below

directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec9interAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat(directory,'MatlabData/Videos/EventExample.avi');

%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};

for i = [1 2 3 4 5 6 9]
    d = dir(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec',int2str(i),'inter*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');


    ACLfp = [];
    DCLfp = [];
    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:255], Par.nChannels,1)';

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

%% Modulation settings
close('all')
FNi = Fs/2;
FNiDC = FsDC/2;

%%% delta-spindle/beta
% fPow = bsxfun(@plus,[23 24],[-10 10]')'/FNi;
% fPh = bsxfun(@plus,[2 3], [-1 1]')'/FNi;

%%% ISA-spindle
% fPow = bsxfun(@plus,[11 12],[-2 2]')'/FNi;
% fPh = bsxfun(@plus,[0.04 0.05], [-0.03 0.03]')'/FNi;

%%% slow-spindle
% fPow = bsxfun(@plus,[11 12],[-2 2]')'/FNi;
% fPh = bsxfun(@plus,[0.1 0.12], [-0.04 0.04]')'/FNi;

%%% ISA-delta

% fPow = bsxfun(@plus,[3 4],[-2 2]')'/FNi;
% fPh = bsxfun(@plus,[0.04 0.05], [-0.03 0.03]')'/FNi;

%%% HVS
fPow = bsxfun(@plus,[16.75 17.5],[-0.75 0.75]')'/FNi;
fPh = bsxfun(@plus,[0.19 0.2], [-0.01 0.01]')'/FNi;

M = [-16 0 16] + [-1 0 1]';
ChIndThetaList= ChIndTheta  + M;

Shuffle.MaxShift= 100*100/10; % after 10 times resampling will have 50sec shift.
Shuffle.Type = 'shift'; Shuffle.nShuffle = 50;

GaussianFact = 0.5;
Coord(:,1) = reshape(repmat([1:16],16,1),[],1);
Coord(:,2) = repmat([1:16]',16,1);
[bX, bY]= meshgrid([1:16]',[1:16]');

%% compute the modulation of all concatenated events
for i = 1:counter
    EventDC2 = Interpolate([TimeDC{sortInd(i)}', LfpGeomEventsDC{sortInd(i)}(:,1)], TimeAC{sortInd(i)},'trim','off');
    indNaNs = find(isnan(EventDC2(:,2)));
    if length(indNaNs)>=1
        indexNaN(i) = indNaNs(1);
    else
        indexNaN(i) = length(EventDC2(:,2));
    end
end
indexNaNSum = sum(indexNaN);
EventStackDC = zeros(indexNaNSum-1,NumCh);
EventStackAC = zeros(indexNaNSum-1,NumCh);

figure
lengthcounter = 1;
for i = 1:counter
    for ii = 1:NumCh
        EventDC2 = Interpolate([TimeDC{sortInd(i)}', LfpGeomEventsDC{sortInd(i)}(:,ii)], TimeAC{sortInd(i)},'trim','off');
        indNaNs = find(isnan(EventDC2(:,2)));
        if length(indNaNs)>=1
            indexNaN = indNaNs(1);
        else
            indexNaN = length(EventDC2(:,2));
        end
        EventDC = EventDC2(1:indexNaN-1,2);
        display(i);display(ii);
        
        EventStackDC(lengthcounter:lengthcounter+indexNaN-2,ii) = EventDC';
        EventStackAC(lengthcounter:lengthcounter+indexNaN-2,ii)= LfpGeomEventsAC{sortInd(i)}(1:indexNaN-1,ii)';
        
    end
    lengthcounter = lengthcounter + indexNaN;
end
% sz = size(EventStackDC);
% for k =1:sz(1)
%     map = EventStackDC(k,:);   
%     F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
%     bF = F(bX,bY);
%     bF = imgaussfilt(bF,GaussianFact);
% end
%% 
for i =1:NumCh
    display(i);
    [out4sh{i},PhMat{i}, PowMat{i}] = PowerPhasePairsR(EventStackDC(:,i),  fPh, EventStackAC(:,i) , fPow, 1,'but',@PowerModulation,Shuffle);
end

%% Plot modulation strength geometrical map
for i =1:NumCh
    ModStr(i) = out4sh{i}.Ramp(1,1);
    ModStrPval(i) = out4sh{i}.Ramp(1,1);
    ModPval(i) = out4sh{i}.Rpval(1,1);
end


inds = find(ModPval ~= 0);
ModStrPval(inds) = NaN;

figure(111)
map = ModStr;   
F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
bF = F(bX,bY);
imagesc(bF);


figure(112)
map = ModStrPval;   
F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
bF = F(bX,bY);
% p = imagesc(bF);
% colorbar
% set(gca,'color',uint8([170 170 170]))

p=pcolor([1:16], [1:16], bF); p.FaceColor='interp'; shading flat; 
colorbar
% caxis([max(max(bF))/7 max(max(bF))])
set(gca,'Color',uint8([170 170 170]))
set(gca,'YDir','reverse')