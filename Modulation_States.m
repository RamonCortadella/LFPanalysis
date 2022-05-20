% Here the modulation of LFP power by ISA phase is computed for a selected
% brain state
close('all')
clearvars
%% Define settings
stats = false;

LFPfs = 651.04166667;
DownSampDC = 100;

MinDurationState = 100; %Min duration SWS 200s %Min duration REM 55s
MinDurationStateREM = 55;

timeRes = 10;
timeResDC = 100;

Fs=LFPfs;
FsDC = Fs/DownSampDC;
SpeedDCFact=100;
SpeedACFact = 1;
OrderFiltAC = 4;
OrderFiltDC = 2;

ChIndHVS= 229; % centre of the 9 channels averaged
ChIndTheta = 93; % 93 initially used (211 to detect modulation of spindle by ISA)
ChIndHVS = 229;
% ChIndTheta = 90;
FminAC = 8;
FmaxAC = 16;
FminDC = 0.02;
FmaxDC = 0.08;
SlowDown = 30;
CnstScale = false;
GaussianFact= false;

% directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
% Par =  LoadXml(strcat(directory,'DatData/Clipped/B13289O14-DH1-Rec9interAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
% videoOutPath = strcat(directory,'MatlabData/Videos/EventExample.avi');


directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/ClippedMapped/B13289O14-DH1-Rec3interAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat(directory,'MatlabData/Videos/SWS-ISA-spindle.avi');
%% Get ephys and states data and organize events ephys 

counter = 0;
counterREM = 0;
LfpGeomEvents = {};

for i = [3]
    d = dir(strcat(directory,'DatData/ClippedMapped/B13289O14-DH1-Rec',int2str(i),'inter*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');


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
    
    
    for per = Periods.PerStates.PerSWS' %THE' %SWS
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

%%% THETA modulation of GAMMA
fPh = [bsxfun(@plus,[1:0.3:10],[-0.3 0.3]')']/FNi;
fPow = bsxfun(@plus, [15:3:100],[-5 5]')'/FNi;

%%% REM 
% fPow = [bsxfun(@plus,[1:0.5:15],[-0.3 0.3]')']/FNi;
% fPh1 = bsxfun(@plus,[0.03:0.01:0.2], [-0.02 0.02]')'/FNi;
%%% SWS
% fPow = [bsxfun(@plus,[6:0.3:20],[-1.5 1.5]')']/FNi;
% fPh = bsxfun(@plus,[0.03:0.01:0.2], [-0.02 0.02]')'/FNi;
%%% SWS
% fPow = [bsxfun(@plus,[1:0.3:15],[-0.3 0.3]')']/FNi;
% fPh = bsxfun(@plus,[0.03:0.01:0.2], [-0.02 0.02]')'/FNi;
%%% ISA-gamma
% fPow = bsxfun(@plus,[30:2:70],[-6 6]')'/FNi;

% fPh = [bsxfun(@plus,[1.0:0.3:8],[-0.3 0.3]')']/FNi;
% fPow = bsxfun(@plus,[2.5:0.3:15],[-1.5 1.5]')'/FNi;
% fPh = bsxfun(@plus,[3:0.5:17], [-2 2]')'/FNi;

%%%HVS
% fPow = [bsxfun(@plus,[3:0.5:25],[-1 1]')']/FNi;
% fPh = bsxfun(@plus,[0.05:0.01:0.3], [-0.02 0.02]')'/FNi;

M = [-32 -16 0 16 32] + [-2 -1 0 1 2]';
ChIndThetaList= ChIndTheta  + M;

%% independently compute modulation for all events
if stats == false
    figure
    for i = 1:counter
        EventDC2 = Interpolate([TimeDC{sortInd(i)}', mean(LfpGeomEventsDC{sortInd(i)}(:,ChIndThetaList),2)], TimeAC{sortInd(i)},'trim','off');
        indexNaN = find(isnan(EventDC2(:,2)));
        EventDC = EventDC2(1:indexNaN-1,2);

        [out, PhMat, PowMat] = PowerPhasePairsR(EventDC,  fPh, mean(LfpGeomEventsAC{sortInd(i)}(1:indexNaN-1,ChIndThetaList),2) , fPow, 1,'but',@PowerModulation);
        if i == 1
            RampStack = out.Ramp;
        else
            RampStack = cat(3,RampStack,out.Ramp);
        end

        ax = subplot(floor(sqrt(counter))+1,floor(sqrt(counter))+1,i);
        pcolor(out.fPh(1:end)*FNi, out.fPow*FNi, out.Ramp(1:end,:)'); shading flat; 
        colorbar
    end
    figure
    p = pcolor(out.fPh(1:end)*FNi, out.fPow*FNi, mean(RampStack,3)'); p.FaceColor='interp'; shading flat; 
    colorbar

    ECoG2Video(LfpGeomEventsAC{sortInd(i)}, LfpGeomEventsDC{sortInd(i)}, videoOutPath, Fs, 100, SpeedDCFact,SpeedACFact, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,CnstScale, GaussianFact)
else
%% compute the modulation of all concatenated events
    Shuffle.MaxShift= 100*100/10; % after 10 times resampling will have 50sec shift.
    Shuffle.Type = 'shift'; Shuffle.nShuffle = 50;

    figure
    for i = 1:counter
        EventDC2 = Interpolate([TimeDC{sortInd(i)}', mean(LfpGeomEventsDC{sortInd(i)}(:,ChIndThetaList),2)], TimeAC{sortInd(i)},'trim','off');
        indexNaN = find(isnan(EventDC2(:,2)));
        EventDC = EventDC2(1:indexNaN-1,2);
        display(i);
        if i == 1
            EventStackDC = EventDC';
            EventStackAC = mean(LfpGeomEventsAC{sortInd(i)}(1:indexNaN-1,ChIndThetaList),2)';
        else
            EventStackDC = [EventStackDC EventDC'];
            EventStackAC = [EventStackAC mean(LfpGeomEventsAC{sortInd(i)}(1:indexNaN-1,ChIndThetaList),2)'];
        end

    end

    [out4sh,PhMat, PowMat] = PowerPhasePairsR(EventStackDC',  fPh, EventStackAC' , fPow, 1,'but',@PowerModulation,Shuffle);
end
% 
figure
p = pcolor(out4sh.fPh(2:end)*FNi, out4sh.fPow*FNi, out4sh.Ramp(2:end,:)'); p.FaceColor='interp'; shading flat; 
colorbar

figure
p = pcolor(out4sh.fPh(2:end)*FNi, out4sh.fPow*FNi, out4sh.Rpval(2:end,:)'); p.FaceColor='interp'; shading flat; 
colorbar
caxis([0 0.01])

figure
inds = find(out4sh.Rpval ~= 0);
out4sh.Ramp(inds) = NaN;
p=pcolor(out4sh.fPh(2:end)*FNi, out4sh.fPow*FNi, out4sh.Ramp(2:end,:)'); p.FaceColor='interp'; shading flat; 
colorbar
set(gca,'Color',uint8([170 170 170]))

%% plot polar distribution and phase
figure(3)
polarplot(out4sh.phbins(:,1,1),mean(mean(sq(out4sh.pow_dens(:,14:15,27:29)),2),3)'); axis tight; box off

figure

% for i = 13:19
    
polarplot(out4sh.phbins(:,1,1),sq(out4sh.pow_dens(:,14,27:29))'); axis tight; box off

[MaxMod, phInd] = max(out4sh.pow_dens);

for iph =1:length(out4sh.pow_dens(1,:,1))
    for ipow = 1:length(out4sh.pow_dens(1,1,:))
        PrefPh(iph,ipow) = out4sh.phbins(phInd(1,iph,ipow));
    end
end

figure
p=pcolor(out4sh.fPh*FNi, out4sh.fPow*FNi, PrefPh'); p.FaceColor='interp'; shading flat; 
colormap jet
caxis([-3 3])
colorbar

% end
% ph = angle(hilbert(ButFilter(EventDC, 2, fPh(3,:), 'bandpass')));
% hold on
% plot(EventDC2(1:indexNaN-1,1),ph)