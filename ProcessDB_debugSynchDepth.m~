close('all')
%% Define Queries
Queries={};
% Queries.depth=1;
Queries.coverage=1;
Queries.flag=1;
Queries.CoupledDC=1;
Queries.mocap=1;
Queries.SingleShank=1;
Queries.LostMocapSamples=0;
Queries.depth=1;
Queries.drifts=0;
Queries.recordingSystem = 0;
Queries.Degraded = 0;
Queries.DeviceName = 'B13289O14O23-DH3SL5';%{'B13289O14-DH1','B13289O14-DH2','B13289O14-DH3','B13289O14-DH5','B13289O24-DH1SL7'};%{'B13289O14-DH1','B13289O14-DH2','B13289O14-DH3','B13289O14-DH5','B13289O24-DH1SL7'};%{-1,'B13289O14-DH1'}; % if cell array with -1 in first field, the opposite of the query is used


extension = '.interp.lfp';
%% Search in DB
FileNameDB = 'RecordingsDB';
[recs, val, ind, T] = SearchDB(['/storage2/ramon/data/DB/' FileNameDB, '.xlsx'], Queries);

%% Apply processing function on retrieved filenames
zsRec = size(recs);

for iFile = 1:zsRec(1)

    display(iFile,'Processing File')

    if T.CoupledAC(ind(iFile)) == 1 & T.CoupledDC(ind(iFile)) ==0
        Coupling ='.AC';
    elseif T.CoupledAC(ind(iFile)) == 0 & T.CoupledDC(ind(iFile)) ==1
        Coupling ='.DC';
    else
        Coupling ='';
    end
        
        
    InputPathMotor = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/Mocap/'];
    
    OutputPathTable = '../../../data/DB/';
    
    InputPath = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/' T.RecordingId{ind(iFile)}];
    OutPathStates = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/' T.RecordingId{ind(iFile)} '/States/'];
    
    FileName = [InputPath,'/',T.FileName{ind(iFile)}(1:end-4),extension];
    MocapFile = [InputPathMotor, T.RecordingId{ind(iFile)},'.csv'];
    
end

% %% compute triggered average
BrainState = 'PerSWS';
InputPath = strcat('../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ECoG','.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
OutputPath = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';

BrainState='PerSWS';
FreqBandLow={'delta'};
FreqBandsHigh={'beta','lowgamma','highgamma'};
HighFreqCoupling='.AC';
FMinSpec =400;
DownFact=[];
Twindow=1;
TimeScaleGlue=0.2;
indDB = ind(iFile);
nCh = T.NumCh(indDB);
nRowsD = T.nRows(indDB);
nColsD = T.nCols(indDB);
nRowsS = T.nRows(indDB-1);
nColsS = T.nCols(indDB-1);
Fs = T.Fs(indDB);
GaussianFact = [];
MedFact = true; % or []

load(InputPath,'out')
fB.isa = [0.02 0.5];
fB.delta = [0.5 5];
fB.theta = [5 9];
fB.spindle = [9 20];
fB.beta = [20 30];
fB.lowgamma = [30 50];
fB.highgamma = [50 200];


% The coupling of the data in the DB search, that must match that in
% InputPath is used to compute the Coupling. If it is AC coupled the Fs is
% corrected as in the ICA computation

if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
    Coupling ='.AC';
elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
    Coupling ='.DC';
else
    Coupling ='';
end

if strcmp(Coupling, '.AC') | strcmp(Coupling, '')
    if T.depth(indDB) == 1 & T.SingleShank(indDB) == 0
        DownFact = 6;
    elseif T.depth(indDB)==0 
        DownFact = 3;
    end
    Fs = Fs/DownFact;

end

%indicate a single frequency band that defines the IC used for the
%triggered average
A = out.(FreqBandLow{1}).A;
display(length(A(1,:)),'Index of IC in the displayed computation: CHECK FOR INCONSISTENCIES!!!')

%display components map
[bX, bY]= meshgrid([1:nColsS]',[1:nRowsS]');
v = reshape(bX,[],1);
w = reshape(bY,[],1);
figure()
title(FreqBandLow{1})
NumICs = length(A(1,:));
for i = 1:NumICs
    subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
    F=  scatteredInterpolant(v,w,A(:,i));
    bF = F(bX,bY);
    imagesc(bF)
    colorbar

    title(['index=' num2str(i)])
end

display('ready*******')



InputPathState = '/storage2/ramon/data/DB/Files/';
Lfp=[];
val = out.val;
if strcmp(Coupling,'.DC') & strcmp(HighFreqCoupling,'.AC')
   FsHigh = T.Fs(indDB-2);
else
   FsHigh = T.Fs(indDB);
end

if DownFact
    
    FsHigh = FsHigh/DownFact;
end

SelectedPeriods = [];
for i = 1:length(val)

    Tend = 0;
    if contains(extension,'interp') %force use of interpolated files
        dn = split(val{1},'-');
        Device = [dn{1} '-' dn{2}];
        if contains(T.FileName(indDB),'ECoG')
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.ECoG' HighFreqCoupling extension];
        elseif  contains(T.FileName(indDB),'Depth')
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.Depth' HighFreqCoupling extension];
        else
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} HighFreqCoupling extension];
        end

        display(FileName, 'loading file for Lfp concatenation across recordings')
        LfpT = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
        
        if DownFact
            LfpT = LfpT(1:DownFact:end,:);
        end
        szLfp = size(Lfp);
        Tend = szLfp(1)/FsHigh;
        display(FsHigh,'*****************')
        Lfp = [Lfp; LfpT];

        if BrainState
            display(indDB)
            display(T.RecordingId{indDB})

            display(strcat(InputPathState,Device,'/',val{i},'/States/',val{i},'.PerStates.mat'))

            PerStates = load(strcat(InputPathState,Device,'/',val{i},'/States/',val{i} ,'.PerStates.mat'));
            SelectedPeriodsTemp = PerStates.PerStates.(BrainState);
            DurationStates = SelectedPeriodsTemp(:,2)-SelectedPeriodsTemp(:,1);
            SelectedPeriodsTemp = SelectedPeriodsTemp(find(DurationStates >=out.MinDurationStates),:); 
            SelectedPeriodsTemp = SelectedPeriodsTemp + Tend;

            SelectedPeriods = [SelectedPeriods; SelectedPeriodsTemp];
%             DurationStates = DurationStates(find(DurationStates >=MinDurationStates));
        end
    else
        display('ATTENTION! : Extension must contain "interp" processing step for comparable recordings!!!!! This process will crash in 100s')
        pause(100)
        break
    end
end

         
LfpState = [];
%concatenate Lfp for selected periods using a zeroing
for i = 1:length(SelectedPeriods(:,1))
    display(size(Lfp),'size Lfp')
    display(floor(FsHigh*SelectedPeriods(i,:)), 'min and max samp Period')
    LfpStateTemp = Lfp(floor(FsHigh*SelectedPeriods(i,1)+1):floor(FsHigh*SelectedPeriods(i,2))-1,:); 

    WindowTemp = ones([length(LfpStateTemp(:,1)),1]);
    Window1 = sigmoid([1:length(LfpStateTemp(:,1))]/FsHigh-10*TimeScaleGlue,TimeScaleGlue);
    Window= WindowTemp.*Window1'.*flip(Window1');
    Window = repmat(Window,1,nCh);
    display(size(LfpStateTemp))
    display(size(Window))
    LfpStateTemp = LfpStateTemp.*Window;


    LfpState = [LfpState; LfpStateTemp];
end
if strcmp(HighFreqCoupling,'.AC')
    display([49 51]/(FsHigh)*2)
    display(FsHigh)
    LfpState = ButFilter(LfpState,2,[49 51]/(FsHigh)*2,'stop');
    LfpState = ButFilter(LfpState,2,[99 101]/(FsHigh)*2,'stop');
    LfpState2 = ButFilter(LfpState,2,[20]/(FsHigh)*2,'low');
    
%     LfpState = ButFilter(LfpState,2,[149 151]/(FsHigh)*2,'stop');
end 



SpatialFilt = false;
if GaussianFact 
    SpatialFilt = true;
elseif MedFact
    SpatialFilt = true;
end
if SpatialFilt
    [bX, bY]= meshgrid([1:3]',[1:nRowsD]');
    LfpState2 = [LfpState,LfpState,LfpState];

    v = reshape(bX,[],1);
    w = reshape(bY,[],1); 
    LfpStateSmooth = [];
    lLfp =length(LfpState2(:,1));
    LfpInterp = zeros([length(Lfp(:,1)), nRowsD*1]);
    [bX, bY]= meshgrid([1:1]',[1:nRowsD]');
    for k =1:lLfp
        display(k)
%         display(lLfp)
        map = LfpState2(k,:);
        F=  scatteredInterpolant(v,w,map');
        bF = F(v,w);
        if GaussianFact
            bF = imgaussfilt(bF,GaussianFact);
        elseif MedFact
            bF = medfilt2(bF);
        end
        LfpIntT = reshape(bF,[],1);
        LfpInterp(k,:) = LfpIntT(1:32);
        display(LfpInterp(k,1:3))
    end
end
LfpState =LfpInterp;

nFFT = 64;
if T.depth(indDB)==1 & T.SingleShank(indDB)==0
    nFFT = 128;
end
if DownFact
    nFFT = nFFT/2;
end
sStruct = struct();
pow=struct();
counter = 0;
Channels = input('ch indices for spectrogram//(all) for video');
if strcmp(Channels,'all')
    Channels = [1:nCh];
end

for iCh = Channels
    counter = counter+1;
    display(iCh,'Computing spectrogram for Channel')

%     display(size(LfpState))
%     display(LfpState(1:4,1:4))
    [s,f,t] = mtcsglong(LfpState(:,iCh),nFFT,FsHigh,nFFT);%,1024);

    if counter == 1
        display(t(2)-t(1),'Time resolution of spectrogram')
        TsSpec = t(2)-t(1);
        FsSpec = 1/(t(2)-t(1));

    end
    sStruct.(['Channel' num2str(iCh)]) = s;
    AvSpec.f = f;
    sStruct.t = t;
    AvSpec.t = t;
    
    for iFreq = 1:length(FreqBandsHigh)%compute power for freq. bands of interest
        indFreq = find(f>=fB.(FreqBandsHigh{iFreq})(1) & f<=fB.(FreqBandsHigh{iFreq})(2));
%         display(indFreq)
%         display(sum(s(:,indFreq),2))
%         
        pow.(FreqBandsHigh{iFreq})(:,iCh) = sum(s(:,indFreq),2);
    end
end


SelectedICs = input('indices of selected ICs to use for triggered average');
SignCorr = input('sign of peaks to detect ([1]/[-1]');
AvPower = struct();
AvLfp = struct();
PowEvent=struct();
LfpEvent=struct();
PowerDistSurf = struct();
LfpDistSurf = struct();

PowerDistDepth = struct();
LfpDist.(FreqBandLow{1}) = struct();
for iSelectICs = SelectedICs %= input('index of IC to further process');
    %input('multiply activation by ([1]/[-1]) // inverse sign if loadings are not correct (i.e. negative for delta waves)');
    % find selected IC peaks

    ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:)*SignCorr;
    ICt = [1:length(ICsig)]/Fs;

    Threshold = 0.1;%input('input the threshold for detection of peaks (it must be positive)');

    [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',Threshold);
    if length(pks)<=50
        continue
        display(['No pks above threshold for IC' num2str(iSelectICs)])
    end
    p = prctile(pks,97);
    figure()
    histogram(pks,100)
    vline(p)
    [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',p);
    pks = pks(1:end-2);
    locs = locs(1:end-2);
    
    [sPks, sLocs] = sort(pks);
    sLocs2 = locs(sLocs);
    % use val of loaded ICs to compute spectrograms
    
    % compute triggered average of the power bands
    SampWidthIC = floor(Twindow*Fs);
    AvIC = zeros([length(-SampWidthIC:SampWidthIC),1]);
    for iE = 3:length(pks)
        if sLocs2(iE) == (locs(1)|locs(2)|locs(3))  
            continue
        end
       
    %if all channels selected compute map metrics for videos
    % define samples for power spectrum corresponding IC activation peak
        SampsTriggerPower = floor(ICt(sLocs2(iE))*FsSpec);
        display(iE)
        display(sPks(iE))
        display(pks(sLocs(iE)))
        SampWidth = floor(Twindow*FsSpec);
        SampsWindow= [SampsTriggerPower-SampWidth:SampsTriggerPower+SampWidth];
    % define samples for Lfp  corresponding IC activation peak
        SampsTriggerLfp = floor(ICt(sLocs2(iE))*FsHigh);
        SampWidthLfp = floor(Twindow*FsHigh);
        SampsWindowLfp= [SampsTriggerLfp-SampWidthLfp:SampsTriggerLfp+SampWidthLfp];
    
        if (sLocs2(iE)-SampWidthIC)<1
            continue
        end
        AvIC = AvIC + ICsig(sLocs2(iE)-SampWidthIC:sLocs2(iE)+SampWidthIC)';

        if iE ==3
            AvLfp.(FreqBandLow{1}) = LfpState2(SampsWindowLfp,:);
        elseif iE == length(pks)
            LfpEvent.(FreqBandLow{1}) = LfpState2(SampsWindowLfp,:);
            ICEvent = ICsig(sLocs2(iE)-SampWidthIC:sLocs2(iE)+SampWidthIC);
        else
            AvLfp.(FreqBandLow{1}) = AvLfp.(FreqBandLow{1}) + LfpState2(SampsWindowLfp,:);
        end

        for iFreq= 1:length(FreqBandsHigh)
%             zeroedPow = repmat(mean(pow.(FreqBandsHigh{iFreq})(SampsWindow,:)),[length(SampsWindow),1]); 
%             NormalizedPow = (pow.(FreqBandsHigh{iFreq})(SampsWindow,:)-zeroedPow)./repmat(max(pow.(FreqBandsHigh{iFreq})(SampsWindow,:)-zeroedPow),[length(SampsWindow),1]);%normalie power bands
            NormalizedPow = zscore(log10(pow.(FreqBandsHigh{iFreq})(SampsWindow,:)));
            if iE ==3
                AvPower.(FreqBandsHigh{iFreq}) = NormalizedPow;
                PowerDistSurf.(FreqBandsHigh{iFreq})= mean(NormalizedPow(:,1:3)')';
                PowerDistDepth.(FreqBandsHigh{iFreq})= mean(NormalizedPow(:,24:27)')';
                if iFreq==1
                    LfpDistSurf = LfpState2(SampsWindowLfp,1)-LfpState2(1,1);
                end
            elseif iE == length(pks)
                PowEvent.(FreqBandsHigh{iFreq}) = NormalizedPow;
            else
                AvPower.(FreqBandsHigh{iFreq}) = AvPower.(FreqBandsHigh{iFreq}) + NormalizedPow;%-repmat(pow.(FreqBandsHigh{iFreq})(SampsWindow(1),:),[length(pow.(FreqBandsHigh{iFreq})(SampsWindow(1),:)),1]);
                PowerDistSurf.(FreqBandsHigh{iFreq}) = [PowerDistSurf.(FreqBandsHigh{iFreq}),mean(NormalizedPow(:,1:3)')'];
                PowerDistDepth.(FreqBandsHigh{iFreq}) = [PowerDistDepth.(FreqBandsHigh{iFreq}),mean(NormalizedPow(:,24:27)')'];
                if iFreq == 1
                    LfpDistSurf = [LfpDistSurf, LfpState2(SampsWindowLfp,1)-LfpState2(1,1)];
                end
            end
        end
        
        for iCh = Channels
%             zeroedEvent = repmat(mean(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)),[length(sStruct.(['Channel' num2str(iCh)])(SampsWindow,1)),1]);
%             whiteningEvent = repmat(max(abs(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)-zeroedEvent)),[length(sStruct.(['Channel' num2str(iCh)])(SampsWindow,1)),1]);
%             whitenedEvent = (sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)-zeroedEvent)./whiteningEvent;
            whitenedEvent = zscore(log10(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)));
            if iE ==3
                AvSpec.(['Channel' num2str(iCh)]) = whitenedEvent;
                LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]) = LfpState2(SampsWindowLfp,iCh)-LfpState2(SampsTriggerLfp,iCh);
            else
                AvSpec.(['Channel' num2str(iCh)]) =  AvSpec.(['Channel' num2str(iCh)]) + whitenedEvent;
                LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]) = [LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]),LfpState2(SampsWindowLfp,iCh)-LfpState2(SampsTriggerLfp,iCh)];
            end
        end
        
    end

    
    AvIC = AvIC/length(pks);
    %display the time evolution of AvPOwer for the different freq bands usign a
    %video
    out.AvIC = AvIC;
    out.ICt = [-SampWidthIC:SampWidthIC]/Fs;
    out.AvPower = AvPower;
    out.PowerT = [-SampWidth:SampWidth]/FsSpec;
    LfpT = [-SampWidthLfp:SampWidthLfp]/FsHigh;
%     save([OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-AvPower' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ').mat'],'out')
    if length(Channels)==32
%         for iFreq= 1:length(FreqBandsHigh)
            figure(200)
            maxAvPow=0;
            maxAvLfp =0;
            hold on
%             subplot(33,1,1)
%             plot(out.ICt,AvIC,'k')
%             legend(['average IC' num2str(iSelectICs) 'activation'])
            for iCh = Channels
%                 subplot(33,1,iCh+1)
                yyaxis left
                maxAvPow = maxAvPow+0.6*max(max(abs(AvPower.(FreqBandsHigh{2}))));
                plot(out.PowerT, AvPower.(FreqBandsHigh{3})(:,iCh)-maxAvPow,'k-')
                ylim([-1.05*maxAvPow 0])
                xlim([out.PowerT(1) out.PowerT(end)])
                yyaxis right
                maxAvLfp = maxAvLfp+0.6*max(max(abs(AvLfp.(FreqBandLow{1}))));
                plot(LfpT, AvLfp.(FreqBandLow{1})(:,iCh)-maxAvLfp,'r-')
                ylim([-1.05*maxAvLfp 0])
                xlim([out.PowerT(1) out.PowerT(end)])
            end
%             ch = [1+1:length(Channels)-1];
%             csd1 = (-1/1)*[ AvLfp.(FreqBandLow{1})(:,ch-1) - 2*AvLfp.(FreqBandLow{1})(:,ch) + AvLfp.(FreqBandLow{1})(:,ch+1) ];
%             ch = [2+1:length(Channels)-2];
%             csd2 = (-1/4)*[ AvLfp.(FreqBandLow{1})(:,ch-2) - 2*AvLfp.(FreqBandLow{1})(:,ch) + AvLfp.(FreqBandLow{1})(:,ch+2) ];
            ch = [2+1:length(Channels)-2];
            csd3 = (-1/7)*[ 2*AvLfp.(FreqBandLow{1})(:,ch-2) -  AvLfp.(FreqBandLow{1})(:,ch-1) - 2*AvLfp.(FreqBandLow{1})(:,ch) - AvLfp.(FreqBandLow{1})(:,ch+1) + 2* AvLfp.(FreqBandLow{1})(:,ch+2) ];
            figure()
            imagesc(csd3')
            colormap('jet')
            colorbar()
            clim([-0.8*max(max(csd3)) 0.8*max(max(csd3))])

%         end
%         for iFreq= 1:length(FreqBandsHigh)
%             figure(201)
%             hold on
% %             subplot(33,1,1)
% %             plot(out.ICt,AvIC,'k')
% %             legend(['average IC' num2str(iSelectICs) 'activation'])
%             for iCh = Channels
% %                 subplot(33,1,iCh+1)
%                 
%             end
%         end
%         for iFreq= 1:length(FreqBandsHigh)
        figure()
        imagesc(PowerDistSurf.(FreqBandsHigh{3}))
        clim([0.7*min(min(PowerDistSurf.(FreqBandsHigh{3}))) 0.7*max(max(PowerDistSurf.(FreqBandsHigh{3})))])
        hline(floor(length(SampsWindow)/2))
        hline(floor(length(SampsWindow)/2)+7)
        figure()
        imagesc(PowerDistDepth.(FreqBandsHigh{3}))
        clim([0.7*min(min(PowerDistDepth.(FreqBandsHigh{3}))) 0.7*max(max(PowerDistDepth.(FreqBandsHigh{3})))])
        hline(floor(length(SampsWindow)/2))
        hline(floor(length(SampsWindow)/2)+7)
        figure()
        imagesc(LfpDistSurf)
        clim([0.7*min(min(LfpDistSurf)) 0.7*max(max(LfpDistSurf))])
        hline(SampWidthLfp)
        hline(SampWidthLfp+200)
        figure()
        hold on
        scatter(LfpDistSurf(SampWidthLfp,:),PowerDistSurf.(FreqBandsHigh{3})(floor(length(SampsWindow)/2),:))
        
        scatter(LfpDistSurf(SampWidthLfp+200,:),PowerDistSurf.(FreqBandsHigh{3})(floor(length(SampsWindow)/2)+7,:),'r')
%       
%         end
    end
    
    
    for iCh = Channels
        pause(2)
        figure()
        ax1 = subplot(5,1,1);
        plot(out.ICt,AvIC)
        legend(['averaged IC' num2str(iSelectICs)])
        ax2 = subplot(5,1,2);
        plot(LfpT,LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]))
        legend(['LFP across events' num2str(iSelectICs)])
        ax3 = subplot(5,1,3);
        plot(LfpT,mean(LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)])'))
        legend(['averaged LFP' num2str(iSelectICs)])
        ax4 = subplot(5,1,4);
        indFreqSpec = find(AvSpec.f<=FMinSpec);
        whitening = repmat(max(abs(AvSpec.(['Channel' num2str(iCh)]))),[length(AvSpec.(['Channel' num2str(iCh)])(:,1)),1]);
        whitened = AvSpec.(['Channel' num2str(iCh)])'./whitening';
        imagesc(out.PowerT,AvSpec.f(indFreqSpec),whitened(indFreqSpec,:))

        title(['Av. Spec. channel' num2str(iCh) '-IC' num2str(iSelectICs)])
        ax5 = subplot(5,1,5);
        hold on
        for iFreq = 1:length(FreqBandsHigh)
            PowTrace=out.AvPower.(FreqBandsHigh{iFreq})(:,iCh);
            plot(out.PowerT,PowTrace/max(abs(PowTrace)),'DisplayName',FreqBandsHigh{iFreq})
        end
        linkaxes([ax1,ax2,ax3,ax4],'x')
    end
end


