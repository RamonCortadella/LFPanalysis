close('all')
%% Define Queries
Queries={};
% Queries.depth=1;
Queries.coverage=1;
Queries.flag=1;
Queries.CoupledDC=1;
Queries.mocap=1;
Queries.SingleShank=0;
Queries.LostMocapSamples=0;
Queries.depth=0;
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
% InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'','.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';

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
nRows = T.nRows(indDB);
nCols = T.nCols(indDB);
Fs = T.Fs(indDB);
GaussianFact = [];
MedFact = true; % or []

load(InputPath,'out')
fB.isa = [0.02 0.5];
fB.delta = [0.5 5];
fB.theta = [5 9];
fB.spindle = [9 20];
fB.beta = [20 30];
fB.lowgamma = [30 60];
fB.highgamma = [60 200];


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
[bX, bY]= meshgrid([1:nCols]',[1:nRows]');
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
   FsHigh = T.Fs(indDB-1);
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
%     LfpState = ButFilter(LfpState,2,[149 151]/(FsHigh)*2,'stop');
end

SpatialFilt = false;
if GaussianFact 
    SpatialFilt = true;
elseif MedFact
    SpatialFilt = true;
end
if SpatialFilt
    [bX, bY]= meshgrid([1:nCols]',[1:nRows]');

    v = reshape(bX,[],1);
    w = reshape(bY,[],1); 
    LfpStateSmooth = [];
    lLfp =length(LfpState(:,1));
    for k =1:lLfp
        display(k)
%         display(lLfp)
        map = LfpState(k,:);
        F=  scatteredInterpolant(v,w,map');
        bF = F(bX,bY);
        if GaussianFact
            bF = imgaussfilt(bF,GaussianFact);
        elseif MedFact
            bF = medfilt2(bF);
        end
        LfpState(k,:) = reshape(bF,[],1);
    end
end

% compute spectrogram on concatenated Lfp
% Parameters spectrogram
nFFT = 128;
if T.depth(indDB)==1 & T.SingleShank(indDB)==0
    nFFT = 256;
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
ChInd = input('channel index for events distribution');
AvPower = struct();
AvLfp = struct();
PowEvent=struct();
LfpEvent=struct();
LfpDistCh = struct();
PowerDist = struct();
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
    p = prctile(pks,95);
    figure()
    histogram(pks,100)
    vline(p)
    [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',p);
    pks = pks(1:end-2);
    locs = locs(1:end-2);
    
    [sPks, sLocs] = sort(pks);
    sLocs2 = locs;%
    %display peak amplitude vs power at peak
    if length(Channels)~=nCh
        for iCh = Channels
            figure()
            title(['Channel' num2str(iCh) '// IC' num2str(iSelectICs)])
            for iFreq  = 1:length(FreqBandsHigh)
                PowPeaks.(FreqBandsHigh{iFreq}) = pow.(FreqBandsHigh{iFreq})(floor(ICt(sLocs2)*FsSpec),iCh);

                subplot(floor(sqrt(length(FreqBandsHigh)))+1,floor(sqrt(length(FreqBandsHigh)))+1,iFreq)
                hist3([pks',PowPeaks.(FreqBandsHigh{iFreq})],'Nbins',[30,30],'CdataMode','auto')
    %             set(gca,'YScale','log')
                set(gca,'ColorScale','log')
                view(2)
                title(FreqBandsHigh{iFreq})
            end
        end

%         Chs = input('list of channels indices for power histogram comparison // or [] to skip')
        Chs = [];
        if Chs
            for DeltaT = [-0.1,0,0.1]
                for iFreq  = 1:length(FreqBandsHigh)
                    figure()
                    title(['Power' FreqBandsHigh{iFreq} '//time from peak' num2str(DeltaT) 's'])
                    hold on
                    for iCh = Chs
                        PowPeaksHist.(FreqBandsHigh{iFreq}) = pow.(FreqBandsHigh{iFreq})(floor((ICt(sLocs2)+DeltaT)*FsSpec),iCh);            
                        histogram(PowPeaksHist.(FreqBandsHigh{iFreq}),100)
    %                 PowPeaksHistNorm.(FreqBandsHigh{iFreq}) = PowPeaksHist.(FreqBandsHigh{iFreq})-mean(PowPeaksHist.(FreqBandsHigh{iFreq}))
    %                 subplot(floor(sqrt(length(FreqBandsHigh)))+1,floor(sqrt(length(FreqBandsHigh)))+1,iFreq)
                    end
                    legend(['Channel' num2str(Chs(1))],['Channel' num2str(Chs(2))])
                    if iFreq ==1
                        figure()
                        title(['LFP amp ' '//time from peak' num2str(DeltaT) 's'])
                        hold on
                        for iCh = Chs
                            FiltLfpState = ButFilter(LfpState(:,iCh),2,10/(FsHigh)*2,'low');
                        
                            histogram(FiltLfpState(floor((ICt(sLocs2)+DeltaT)*FsHigh)),100)
                        end
                        legend(['Channel' num2str(Chs(1))],['Channel' num2str(Chs(2))])
                    end

                end
            end
        end
    end
    
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
        SampsTriggerPower = floor(ICt(locs(iE))*FsSpec);
        display(iE)
        SampWidth = floor(Twindow*FsSpec);
        SampsWindow= [SampsTriggerPower-SampWidth:SampsTriggerPower+SampWidth];
    % define samples for Lfp  corresponding IC activation peak
        SampsTriggerLfp = floor(ICt(locs(iE))*FsHigh);
        SampWidthLfp = floor(Twindow*FsHigh);
        SampsWindowLfp= [SampsTriggerLfp-SampWidthLfp:SampsTriggerLfp+SampWidthLfp];

        if (locs(iE)-SampWidthIC)<1
            continue
        end
        AvIC = AvIC + ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC)';

        if iE ==3
            AvLfp.(FreqBandLow{1}) = LfpState(SampsWindowLfp,:);
        elseif iE == length(pks)
            LfpEvent.(FreqBandLow{1}) = LfpState(SampsWindowLfp,:);
            ICEvent = ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC);
        else
            AvLfp.(FreqBandLow{1}) = AvLfp.(FreqBandLow{1}) + LfpState(SampsWindowLfp,:);
        end

        for iFreq= 1:length(FreqBandsHigh)
%             zeroedPow = (pow.(FreqBandsHigh{iFreq})(SampsWindow,:)-repmat(mean(pow.(FreqBandsHigh{iFreq})(SampsWindow,:)),[length(SampsWindow),1])); 
%             NormalizedPow = (zeroedPow)./repmat(max(pow.(FreqBandsHigh{iFreq})(SampsWindow,:)),[length(SampsWindow),1]);%normalie power bands
            NormalizedPowCh = zscore(log10(pow.(FreqBandsHigh{iFreq})(SampsWindow,ChInd)));
            NormalizedPow = zscore(log10(pow.(FreqBandsHigh{iFreq})(SampsWindow,:)));
            if iE ==3
                AvPower.(FreqBandsHigh{iFreq}) = NormalizedPow;
                PowerDist.(FreqBandsHigh{iFreq})= NormalizedPowCh;
                if iFreq==1
                    LfpDistCh = LfpState(SampsWindowLfp,ChInd)-LfpState(SampsWindowLfp(1),ChInd);
                end
            elseif iE == length(pks)
                PowEvent.(FreqBandsHigh{iFreq}) = NormalizedPowCh;
            else
                AvPower.(FreqBandsHigh{iFreq}) = AvPower.(FreqBandsHigh{iFreq}) + NormalizedPow;%-repmat(pow.(FreqBandsHigh{iFreq})(SampsWindow(1),:),[length(pow.(FreqBandsHigh{iFreq})(SampsWindow(1),:)),1]);
                PowerDist.(FreqBandsHigh{iFreq})= [PowerDist.(FreqBandsHigh{iFreq}), NormalizedPowCh];
                if iFreq==1
                    LfpDistCh = [LfpDistCh,LfpState(SampsWindowLfp,ChInd)-LfpState(SampsWindowLfp(1),ChInd)];
                end
            end
        end
        if length(Channels)~=nCh
            for iCh = Channels
                
%                 zeroedEvent = repmat(abs(mean(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:))),[length(sStruct.(['Channel' num2str(iCh)])(SampsWindow,1)),1]);
%                 whiteningEvent = repmat(max(abs(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)-zeroedEvent)),[length(sStruct.(['Channel' num2str(iCh)])(SampsWindow,1)),1]);
%                 whitenedEvent = (sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)-zeroedEvent)./whiteningEvent;
                whitenedEvent = zscore(log10(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)));
                if iE ==3
                    AvSpec.(['Channel' num2str(iCh)]) = whitenedEvent;
                    LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]) = LfpState(SampsWindowLfp,iCh)-LfpState(SampsTriggerLfp,iCh);
                else
                    AvSpec.(['Channel' num2str(iCh)]) =  AvSpec.(['Channel' num2str(iCh)]) + whitenedEvent;
                    LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]) = [LfpDist.(FreqBandLow{1}).(['Ch' num2str(iCh)]),LfpState(SampsWindowLfp,iCh)-LfpState(SampsTriggerLfp,iCh)];
                end
            end
        end
    end

%     for iFreq= 1:length(FreqBandsHigh)
%         AvPower.(FreqBandsHigh{iFreq}) = (AvPower.(FreqBandsHigh{iFreq})-repmat(AvPower.(FreqBandsHigh{iFreq})(1,:),[length(AvPower.(FreqBandsHigh{iFreq})(:,1)),1]))/length(pks); 
%     end
    
    
    AvIC = AvIC/length(pks);
    %display the time evolution of AvPOwer for the different freq bands usign a
    %video
    out.AvIC = AvIC;
    out.ICt = [-SampWidthIC:SampWidthIC]/Fs;
    out.AvPower = AvPower;
    out.PowerT = [-SampWidth:SampWidth]/FsSpec;
    LfpT = [-SampWidthLfp:SampWidthLfp]/FsHigh;
%     save([OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-AvPower' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ').mat'],'out')

    figure()
    imagesc(PowerDist.(FreqBandsHigh{2}))
    clim([0.7*min(min(PowerDist.(FreqBandsHigh{2}))) 0.7*max(max(PowerDist.(FreqBandsHigh{2})))])
    hline(floor(length(SampsWindow)/2)+2)
    hline(floor(length(SampsWindow)/2)+7)

    figure()
    imagesc(LfpDistCh)
    clim([0.7*min(min(LfpDistCh)) 0.7*max(max(LfpDistCh))])
    hline(SampWidthLfp+10)
    hline(SampWidthLfp+200)
    figure()
    hold on
    scatter(LfpDistCh(SampWidthLfp+10,:),PowerDist.(FreqBandsHigh{2})(floor(length(SampsWindow)/2)+2,:))

    scatter(LfpDistCh(SampWidthLfp+200,:),PowerDist.(FreqBandsHigh{2})(floor(length(SampsWindow)/2)+7,:),'r')
%       
    if length(Channels)==nCh
        %video of AvPower
        OutputVideo = [OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-AvPower' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ')-MedFilt-95prctile-spec.avi'];
        Video_powerbands(out.ICt,out.AvIC',out.PowerT,out.AvPower, OutputVideo, T, indDB)
%         
        %video of AvLfp
        OutputVideo = [OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-AvLfp' Coupling 'PeaksSign(' num2str(SignCorr) ')-MedFilt-95prctile-spec.avi'];
        Video_powerbands(out.ICt,out.AvIC',LfpT,AvLfp, OutputVideo, T, indDB)        
        %video of power single event
        OutputVideo = [OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-PowerSingleEvent' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ')-MedFilt-95prctile-spec.avi'];
        Video_powerbands(out.ICt,ICEvent,out.PowerT,PowEvent, OutputVideo, T, indDB)
        %video of Lfp single event
        OutputVideo = [OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-LfpSingleEvent' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ')-MedFilt-95prctile-spec.avi'];
        Video_powerbands(out.ICt,ICEvent,LfpT,LfpEvent, OutputVideo, T, indDB)
    else
        for iCh = Channels
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
            
            figure()
            
            ax1 = subplot(4,1,1);
            plot(out.ICt,ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC))
            legend(['last event IC' num2str(iSelectICs)])
            ax2 = subplot(4,1,2);  
            plot(LfpT,LfpEvent.(FreqBandLow{1}))
            legend(['last event IC' num2str(iSelectICs)])
            
            ax3 = subplot(4,1,3);            
            indFreqSpec = find(AvSpec.f<=FMinSpec);
            whiteningEvent = repmat(max(abs(sStruct.(['Channel' num2str(iCh)])(SampsWindow,:))),[length(sStruct.(['Channel' num2str(iCh)])(SampsWindow,1)),1]);
            whitenedEvent = sStruct.(['Channel' num2str(iCh)])(SampsWindow,:)'./whiteningEvent';
            imagesc(out.PowerT,AvSpec.f(indFreqSpec),whitenedEvent(indFreqSpec,:))
  
            title(['Last event spec. channel' num2str(iCh) '-IC' num2str(iSelectICs)])
            ax4 = subplot(4,1,4);
            hold on
            for iFreq = 1:length(FreqBandsHigh)
                PowTrace = pow.(FreqBandsHigh{iFreq})(SampsWindow,iCh);
                plot(out.PowerT,PowTrace/max(abs(PowTrace)),'DisplayName',FreqBandsHigh{iFreq})
            end
            linkaxes([ax1,ax2,ax3,ax4],'x')
        end
    end
end




% BrainState = 'PerSWS';
% InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% out = TriggeredAmplitude(InputPath,OutputPath,T,ind(iFile),extension);

% BrainState = 'PerSWS';
% DownFact = 2;
% InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% out = TriggeredSpectrogram(InputPath,OutputPath,T,ind(iFile),extension,BrainState,{'isa'},{'spindle','beta','lowgamma','highgamma'},'.AC',DownFact,25);


% out.PowerR = struct();
% fn=fieldnames(out.AvPower);
% out.PowerR.(fn{1}) = out.AvPower.(fn{1});
% out.PowerR.(fn{2}) = out.AvPower.(fn{2});

% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/IC13-spec.avi';
% Video_powerbands(out.ICt,out.AvIC',out.PowerT,out.AvPower, OutputPath, T, ind(iFile))


%% signular spectrum analysis of signal amplitude
% BrainState = 'PerSWS';
% InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% out = SingularSpectrum(InputPath,OutputPath,T,ind(iFile),extension);

%% singular spectrum power
% % 
% BrainState = 'PerSWS';
% InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% out = SingularSpectrumPower(InputPath,OutputPath,T,ind(iFile),extension);


%%
% % generate video to display coeffSPA with dimensions nCh x 2*NsampsSPA 
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/SingularSpectrum-IC4.avi';
% Video_powerbands(out.SPAtime,out.AvIC(1:end-1)',out.SPAtime,out.SPA, OutputPath, T, ind(iFile))
