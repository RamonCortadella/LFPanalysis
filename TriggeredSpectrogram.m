function out = TriggeredSpectrogram(InputPath,OutputPath,T,indDB,extension,varargin)
%input path specifies the path+filename of the .ICA.mat file containing the
%data of interest. indDB is taken from the searchDB which must contain the
%data used for ICA in order to extract the right metadata (i.e. Fs, nRows,
%etc). 
[BrainState, FreqBandLow, FreqBandsHigh, HighFreqCoupling,DownFact, Twindow,TimeScaleGlue] = DefaultArgs(varargin,{'PerSWS',{'delta'},{'spindle','beta','lowgamma','highgamma'},'.AC',[],1,0.2}); %'isa','delta','theta','spindle','beta','lowgamma','highgamma'
nCh = T.NumCh(indDB);
nRows = T.nRows(indDB);
nCols = T.nCols(indDB);
Fs = T.Fs(indDB);

load(InputPath,'out')
fB.isa = [0.02 0.5];
fB.delta = [0.5 5];
fB.theta = [5 9];
fB.spindle = [9 20];
fB.beta = [20 30];
fB.lowgamma = [30 50];
fB.highgamma = [50 100];


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
counter = 0;
Channels = input('channel indices to compute spectrogram from // if string all it computes spectrogram for all channels and produces video');
if strcmp(Channels,'all')
    Channels = [1:nCh];
end
for iCh = Channels
    counter = counter+1;
    display(iCh,'Computing spectrogram for Channel')

%     display(size(LfpState))
%     display(LfpState(1:4,1:4))
    [s,f,t] = mtcsglong(LfpState(:,iCh),nFFT,FsHigh,nFFT,3*nFFT/4);%,1024);
    display(size(s))
    display(size(f))
    display(f)
    display(s(1:10,1:10))
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
        pow.(FreqBandsHigh{iFreq})(:,iCh) = sum(s(:,indFreq),2);
    end
end


SelectedICs = input('indices of selected ICs to use for triggered average');
SignCorr = input('sign of peaks to detect ([1]/[-1]');
for iSelectICs = SelectedICs %= input('index of IC to further process');
    %input('multiply activation by ([1]/[-1]) // inverse sign if loadings are not correct (i.e. negative for delta waves)');
    % find selected IC peaks

    ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:)*SignCorr;
    ICt = [1:length(ICsig)]/Fs;

    Threshold = 2;%input('input the threshold for detection of peaks (it must be positive)');

    [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',Threshold);
    if length(pks)<=50
        continue
        display(['No pks above threshold for IC' num2str(iSelectICs)])
    end
    pks = pks(1:end-2);
    locs = locs(1:end-2);
    % use val of loaded ICs to compute spectrograms


    % compute triggered average of the power bands
    SampWidthIC = floor(Twindow*Fs);
    AvIC = zeros([length(-SampWidthIC:SampWidthIC),1]);
    for iE = 3:length(pks)
        SampsTriggerPower = floor(ICt(locs(iE))*FsSpec);
        display(iE)
        SampWidth = floor(Twindow*FsSpec);
        SampsWindow= [SampsTriggerPower-SampWidth:SampsTriggerPower+SampWidth];
        if (locs(iE)-SampWidthIC)<1
            continue
        end
        AvIC = AvIC + ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC)';
        for iFreq= 1:length(FreqBandsHigh)
            if iE ==3
                AvPower.(FreqBandsHigh{iFreq}) = pow.(FreqBandsHigh{iFreq})(SampsWindow,:);
            else
                AvPower.(FreqBandsHigh{iFreq}) = AvPower.(FreqBandsHigh{iFreq}) + pow.(FreqBandsHigh{iFreq})(SampsWindow,:);
            end
        end
        if ~strcmp(Channels,'all')
            for iCh = Channels
                if iE ==3
                    AvSpec.(['Channel' num2str(iCh)]) = sStruct.(['Channel' num2str(iCh)])(SampsWindow,:);
                else
                    AvSpec.(['Channel' num2str(iCh)]) =  AvSpec.(['Channel' num2str(iCh)]) + sStruct.(['Channel' num2str(iCh)])(SampsWindow,:);
                end
            end
        end
        
    end

    for iFreq= 1:length(FreqBandsHigh)
        AvPower.(FreqBandsHigh{iFreq}) = (AvPower.(FreqBandsHigh{iFreq})-repmat(AvPower.(FreqBandsHigh{iFreq})(1,:),[length(AvPower.(FreqBandsHigh{iFreq})(:,1)),1]))/length(pks); 
    end
    
    
    AvIC = AvIC/length(pks);
    %display the time evolution of AvPOwer for the different freq bands usign a
    %video
    out.AvIC = AvIC;
    out.ICt = [-SampWidthIC:SampWidthIC]/Fs;
    out.AvPower = AvPower;
    out.PowerT = [-SampWidth:SampWidth]/FsSpec;
%     save([OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-AvPower' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ').mat'],'out')

    if strcmp(Channels,'all')
        OutputVideo = [OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectICs) '-AvPower' HighFreqCoupling 'PeaksSign(' num2str(SignCorr) ')-spec.avi'];
        Video_powerbands(out.ICt,out.AvIC',out.PowerT,out.AvPower, OutputVideo, T, indDB)
    else
        for iCh = Channels
            figure()
            ax1 = subplot(3,1,1);
            plot(out.ICt,AvIC)
            legend(['averaged IC' num2str(iSelectICs)])
            ax2 = subplot(3,1,2);
            imagesc(out.PowerT,AvSpec.f,AvSpec.(['Channel' num2str(iCh)]))
            title(['Av. Spec. channel' num2str(iCh) '-IC' num2str(iSelectICs)])
            ax3 = subplot(3,1,3);
            hold on
            for iFreq = 1:length(FreqBandsHigh)
                plot(out.PowerT,out.AvPower.(FreqBandsHigh{iFreq})(:,iCh),'DisplayName',FreqBandsHigh{iFreq})
            end
            linkaxes([ax1,ax2,ax3],'x')
            
            figure()
            
            ax1 = subplot(3,1,1);
            plot(out.ICt,ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC))
            legend(['last event IC' num2str(iSelectICs)])
            ax2 = subplot(3,1,2);
            imagesc(out.PowerT,AvSpec.f,sStruct.(['Channel' num2str(iCh)])(SampsWindow,:))
            title(['Last event spec. channel' num2str(iCh) '-IC' num2str(iSelectICs)])
            ax3 = subplot(3,1,3);
            hold on
            for iFreq = 1:length(FreqBandsHigh)
                plot(out.PowerT,pow.(FreqBandsHigh{iFreq})(SampsWindow,iCh),'DisplayName',FreqBandsHigh{iFreq})
            end
            linkaxes([ax1,ax2,ax3],'x')
        end
    end
end

 