function out = TriggeredAmplitude(InputPath,OutputPath,T,indDB,extension,varargin)
%input path specifies the path+filename of the .ICA.mat file containing the
%data of interest. indDB is taken from the searchDB which must contain the
%data used for ICA in order to extract the right metadata (i.e. Fs, nRows,
%etc). 
[BrainState, FreqBandLow, Twindow,TimeScaleGlue] = DefaultArgs(varargin,{'PerSWS',{'delta'},1,0.2}); %'isa','delta','theta','spindle','beta','lowgamma','highgamma'
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


SelectedPeriods = [];
for i = 1:length(val)

    Tend = 0;
    if contains(extension,'interp') %force use of interpolated files
        dn = split(val{1},'-');
        Device = [dn{1} '-' dn{2}];
        if contains(T.FileName(indDB),'ECoG')
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.ECoG' Coupling extension];
        elseif  contains(T.FileName(indDB),'Depth')
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.Depth' Coupling extension];
        else
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} Coupling extension];
        end

        display(FileName, 'loading file for Lfp concatenation across recordings')
        LfpT = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';

        szLfp = size(Lfp);
        Tend = szLfp(1)/Fs;
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
    display(floor(Fs*SelectedPeriods(i,:)), 'min and max samp Period')
    LfpStateTemp = Lfp(floor(Fs*SelectedPeriods(i,1)+1):floor(Fs*SelectedPeriods(i,2))-1,:); 

    WindowTemp = ones([length(LfpStateTemp(:,1)),1]);
    Window1 = sigmoid([1:length(LfpStateTemp(:,1))]/Fs-10*TimeScaleGlue,TimeScaleGlue);
    Window= WindowTemp.*Window1'.*flip(Window1');
    Window = repmat(Window,1,nCh);
    display(size(LfpStateTemp))
    display(size(Window))
    LfpStateTemp = LfpStateTemp.*Window;


    LfpState = [LfpState; LfpStateTemp];
end
if strcmp(Coupling,'.AC')
    display([49 51]/(Fs)*2)
    display(Fs)
    LfpState = ButFilter(LfpState,2,[49 51]/(Fs)*2,'stop');
    LfpState = ButFilter(LfpState,2,[99 101]/(Fs)*2,'stop');
%     LfpState = ButFilter(LfpState,2,[149 151]/(FsHigh)*2,'stop');
end

SelectICs = input('indices list of ICs to further process');
display(SelectICs)
for SignCorr = [1,-1]
    for iSelectIC = SelectICs %
        %input('multiply activation by ([1]/[-1]) // inverse sign if loadings are not correct (i.e. negative for delta waves)');
        % find selected IC peaks
        display(iSelectIC)
        ICsig = out.(FreqBandLow{1}).ICAsig(iSelectIC,:)*SignCorr;
        ICt = [1:length(ICsig)]/Fs;
    % 
    %     figure()
    %     plot(ICt,ICsig)

        Threshold = 2;%input('input the threshold for detection of peaks (it must be positive)');

        [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',Threshold);
        display(length(pks))
        if length(pks)<=50
            continue
            display(['No pks above threshold for IC' num2str(iSelectIC)])
        end
        pks = pks(1:end-2);
        locs = locs(1:end-2);
        % use val of loaded ICs to compute spectrograms


        % compute triggered average of the power bands
        SampWidthIC = floor(Twindow*Fs);
        AvIC = zeros([length(-SampWidthIC:SampWidthIC),1]);
        for iE = 3:length(pks)
            SampsTrigger = floor(ICt(locs(iE))*Fs);
%             display(iE)
            SampWidth = floor(Twindow*Fs);
            SampsWindow= [SampsTrigger-SampWidth:SampsTrigger+SampWidth];
            if (locs(iE)-SampWidthIC)<1
                continue
            end
            AvIC = AvIC + ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC)';
%             display(size(LfpState))
%             display(SampsWindow(1))
%             display(SampsWindow(end))
            if iE ==3
                AvLfp.(FreqBandLow{1}) = LfpState(SampsWindow,:);
            else
                AvLfp.(FreqBandLow{1}) = AvLfp.(FreqBandLow{1}) + LfpState(SampsWindow,:);
            end

        end

        AvIC = AvIC/(length(pks)-3);
        %display the time evolution of AvPOwer for the different freq bands usign a
        %video
        out.AvIC = AvIC;
        out.ICt = [-SampWidthIC:SampWidthIC]/Fs;
        out.AvLfp = AvLfp;
        out.LfpT = [-SampWidth:SampWidth]/Fs;
        save([OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectIC) '-AvLfp' Coupling 'PeaksSign(' num2str(SignCorr) ').mat'],'out')

        OutputVideo = [OutputPath FreqBandLow{1} '-AvIC' num2str(iSelectIC) '-AvLfp' Coupling 'PeaksSign(' num2str(SignCorr) ')-spec.avi'];
    %     OutputPath = ['../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/' FreqBandsLow{1} '-AvIC' num2str(SelectICs) '-AvPower' HighFreqCoupling '-spec.avi'];
        Video_powerbands(out.ICt,out.AvIC',out.LfpT,out.AvLfp, OutputVideo, T, indDB)
    end
end

 