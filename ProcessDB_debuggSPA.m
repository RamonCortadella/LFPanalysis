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
    
%     display(FileName)
%     factor = 2;
%     out = FindBadCh2(FileName,OutPathStates,T,ind(iFile), OutputPathTable,FileNameDB, 'compute', factor);
%     T = out.T;
% 
%     LoadMocap2(MocapFile,T,ind(iFile),OutPathStates,'display');
%
%     GmSmoothening2(FileName,T,ind(iFile),InputPath,extension);
%
%     SpectrogramSave2(FileName,OutPathStates,T,ind(iFile),'compute');
% 
%     out = BrainStates2(OutPathStates,FileNameDB,T,ind(iFile),ind,'display');
%     T = out.T;        
%      
%     LfpInterp = InterpolateProbes2(FileName, T,ind(iFile),InputPath);
%  
%     ICA2(strcat(FileName(1:end-4),extension),T,ind(iFile),extension,OutPathStates,val,'compute',30,'PerSWS');
%     
%     ICA2(FileName,T,ind(iFile),OutPathStates,extension,'compute',50,'PerSWS');
    
%     detectDelta(FileName,T,ind(iFile),[OutPathStates T.RecordingId{ind(iFile)}],'compute')
%
%     writetable(T,strcat(OutputPathTable,FileNameDB,'.xlsx'))

%     Video2(FileName,[OutPathStates 'Video.avi'], T,ind(iFile))
    
end
%% Group statistics on pre-processed selected files
% 
% % extension = '.interp.lfp';
% % 
% OutPathStates = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/' T.RecordingId{ind(iFile)} '/States/'];
% BrainStates2(OutPathStates,FileNameDB,T,ind(iFile),ind,'d');

% % extension ='.PCA.lfp'; %e.g. "smooth.interp.PCA.lfp"

% % ICA2(FileName,T,ind(iFile),OutPathStates,extension,'compute',50,'PerSWS',val);

%% for ICA, if val is passed, the val should be limited to one device and it will use concatenated data for that device
% OutPathStates = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';
% ICA2('',T,ind(iFile),OutPathStates,extension,'compute',100,'PerSWS',val,{'spindle','beta','lowgamma','highgamma'},'.ECoG');
%
% OutPathStates = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';
% ICA2('',T,ind(iFile),OutPathStates,extension,'display',100,'PerSWS',val,{'isa','delta','spindle'},'.ECoG');
%
% OutPathStates = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/'; %display works on device-specific "brain-state specific" groups file // it takes DC or AC coupling label from val list// the freq. bands should match//also downsampling taken accordingly//Number of IC must also match, otherwise either not all shown or display crashes
% ICA2('',T,ind(iFile),OutPathStates,extension,'display',70,'PerSWS',val,{'theta','spindle','beta','lowgamma','highgamma'});%{'isa','delta','theta','spindle','beta','lowgamma','highgamma'}
%
% OutPathStates = '../../../data/DB/Files/';% path provided must point to generic /Files/ and val will be used to load animal-specific data (val must contain all animals evaluated)// both AC and DC files will be used in discrete electronics recordings//
% ICA2('',T,ind(iFile),OutPathStates,extension,'Significance',100,'PerSWS',val,{'isa','delta','theta','spindle','beta','lowgamma','highgamma'});
%  
% OutPathStates = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% ICA2('',T,ind(iFile),OutPathStates,extension,'PhAmpCoupling',100,'PerSWS',val,{'delta'},'',{'spindle'},'.DC','.AC');


%% compute triggered average
% BrainState = 'PerSWS';
% InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
% OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% out = TriggeredSpectrogram(InputPath,OutputPath,T,ind(iFile),extension);

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
BrainState = 'PerSWS';
InputPath = strcat('../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ECoG','.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
OutputPath = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';
% out = SingularSpectrum(InputPath,OutputPath,T,ind(iFile),extension);
BrainState = 'PerSWS';
FreqBandsIC = {'delta'};
FreqBandsSPA = {'delta'};
SPACoupling = '.DC';
Twindow = 0.3;
TimeScaleGlue =0.2;

indDB = ind(iFile);
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

display(fieldnames(out))

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
A = out.(FreqBandsIC{1}).A;
display(length(A(1,:)),'Index of IC in the displayed computation: CHECK FOR INCONSISTENCIES!!!')

%display components map
[bX, bY]= meshgrid([1:nCols]',[1:nRows]');
v = reshape(bX,[],1);
w = reshape(bY,[],1);
figure()
title(FreqBandsIC{1})
NumICs = length(A(1,:));
for i = 1:NumICs
    subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
    F=  scatteredInterpolant(v,w,A(:,i));
    bF = F(bX,bY);
    imagesc(bF)
    colorbar

    title(['index=' num2str(i)])
end


InputPathState = '/storage2/ramon/data/DB/Files/';
Lfp=[];
val = out.val;
if strcmp(Coupling,'.DC') & strcmp(SPACoupling,'.AC')
   FsSPA = T.Fs(indDB-1);
else
   FsSPA = T.Fs(indDB);
end
SelectedPeriods = [];
for i = 1:length(val)

    Tend = 0;
    if contains(extension,'interp') %force use of interpolated files
        dn = split(val{1},'-');
        Device = [dn{1} '-' dn{2}];
        if contains(T.FileName(indDB),'ECoG')
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.ECoG' SPACoupling extension];
        elseif  contains(T.FileName(indDB),'Depth')
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.Depth' SPACoupling extension];
        else
            FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} SPACoupling extension];
        end

        display(FileName, 'loading file for Lfp concatenation across recordings')
        LfpT = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
        LfpT = ButFilter(LfpT,2,fB.(FreqBandsSPA{1})/(FsSPA)*2,'bandpass');
            
        szLfp = size(Lfp);
        Tend = szLfp(1)/FsSPA;

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



%compute singlular spectrum analysis (SPA)
% concatenate triggered events on target SPA band (i.e.
szLfpState = size(LfpState);
NsampsSPA = floor(Twindow*Fs);
display(NsampsSPA)
LfpSPA = zeros(szLfpState(1),szLfpState(2)*2*NsampsSPA);
for iCh = 1:nCh
    for iSPA = 1:NsampsSPA
        display(iCh)
        display(iSPA)
        LfpSPA(1:szLfpState(1)-iSPA+1,(2*(iCh-1)+1)*NsampsSPA-iSPA+1) = LfpState(iSPA:szLfpState(1),iCh);

        LfpSPA(iSPA:szLfpState(1),(2*(iCh-1)+1)*NsampsSPA+iSPA) = LfpState(1:szLfpState(1)-iSPA+1,iCh);
    end
end

% for SignCorr = [-1,1]
%     for SelectICs = 1:NumICs %= input('index of IC to further process');
SelectICs =input('IC for SPA analysis');
SignCorr= input('sign activation peaks');
%input('multiply activation by ([1]/[-1]) // inverse sign if loadings are not correct (i.e. negative for delta waves)');
% 
%     SelectICs = input('index of IC to further process');
%     SignCorr = input('multiply activation by ([1]/[-1]) // inverse sign if loadings are not correct (i.e. negative for delta waves)');
% find selected IC peaks

ICsig = out.(FreqBandsIC{1}).ICAsig(SelectICs,:)*SignCorr;
ICt = [1:length(ICsig)]/Fs;


Threshold = 0.1;%input('input the threshold for detection of peaks (it must be positive)');

[pks,locs] = findpeaks(ICsig, 'MinPeakHeight',Threshold);
% if length(pks)<=50
%     continue
%     display(['No pks above threshold for IC' num2str(iSelectICs)])
% end
p = prctile(pks,90);
figure()
histogram(pks,100)
vline(p)
[pks,locs] = findpeaks(ICsig, 'MinPeakHeight',p);
pks = pks(1:end-2);
locs = locs(1:end-2);

% if length(pks)<=50
%     continue
% end

% use val of loaded ICs to compute spectrograms


LfpAvPeaks = zeros([2*floor(Twindow*Fs),length(LfpState(1,:))]);
LfpAvSPA = zeros([2*floor(Twindow*FsSPA)+1,length(LfpSPA(1,:))]);
SampWidthIC = floor(Twindow*Fs);
AvIC = zeros([length(-SampWidthIC:SampWidthIC),1]);
LfpEventsSPA = [];
for iE = 1:length(pks)

    AvIC = AvIC + ICsig(locs(iE)-SampWidthIC:locs(iE)+SampWidthIC)';

    SampsTriggerPower = floor(ICt(locs(iE))*FsSPA);

    sampswindow = [SampsTriggerPower-floor(Twindow*FsSPA)+1:SampsTriggerPower+floor(Twindow*FsSPA)];
    display(iE)
    LfpEventsSPA((iE-1)*length(sampswindow)+1:iE*length(sampswindow), :) = LfpSPA(sampswindow,:);
%     LfpAvSPA = LfpSPA(sampswindow,:)+LfpAvSPA;
    LfpAvPeaks = LfpState(sampswindow,:) + LfpAvPeaks;
end
AvIC = AvIC/length(pks);

NumPCs=10;
%compute PCA and recover temporal progression.
% [coeff, score,latent,tsquared,explained,~] = pca(LfpEventsSPA,'NumComponents', NumPCs);
% coeffFact = rotatefactors(coeff,'Method','varimax','Maxit',1000);
[LU,LR,FSr,VT] = erpPCA2(LfpEventsSPA,10);
coeffSPA = struct();
out = struct();
for iPCA = 1:NumPCs
    coeffSPA.(['PC' num2str(iPCA)]) = reshape(coeffFact(:,iPCA),[2*NsampsSPA,nCh]);
    Var = sum(reshape(coeffSPA.(['PC' num2str(iPCA)]),[],1).*reshape(coeffSPA.(['PC' num2str(iPCA)]),[],1));
    VarRCs(iPCA) = Var;
%     if iPCA==1
%         coeffSPA.(['PC' num2str(1)]) = zeros(size(reshape(coeff(:,iPCA),[2*NsampsSPA,nCh])));
%     end
%     coeffSPA.(['PC' num2str(1)]) = coeffSPA.(['PC' num2str(1)])+reshape(coeff(:,iPCA),[2*NsampsSPA,nCh]);
%     display(iPCA)
end
[SVarRCs,ind] = sort(VarRCs);
SelectedCoeffSPA = struct();
for iSPA =  1:4
    SelectedCoeffSPA.(['SPAindex' num2str(iSPA)]) = coeffSPA.(['PC' num2str(ind(end-iSPA))]);
    display(iSPA)
end
out.SPAtime = [1:length(SelectedCoeffSPA.(['SPAindex' num2str(1)])(:,1))]/Fs;
out.SPA = SelectedCoeffSPA;
out.AvIC = AvIC;
out.ICt = [-SampWidthIC:SampWidthIC]/Fs;

%         save([OutputPath FreqBandsIC{1} '-TimeEmbeddedIC' num2str(SelectICs) '-PeaksSign(' num2str(SignCorr) ').mat'],'out')

OutputVideo = [OutputPath FreqBandsIC{1} '-TimeEmbeddedIC' num2str(SelectICs) '-PeaksSign(' num2str(SignCorr) ')-1000PCs-erpPCA-spec.avi'];
%     OutputPath = ['../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/' FreqBandsLow{1} '-AvIC' num2str(SelectICs) '-AvPower' HighFreqCoupling '-spec.avi'];
Video_powerbands(out.SPAtime,out.AvIC(1:end-1)',out.SPAtime,out.SPA, OutputVideo, T, indDB)


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
