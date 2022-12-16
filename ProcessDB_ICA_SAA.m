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
Queries.DeviceName = 'B13289O14-DH3';%{'B13289O14-DH1','B13289O14-DH2','B13289O14-DH3','B13289O14-DH5','B13289O24-DH1SL7'};%{'B13289O14-DH1','B13289O14-DH2','B13289O14-DH3','B13289O14-DH5','B13289O24-DH1SL7'};%{-1,'B13289O14-DH1'}; % if cell array with -1 in first field, the opposite of the query is used


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
% %% Group statistics on pre-processed selected files
% % 
% % % extension = '.interp.lfp';
% % % 
% % OutPathStates = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/' T.RecordingId{ind(iFile)} '/States/'];
% % BrainStates2(OutPathStates,FileNameDB,T,ind(iFile),ind,'d');
% 
% % % extension ='.PCA.lfp'; %e.g. "smooth.interp.PCA.lfp"
% 
% % % ICA2(FileName,T,ind(iFile),OutPathStates,extension,'compute',50,'PerSWS',val);
% 
% %% for ICA, if val is passed, the val should be limited to one device and it will use concatenated data for that device
% % OutPathStates = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';
% % ICA2('',T,ind(iFile),OutPathStates,extension,'compute',100,'PerSWS',val,{'isa','delta','spindle'});
% %
% % OutPathStates = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';
% % ICA2('',T,ind(iFile),OutPathStates,extension,'display',100,'PerSWS',val,{'isa','delta','spindle'},'.ECoG');
% %
% % OutPathStates = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/'; %display works on device-specific "brain-state specific" groups file // it takes DC or AC coupling label from val list// the freq. bands should match//also downsampling taken accordingly//Number of IC must also match, otherwise either not all shown or display crashes
% % ICA2('',T,ind(iFile),OutPathStates,extension,'display',70,'PerSWS',val,{'theta','spindle','beta','lowgamma','highgamma'});%{'isa','delta','theta','spindle','beta','lowgamma','highgamma'}
% %
% % OutPathStates = '../../../data/DB/Files/';% path provided must point to generic /Files/ and val will be used to load animal-specific data (val must contain all animals evaluated)// both AC and DC files will be used in discrete electronics recordings//
% % ICA2('',T,ind(iFile),OutPathStates,extension,'Significance',100,'PerSWS',val,{'isa','delta','theta','spindle','beta','lowgamma','highgamma'});
% %  
% % OutPathStates = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
% % ICA2('',T,ind(iFile),OutPathStates,extension,'PhAmpCoupling',100,'PerSWS',val,{'delta'},'',{'spindle'},'.DC','.AC');
% 
% 
% %% compute triggered average
BrainState = 'PerSWS';
InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';


load(InputPath,'out')

TimeConstantSmooth = 0.1;

BrainState='PerSWS';
FreqBandLow={'delta'};
TimeScaleGlue=0.2;
indDB = ind(iFile);

nCh = T.NumCh(indDB);
Fs = T.Fs(indDB);
nRows = T.nRows(indDB);
nCols = T.nCols(indDB);
Kernel = floor(TimeConstantSmooth*Fs);

if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
    Coupling ='.AC';
elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
    Coupling ='.DC';
else
    Coupling ='';
end

Lfp=[];
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
        display(Fs,'*****************')
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


SelectedICs = 1:NumICs;
SignCorr = input('sign of peaks to detect ([1]/[-1]');
pksConcat=[];
locsConcat=[];
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
    locsConcat = [locsConcat,locs(1:end-2)];
end
locsConcat = sort(unique(locsConcat));

dLocsConcatDiff = diff(locsConcat);
locsConcat = locsConcat(find(dLocsConcatDiff>=4)); % keep only time points corresponding to IC peaks separated by more that 4/Fs
%build data matrix with IC activation values for all components at peaks
%detected in all components
Data=[];
for iSelectICs = SelectedICs %= input('index of IC to further process');
    ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:);
    Data = [Data;ICsig(locsConcat)];
end
Data = Data'; %Rows correspond to observations and columns to variables
% Data(find(Data<=0.2)) = 0;
counter =0;
SumDist=[];

for ik = 2:10:100
    display(ik)
    counter = counter +1;
    [idx,~,sumd] = kmeans(Data,ik,'MaxIter',1000);
    SumDist(counter) = sum(sumd);
end
figure()
plot(2:10:100,SumDist)


% [idx,~,sumd] = kmeans(Data,300,'MaxIter',1000);

kNum = length(unique(idx));
NumICs = length(A(1,:));
figure(210)
% MapKlustK = [];
% KDistSum = [];
% for ik = 1:kNum
%     title(['first 25 events in cluster k=' num2str(ik)])
%     indK = find(idx == ik);
%     KDistSum(ik) = sum(sumd(indK))/length(indK);
%     if length(indK)>=25
%         indK = indK(1:25);
%     else
%         indK = indK;
%     end
%     counter = 0;
%     display(ik)
%     for i = indK' 
%         counter = counter +1;
%         MapKlust = zeros(size(A(:,1)));
%         for iIC = 1:NumICs
%             MapKlust = Data(i,iIC)*A(:,iIC)+MapKlust;
%         end
%         MapKlustK = [MapKlustK,MapKlust];
%     end
%     CovMapKlustK = cov(MapKlustK);
%     CovSumK(ik) = sum(sum(CovMapKlustK))-sum(sum(eye(length(CovMapKlustK)).*CovMapKlustK));
%     figure(210)
%     hold on
%     scatter3(Data(find(idx==ik),1),Data(find(idx==ik),2),Data(find(idx==ik),3)) 
% end

% [SCovSumK,indCov] = sort(CovSumK);
% [SKDistSum,indKDist] = sort(KDistSum);
for ik = 1:kNum
    sumdNorm(ik) = sumd(ik)/length(find(idx==ik));
end
[SSumDist,indSumDist] = sort(sumdNorm);
for ik = 1:10
    figure()
    counter = 0;
    display(ik)
    indK = find(idx == indSumDist(ik));
    display(length(indK))
    display(SSumDist(ik))
    if length(indK)>=25
        indK = indK(1:25);
    else
        indK = indK;
    end
    for i = indK' 

        counter = counter +1;
        MapKlust = zeros(size(A(:,1)));
        for iIC = 1:NumICs
            MapKlust = Data(i,iIC)*A(:,iIC)+MapKlust;
        end

        subplot(floor(sqrt(25)),floor(sqrt(25)),counter)

        F=  scatteredInterpolant(v,w,MapKlust);
        bF = F(bX,bY);
        imagesc(bF)
        colorbar
    end
end
% 
% Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
% 
% figure()
% plot([1:length(Lfp(:,1))]/Fs,Lfp(:,28))
% %         AmpStackT = zeros(size(Lfp));
% %         TimeStackT = zeros(size(Lfp));
% %         AmpStackP = zeros(size(Lfp));
% %         TimeStackP = zeros(size(Lfp));
% 
% PosStackC = zeros(length(Lfp(:,1))-1,length(Lfp(1,:)));
% PosStackT = zeros(length(Lfp(:,1))-1,length(Lfp(1,:)));
% PosStackP = zeros(length(Lfp(:,1))-1,length(Lfp(1,:)));
% 
% LfpFilt = ButFilter(Lfp,2,[0.1 6]/(Fs)*2,'bandpass');
% LfpFiltDown = LfpFilt(1:2:end,:);
% %compute troughs, crossings and peaks for all channels
% %independently. The occurence of these events is summed over all
% %channels later and smoothened to find most coherent events across
% %multiple channels
% for iCh = 1:nCh
%     display(iCh,'Detecting waves in channel')
% 
%     LfpFiltCh = LfpFilt(:,iCh);
% 
%     LfpFiltSign = zeros(size(LfpFiltCh));
%     LfpFiltSign(find(LfpFiltCh>=0)) = 1;
%     LfpFiltSign(find(LfpFiltCh<0)) = 0;
% 
%     crossings = diff(LfpFiltSign); % yields -1 for upward crossings, and 1 for downward crossings
%     PosCrossings = find(crossings == 1); %negative to positive
%     dLfp = diff(LfpFiltCh);
%     dLfpSign = zeros(size(dLfp));
%     dLfpSign(find(dLfp>=0)) = 1;
%     dLfpSign(find(dLfp<0)) = 0;
%     CrosDLfpSign = diff(dLfpSign);
% 
%     PosStackC(:,iCh) =  crossings==1; % ones in samples where there is a negative to positive Lfp crossing
% %         figure()
% %         hold on
%     for i = 1:length(PosCrossings)
% 
%         if PosCrossings(i)>=length(LfpFiltCh)-10*Kernel | PosCrossings(i)<=10*Kernel
%             continue
%         end
% 
%         ii=0;
%         while CrosDLfpSign(PosCrossings(i)-ii) ~= 1
%             ii = ii+1;
%         end  
%         iii=0;
%         while CrosDLfpSign(PosCrossings(i)+iii) ~= -1
%             iii = iii+1;
%         end 
% %                 AmpStackT(i,iCh) = LfpFiltCh(PosCrossings(i)-ii);
% %                 TimeStackT(i,iCh) = ii*Fs;
% %                 AmpStackP(i,iCh) = LfpFiltCh(PosCrossings(i)+iii);
% %                 TimeStackP(i,iCh) = iii*Fs;
% 
%         PosStackT(PosCrossings(i)-ii,iCh) = 1;  
%         PosStackP(PosCrossings(i)+iii,iCh) = 1;
%     end
% end
% 
% SumC = mean(PosStackC')*nCh;
% Ccoocurrence = conv(SumC,ones([Kernel,1])); %number of channels/time (coocurrence) where a negative to possitive 0-crossing is detected
% dCcoocurrence = diff(Ccoocurrence); % calculate the derivative to find local maximums of crossings coocurrence
% dCcoocurrenceSign = zeros(size(dCcoocurrence)); 
% dCcoocurrenceSign(find(dCcoocurrence>=0)) = 1; %if slope is possitive assign 1 else leave 0
% CrosDCcoocurrenceSign = diff(dCcoocurrenceSign); %find transitions in the sign of sope or crossings coocurrence
% PosCArray = CrosDCcoocurrenceSign==-1; % get a 1 when there is a transition from positive to negative slope in crossings coocurrence
% PosGlobalCrossings = find(CrosDCcoocurrenceSign==1); % get indices where there is a "global" crossing (i.e. max in crossing coocurrence) 
% 
% SumT = mean(PosStackT')*nCh;
% Tcoocurrence = conv(SumT,ones([Kernel,1]));
% dTcoocurrence = diff(Tcoocurrence);
% dTcoocurrenceSign = zeros(size(dTcoocurrence));
% dTcoocurrenceSign(find(dTcoocurrence>=0)) = 1;
% CrosDLfpSign = diff(dTcoocurrenceSign);
% PosTArray = CrosDLfpSign==-1; % get a 1 when there is a transition from positive to negative slope in trough coocurrence
% 
% 
% SumP = mean(PosStackP')*nCh;
% Pcoocurrence = conv(SumP,ones([Kernel,1]));
% dPcoocurrence = diff(Pcoocurrence);
% dPcoocurrenceSign = zeros(size(dPcoocurrence));
% dPcoocurrenceSign(find(dPcoocurrence>=0)) = 1;
% CrosDLfpSign = diff(dPcoocurrenceSign);
% PosPArray = CrosDLfpSign==-1; % get a 1 when there is a transition from positive to negative slope in peak coocurrence
% 
% 
% %search "global" troughs and "global" peaks closest to "global"
% %crossings
% AmpGlobT = [];
% AmpGlobP = [];
% PosXGlobT = [];
% PosXGlobP = [];
% PosYGlobT = [];
% PosYGlobP = [];
% 
% PatchGlobT = [];
% PatchGlobP = [];
% 
% VarSpaceGlobT = [];
% VarSpaceGlobP = [];
% for i = 1:length(PosGlobalCrossings)
%     display(PosGlobalCrossings(i),'PosGlobalCrossings')
%     if PosGlobalCrossings(i)>=length(LfpFiltCh)-10*Kernel | PosGlobalCrossings(i)<=10*Kernel
%         continue
%     end
% 
%     ii=0;
%     while PosTArray(PosGlobalCrossings(i)-ii) ~= 1
%         ii = ii+1;
%     end  
%     iii=0;
%     while PosPArray(PosGlobalCrossings(i)+iii) ~=1
%         iii = iii+1;
%     end 
%     AmpGlobT = [AmpGlobT,min(LfpFilt(PosGlobalCrossings(i)-ii,:))];
% %             PatchGlobT = [PatchGlobT, find(LfpFilt(PosGlobalCrossings(i)-ii,:)<= -30)];
%     AmpGlobP = [AmpGlobP,max(LfpFilt(PosGlobalCrossings(i)+iii,:))];
% %             PatchGlobP = [PatchGlobP, find(LfpFilt(PosGlobalCrossings(i)+iii,:)>= -30)];
% 
%     PosXGlobT = [PosXGlobT,floor((find(LfpFilt(PosGlobalCrossings(i)-ii,:) == min(LfpFilt(PosGlobalCrossings(i)-ii,:)))-1)/nRows)+1];
%     PosYGlobT = [PosYGlobT,rem(find(LfpFilt(PosGlobalCrossings(i)-ii,:) == min(LfpFilt(PosGlobalCrossings(i)-ii,:)))-1,nCols)+1];
%     PosXGlobP = [PosXGlobP,floor((find(LfpFilt(PosGlobalCrossings(i)+iii,:) == max(LfpFilt(PosGlobalCrossings(i)+iii,:)))-1)/nRows)+1];
%     PosYGlobP = [PosYGlobP,rem(find(LfpFilt(PosGlobalCrossings(i)+iii,:) == max(LfpFilt(PosGlobalCrossings(i)+iii,:)))-1,nCols)+1];
%     VarSpaceGlobT = [VarSpaceGlobT, min(LfpFilt(PosGlobalCrossings(i)-ii,:))-max(LfpFilt(PosGlobalCrossings(i)-ii,:))];
%     VarSpaceGlobP = [VarSpaceGlobP, max(LfpFilt(PosGlobalCrossings(i)+iii,:))-min(LfpFilt(PosGlobalCrossings(i)+iii,:))];
% 
% %             PosIndexPeak = 
% %             PosIndexTrough = 
% end
% 
% Clu
% IndLargeGlobT = find(AmpGlobT>=30);



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
