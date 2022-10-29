MetaDataBase = '/storage2/ramon/data/DB/RecordingsDB-Temp4.xlsx';
close('all')
%% Functions Settings
%--- bad channels identification
BadChannels = true;
factor = 2.5;

%--- smoothening
smooth = 0;
iterSmoothening = 3;

%--- bad channels interpolation
Interpolation = false;

%--- ICA decomposition
ICA=false;
numSources=100;
%--- splindle factor analysis (sFA) 
sFA = false;

%--- save?
SaveResult = false;

%--- process mocap to find motor states
MotorClass = false;
OutPathMotor = '../../../data/DB/States/';
InputPathMotor = '../../../data/DB/Mocap/';
vThreshold = 0.05; %velocity considered active
duration = 5;%minimum duration of quite state
Fs = 180; %frames per second mocap
FallSleepTime =20; %Fall sleep time must be smaller that MinSleepTime, if a quite period is longer than MinSleepTime, the initial FallSleepTime will be Quite and the rest sleep
MinSleepTime = 40;

%--- compute spectrograms
SpecSave = false;
OutPathSpectrogram = '../../../data/DB/States/';

%--- Compute Brain State Classification
BrainClass = false;
OutPathBrainStates = '../../../data/DB/States/';
MinThetaDeltaRatio = 0.8;
MinTimeRatio = 10; %minimum duration of theta ON state
MinSpindleDeltaRatio = 4.5;
MinDurationHVS_sec = 1;
ThetaBand = [[3,5];[6,9]];
HVSband = [[0,4];[11 30]];

%% define queries
Queries={};
% Queries.depth=1;
% Queries.coverage=1;
Queries.flag=1;
Queries.CoupledAC=1;
Queries.mocap=1;
Queries.SingleShank=0;
Queries.LostMocapSamples=0;
% Queries.depth=0;
fn = fieldnames(Queries);

%% load MetaDataBase  and apply queries/return selected files
T = readtable(MetaDataBase);
T = rmmissing(T);

for i = 1:length(fn)
    ind = find(table2array(T(:,fn(i)))==Queries.(fn{i}));
    st = table2array(T(ind,'recordingId'));
    if i >= 2
        [st,~] = intersect(st,st0,'stable');
        [ind,~] = intersect(ind,ind0);
    end
    st0=st;
    ind0=ind;
end
val = st;

%% apply functions on loaded recordings
for iFile = 1:length(val)
    if iFile<=12
        continue
    end
    if smooth==1 | BadChannels == true | SpecSave == true| BrainClass == true | ICA==true
        nCh = T.NumCh(ind(iFile));
        InputPath = '../../../data/DB/MappedData/';
        FileName = [InputPath, val{iFile}];
        Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';

        smooth = T.Smooth(ind(iFile));
        nRows = T.nRows(ind(iFile));
        nCols = T.nCols(ind(iFile));
        Chs = [1:nRows*nCols];
        Fs = T.Fs(ind(iFile));
    end

%---------------    
    if BadChannels == true
        OutputPath = '../../../data/DB/';
        if smooth==1
            factor = factor+3;
        end
        [badChs,RMSarray] = FindBadCh(Lfp, Fs, nRows, nCols, factor);
       
        T.badChannels(ind(iFile))= {mat2str(badChs)};
        writetable(T,strcat(OutputPath,'RecordingsDB-Temp4','.xlsx'))
        
    end
%---------------    
    if smooth == 1
        badindex = str2num(T.badChannels{ind(iFile)});
        Lfp = GmSmoothening(Lfp, [0.5, 10], Fs, nRows, nCols, Chs, iterSmoothening,badindex);
        %repeat bad ch search
        OutputPath = '../../../data/DB/';
        
        
        factor = factor-3;
        [badChs,RMSarray] = FindBadCh(Lfp, Fs, nRows, nCols, factor);
        T.badChannels(ind(iFile))= {mat2str(badChs)};
        writetable(T,strcat(OutputPath,'RecordingsDB-Temp4','.xlsx'))
    end   
%---------------    
    if Interpolation == true
        badindex = str2num(T.badChannels{ind(iFile)});
        nRows = T.nRows(ind(iFile));
        nCols = T.nCols(ind(iFile));
        depth = T.depth(ind(iFile));
        SingleShank = T.SingleShank(ind(iFile));
        Lfp = InterpolateProbes(Lfp, badindex,nRows,nCols,depth,SingleShank);
    end    
%---------------
    if MotorClass == true     
        fn = split(val{iFile},'-');
        rn = split(fn{3},'.');
        MocapFile = strcat(fn{1},'-',fn{2},'-rec',rn{1}(4:end),'.csv');
        FsMocap = T.FsMocap(ind(iFile));
        Fs = T.Fs(ind(iFile));
        [Periods,Statemap,Broken] = LoadMocap(MocapFile,InputPathMotor,OutPathMotor,vThreshold,duration,Fs,FsMocap,FallSleepTime,MinSleepTime);
        T.LostMocapSamples(ind(iFile)) = Broken;
        
        writetable(T,strcat(OutputPath,'RecordingsDB-Temp4','.xlsx'))
        
    end
%---------------        
    if SpecSave == true
        ChHVS = T.ChHVS(ind(iFile));
        ChTheta = T.ChTheta(ind(iFile));
        sLfp = size(Lfp);
        if sLfp(2)==1024
            ChHVS=ChHVS*2;
            ChTheta = ChTheta*2;
            nRows=32;
        end
        Lfp(isnan(Lfp))=0;
        SpectrogramSave(Lfp,Fs,ChHVS,ChTheta,OutPathSpectrogram,val{iFile},nRows);
    end
%---------------
    if BrainClass == true
        fn = split(val{iFile},'-');
        rn = split(fn{3},'.');
        MotorState = load(strcat(OutPathBrainStates,fn{1},'-',fn{2},'-',rn{1},'-MotorState.mat'));
        Periods = load(strcat(OutPathBrainStates,fn{1},'-',fn{2},'-',rn{1},'-Periods.mat'));
        MocapDuration = load(strcat(OutPathBrainStates,fn{1},'-',fn{2},'-',rn{1},'-MocapDuration.mat'));
        SpecHVS = load(strcat(OutPathBrainStates,fn{1},'-',fn{2},'-',rn{1},'-SpecHVS.mat'));
        SpecStates = load(strcat(OutPathBrainStates,fn{1},'-',fn{2},'-',rn{1},'-SpecStates.mat'));
        Trigger = str2num(T.TriggerTimes{ind(iFile)});
        FsMocap = T.FsMocap(ind(iFile));
        MinTimeRatio = 20;
    
        ReturnTrigger = BrainStates(OutPathBrainStates,val{iFile},MotorState,SpecHVS,SpecStates,Trigger,Periods,MinThetaDeltaRatio,MinTimeRatio,MinSpindleDeltaRatio,MinDurationHVS_sec,ThetaBand,HVSband);

        T.TimeExpectedTrigger(ind(iFile))= ReturnTrigger(1);
        T.TimeEndTrigger(ind(iFile))= ReturnTrigger(2);
        T.SynchMocap(ind(iFile))= ReturnTrigger(3);
        OutputPath = '../../../data/DB/';
        writetable(T,strcat(OutputPath,'RecordingsDB-Temp4','.xlsx'))
    end
%---------------
    if ICA == true
        [~, A, W]=fastica(Lfp','numOfIC',12);
        
        nRows = T.nRows(ind(iFile));
        nCols = T.nCols(ind(iFile));
        [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
        v = reshape(bX,[],1);
        w = reshape(bY,[],1);
        for i = 1:12
            subplot(10,10,i)
            F=  scatteredInterpolant(v,w,A(:,i));
            bF = F(bX,bY);
            pcolor(bF)
            caxis([-0.3 0.3])
            colorbar
        end

    end
    
%---------------
    if sFA == true
        fn = split(val{iFile},'-');
        rn = split(fn{3},'.');
        SpecHVS = load(strcat(OutPathBrainStates,fn{1},'-',fn{2},'-',rn{1},'-SpecHVS.mat'));
%         [coeff, Data_PCA,latent,tsquared,explained,mu] = pca(LfpGeom, 'NumComponents', q);
% 
%         B = ROTATEFACTORS(coeff, 'Method','varimax');
%         nRows = T.nRows(ind(iFile));
%         nCols = T.nCols(ind(iFile));
%         [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
%         v = reshape(bX,[],1);
%         w = reshape(bY,[],1);
%         for i = 1:12
%             subplot(10,10,i)
%             F=  scatteredInterpolant(v,w,A(:,i));
%             bF = F(bX,bY);
%             pcolor(bF)
%             caxis([-0.3 0.3])
%             colorbar
%         end

    end
%--------------- 
    if SaveResult == true
        
        OutputPath = '../../../data/DB/ProcessedFiles/';
        SaveBinary(strcat(OutputPath,val{iFile}(1:end-4),'Interp','.lfp'), Lfp)
    end    


%     LfpCorr = ButFilter(LfpCorr,4,[0.5, 30]/(Fs/2),'bandpass');

%     [xlsxFields,indRec,values] = FindBadCh(Lfp, Fs);
%     T.badChannels{ind(iFile)}=mat2str(badChs);%generalize to different
%     outputs from my fun
    

end
%% write a log file with queries and functions used on which date