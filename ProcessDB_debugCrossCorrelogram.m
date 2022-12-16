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
InputPath = strcat('../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ECoG','.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
OutputPath = '../../../data/DB/Files/B13289O14O23-DH3SL5/GroupsAnalysis/';

BrainState='PerSWS';
FreqBandLow={'delta'};
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

SelectedICs = input('indices of selected ICs to use for triggered average');
SignCorr = input('sign of peaks to detect ([1]/[-1]/both)');
EventsTrain = struct();
RandEventsTrain= struct();
T= [];
G= [];
TimeBin = 1;
MaxLag = 50;

sampsBin = floor(TimeBin*Fs);
if sampsBin ==0
    sampsBin=1;
end

for iSelectICs = SelectedICs %= input('index of IC to further process');
    %input('multiply activation by ([1]/[-1]) // inverse sign if loadings are not correct (i.e. negative for delta waves)');
    % find selected IC peaks
    ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:);
    if strcmp(SignCorr,'both')
        ICsig = abs(ICsig);
    else
        if iSelectICs  == 5
            ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:)*1;
        elseif iSelectICs  == 6
            ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:)*(-1);  
        else
            ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:)*SignCorr;
        end
    end
    ICt = [1:length(ICsig)]/Fs;

    Threshold = 0.1;%input('input the threshold for detection of peaks (it must be positive)');

    [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',Threshold);
    if length(pks)<=50
        continue
        display(['No pks above threshold for IC' num2str(iSelectICs)])
    end
    p = prctile(pks,90);
    
    [pks,locs] = findpeaks(ICsig, 'MinPeakHeight',p);
    pks = pks(1:end-2);
    locs = locs(1:end-2);
    if iSelectICs == 1
        T= locs;
        G = iSelectICs*ones(size(locs));
        size(T)
    else
        T= [T,locs];
        G= [G,iSelectICs*ones(size(locs))];
        size(T)
    end
    EventsTrain.(['IC' num2str(iSelectICs)]) = zeros([1,floor(length(ICsig)/sampsBin)]);
    for iSamps = 1:floor(length(ICsig)/sampsBin)
        EventsTrain.(['IC' num2str(iSelectICs)])(iSamps) = any((locs>=(iSamps-1)*sampsBin+1) & locs<(iSamps)*sampsBin);
        
    end
    display(length(locs))
    RandEventsTrain.(['IC' num2str(iSelectICs)]) = EventsTrain.(['IC' num2str(iSelectICs)])(randperm(length(EventsTrain.(['IC' num2str(iSelectICs)]))));   
end
CCG(T, G, sampsBin, floor(MaxLagSamp/2), Fs)

FsLags = Fs/sampsBin;
MaxLagSamp = floor(MaxLag*FsLags);

for iSelectICs = SelectedICs
    figure()
    counter2=0;
    for iSelectICs2 = SelectedICs
        display(counter2)
        counter2=counter2+1;
        [c,lags] = xcorr(EventsTrain.(['IC' num2str(iSelectICs)]),EventsTrain.(['IC' num2str(iSelectICs2)]),MaxLagSamp);
        lags = lags/(FsLags);
        [cRand,lagsRand] = xcorr(RandEventsTrain.(['IC' num2str(iSelectICs)]),RandEventsTrain.(['IC' num2str(iSelectICs2)]),MaxLagSamp);
        lagsRand = lagsRand/(FsLags);
        
        subplot(floor(sqrt(length(SelectedICs)))+1,floor(sqrt(length(SelectedICs)))+1,counter2)
        if iSelectICs == iSelectICs2
            Color='r';
        else
            Color='k';
        end
        b = bar(lags,c,'FaceColor','r','BarWidth',1);
%         b.BarWidth(1);
%         b.color('r')
        if iSelectICs2 == SelectedICs(end)
            if length(SelectedICs)>= 4
                [cRand1,lagsRand1] = xcorr(RandEventsTrain.(['IC' num2str(SelectedICs(counter2))]),EventsTrain.(['IC' num2str(SelectedICs(counter2-1))]),MaxLagSamp);
                lagsRand1 = lagsRand1/(FsLags);
                [cRand2,lagsRand2] = xcorr(RandEventsTrain.(['IC' num2str(SelectedICs(counter2-2))]),RandEventsTrain.(['IC' num2str(SelectedICs(counter2-3))]),MaxLagSamp);
                lagsRand2 = lagsRand2/(FsLags);
                subplot(floor(sqrt(length(SelectedICs)))+1,floor(sqrt(length(SelectedICs)))+1,counter2+2)
                b = bar(lagsRand1,cRand2,'FaceColor',Color,'BarWidth',1);
%                 b.BarWidth=1;
%                 b.color('k')
                title('random crosscorrelogram')
            end
        end
    end
end
