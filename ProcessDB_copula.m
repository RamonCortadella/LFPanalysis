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
      
end

% %% compute triggered average
BrainState = 'PerSWS';
InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'','.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';

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
hValues = struct();
counter1=0;
for iSelectICs = SelectedICs
    counter1=counter1+1;
    ICsig = out.(FreqBandLow{1}).ICAsig(iSelectICs,:);
    if counter1~=2
        continue
    end
    figure()
    counter2=0;
    for iSelectICs2 = SelectedICs
        counter2 = counter2 +1;
        ICsig2 = out.(FreqBandLow{1}).ICAsig(iSelectICs2,:);
        
        subplot(floor(sqrt(length(SelectedICs)))+1,floor(sqrt(length(SelectedICs)))+1,counter2)
        h = histogram2(ICsig,ICsig2,linspace(-3,3,21),linspace(-3,3,21),'DisplayStyle','tile');
        colorbar
        hValues.(['IC' num2str(iSelectICs) 'to' num2str(iSelectICs2)]) = h.Values;
%         h=histogram2('XBinEdges',edges(1),'YBinEdges',edges(2),'BinCounts',counts);
    end
    
    figure()
    counter2=0;
    for iSelectICs2 = SelectedICs
        counter2 = counter2 +1;
        ICsig2 = ICsig2(randperm(length(ICsig2)));
        
        subplot(floor(sqrt(length(SelectedICs)))+1,floor(sqrt(length(SelectedICs)))+1,counter2)
        h2 = histogram2(ICsig,ICsig2,linspace(-3,3,21),linspace(-3,3,21),'DisplayStyle','tile');
        colorbar
    end
    figure()
    counter3=0;
    values=struct();
    for iSelectICs2 = SelectedICs
        counter3 = counter3 +1;
        values.(['IC' num2str(iSelectICs) 'to' num2str(iSelectICs2)]) = hValues.(['IC' num2str(iSelectICs) 'to' num2str(iSelectICs2)])-h2.Values;
%         values =abs(values);
%         valuesNeg = -values;
%         values(find(values<=0))=max(max(values))*2;
%         valuesNeg(find(valuesNeg<=0))=max(max(valuesNeg))*2;
        subplot(floor(sqrt(length(SelectedICs)))+1,floor(sqrt(length(SelectedICs)))+1,counter3)
        im = imagesc(values.(['IC' num2str(iSelectICs) 'to' num2str(iSelectICs2)]));
%         h3 = histogram2('XBinEdges',h.XBinEdges,'YBinEdges',h.YBinEdges,'BinCounts',values,'DisplayStyle','tile');%negative saturated (shows positive correlations)
%         h3 = histogram2('XBinEdges',h.XBinEdges,'YBinEdges',h.YBinEdges,'BinCounts',valuesNeg,'DisplayStyle','tile'); %positive saturated (shows negative correlations
        colorbar
        colormap jet
        clim([-0.5*max(max(abs(values.(['IC' num2str(iSelectICs) 'to' num2str(iSelectICs2)])))) 0.5*max(max(abs(values.(['IC' num2str(iSelectICs) 'to' num2str(iSelectICs2)]))))]);
    end
    set(im,'XData',[min(h.XBinEdges) max(h.XBinEdges)])
    set(im,'YData',[min(h.YBinEdges) max(h.YBinEdges)])
    
    
end
