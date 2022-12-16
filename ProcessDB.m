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
InputPath = strcat('../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'); %the DB search must point to the files used for ICA computation to access associated metadata
OutputPath = '../../../data/DB/Files/B13289O14-DH3/GroupsAnalysis/';
out = SingularSpectrum(InputPath,OutputPath,T,ind(iFile),extension);

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
