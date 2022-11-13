close('all')
%% Define Queries
Queries={};
% Queries.depth=1;
Queries.coverage=1;
% Queries.flag=1;
Queries.CoupledDC=1;
Queries.mocap=1;
Queries.SingleShank=0;
Queries.LostMocapSamples=0;
Queries.depth=0;
Queries.recordingSystem = 0;
Queries.DeviceName = 'B13289O14-DH3';%{-1,'B13289O14-DH1'}; % if cell array with -1 in first field, the oposite of the query is used

%% Search in DB
FileNameDB = 'RecordingsDB';
[recs, val, ind, T] = SearchDB(['/storage2/ramon/data/DB/' FileNameDB, '.xlsx'], Queries);

%% Apply processing function on retrieved filenames
for iFile = 1:length(recs)
    if iFile~=1
        continue
    end
    
    display(iFile,'Processing File')

    InputPath = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/' T.RecordingId{ind(iFile)}];
    InputPathMotor = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/Mocap/'];
    
    OutputPathTable = '../../../data/DB/';
    OutPathStates = ['../../../data/DB/Files/' T.DeviceName{ind(iFile)} '/' T.RecordingId{ind(iFile)} '/States/'];
    
    FileName = [InputPath,'/',recs{iFile}];

%     factor = 2;
%     out = FindBadCh2(FileName,InputPath,T,ind(iFile), OutputPathTable,FileNameDB, 'compute', factor);
%     T = out.T;
% 
%     LoadMocap2(recs{iFile},T,ind(iFile),InputPathMotor,OutPathStates,OutputPathTable, FileNameDB,'compute');
%
%     GmSmoothening2(FileName,T,ind(iFile),InputPath);
%
%     SpectrogramSave2(strcat(FileName(1:end-4),'.lfp'),OutPathStates,T,ind(iFile),'compute');
% 
%     out = BrainStates2(OutPathStates,OutputPathTable,FileNameDB,T,ind(iFile),ind,'compute');
%     T = out.T;        
%      
%         LfpInterp = InterpolateProbes2(FileName, T,ind(iFile),InputPath);
%
%     LfpReduced = PCA_ICA2(strcat(FileName(1:end-4),'.interp.lfp'),T,ind(iFile),InputPath,'compute',50);
%     
%     extension = '.lfp';
%     ICA2(strcat(FileName(1:end-4),extension),T,ind(iFile),extension,OutPathStates,val,'compute',30,'PerSWS');
    extension = '.ICA.lfp';
    ICA2(strcat(FileName(1:end-4),extension),T,ind(iFile),extension,OutPathStates,val,'display',30,'PerSWS');
    
%     ICA2(strcat(FileName(1:end-4),extension),T,ind(iFile),extension,InputPath,OutPathStates,val,'display',27,'PerSWS');
%     
%     detectDelta(strcat(FileName(1:end-4),'.interp.PCA.lfp'),T,ind(iFile),[OutPathStates T.RecordingId{ind(iFile)}],'display')
%
%     writetable(T,strcat(OutputPathTable,FileNameDB,'.xlsx'))




%     Video2(FileName,[OutPathStates 'Video.avi'], T,ind(iFile))

end
%% Group statistics on pre-processed selected files
 
% BrainStates2(OutPathStates,OutputPathTable,FileNameDB,T,ind(iFile),ind,'groupstats');
% extension ='.PCA.lfp'; %e.g. "smooth.interp.PCA.lfp"
% ICA2(strcat(FileName(1:end-4),'.PCA.lfp'),T,ind(iFile),InputPath,OutPathStates,val,extension,'groupstats');
%   