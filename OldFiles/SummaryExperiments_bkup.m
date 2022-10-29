MetaDataBase = '/storage2/ramon/data/DB/RecordingsDB.xlsx';
InputPath = '../../../data/DB/';
OutputPath = '../../../data/DB/';
%% Functions Settings
%--- save?
SaveResult = false;
%--- smoothening
smooth = false;
iterSmoothening = 5;
%--- bad channels identification
BadChannels = false;
sigmas = 2;
%--- find trigger in signal
FindTrigger = true;

%% define queries
Queries={};
% Queries.depth=1;
% Queries.coverage=1;
Queries.flag=1;
fn = fieldnames(Queries);

%% load MetaDataBase  and apply queries/return selected files
T = readtable(MetaDataBase);
T = rmmissing(T);

for i = 1:length(fn)
    ind = find(table2array(T(:,fn(i)))==Queries.(fn{i}));
    st = table2array(T(ind,'recordingId'));
    if i >= 2
        [st,~] = intersect(st,st0);
        [ind,~] = intersect(ind,ind0);
    end
    st0=st;
    ind0=ind;
end
val = st;

%% apply functions on loaded recordings
for iFile = 1:length(val)
    nCh = T.NumCh(ind(iFile));
    Lfp = LoadBinaryDAT([InputPath, val{iFile}], [0:nCh-1], nCh,1)';
    
    
    nRows = T.nRows(ind(iFile));
    nCols = T.nCols(ind(iFile));
    Chs = [1:nRows*nCols];
    
    Fs = T.Fs(ind(iFile));
%---------------    
    if smooth == true
        LfpCorr = GmSmoothening(Lfp, [0.5, 10], Fs, nRows, nCols, Chs, iterSmoothening);
    end
%---------------   
    if BadChannels == true
        badChs = FindBadCh(Lfp, Fs, nRows, nCols, sigmas);
        if length(find(ismember(T.Properties.VariableNames,'badChannels')==1))>=1
            T.badChannels(ind(iFile))= {mat2str(badChs)};
            writetable(T,strcat(OutputPath,'RecordingsDB-Temp','.xlsx'))
        else
            display('Bad channels variable not present in table')
        end
    end
%--------------- 
    if SaveResult == true
        SaveBinary(strcat(OutputPath,val{iFile}(1:end-4),'-Smooth','.lfp'), LfpCorr)
    end    


%     LfpCorr = ButFilter(LfpCorr,4,[0.5, 30]/(Fs/2),'bandpass');

%     [xlsxFields,indRec,values] = FindBadCh(Lfp, Fs);
%     T.badChannels{ind(iFile)}=mat2str(badChs);%generalize to different
%     outputs from my fun
    

end
%% write a log file with queries and functions used on which date