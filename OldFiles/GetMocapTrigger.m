close('all')
%% Define settings

directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/Interpolated_NotFilteredLFP_HP_1HzAC_0p005HzDC/NotClippedMap/B13289O14-DH1-Rec2AC.xml'));
FileName = strcat(directory,'DatData/Interpolated_NotFilteredLFP_HP_1HzAC_0p005HzDC/NotClippedMap/B13289O14-DH1-Rec2DC.dat');

Fs = 651.04166667;
DownSamp = 1;

Fs = Fs/DownSamp;

Tstab = 60; %Time removed from the beggining to remove artifactual peaks not related to trigger
TmaxInit = 200; % Time from the beggining taken to look for start trigger
TmaxEnd = 200; % Time from the end taken to look for stop trigger
%% get ephys data

% FileName = [d.folder ,'/',d.name];
LfpDC = LoadBinaryDAT(FileName, [0:255], Par.nChannels,1)';
SplitName = split(FileName,'-');

%% remap geometrically
LfpGeomDC = zeros([size(LfpDC)]);
for i = 1:length(Par.AnatGrps)
    for ii = 1:length(Par.AnatGrps(i).Channels)
        LfpGeomDC(:,ii+(i-1)*16) = LfpDC(:,Par.AnatGrps(i).Channels(ii)+1);
    end
end

%% plot trigger Mocap
figure()
hold on


time = [1:1:length(LfpGeomDC(:,1))]/Fs;
plot(time(Tstab*floor(Fs):TmaxInit*floor(Fs)), LfpGeomDC(Tstab*floor(Fs):TmaxInit*floor(Fs),81))
plot(time(end-TmaxEnd*floor(Fs):end), LfpGeomDC(end-TmaxEnd*floor(Fs):end,81))
hold off

SigStart = LfpGeomDC(Tstab*floor(Fs):TmaxInit*floor(Fs),81);


indexTrigger = find(abs(SigStart)==max(abs(SigStart)));
indexTrigger = indexTrigger(1);

MocapStart = time(Tstab*floor(Fs)+indexTrigger);

display(LfpGeomDC(Tstab*floor(Fs)+indexTrigger-1,81))

indexTrigger = find(abs(LfpGeomDC(end-TmaxEnd*floor(Fs):end,81))==max(abs(LfpGeomDC(end-TmaxEnd*floor(Fs):end,81))));
indexTrigger = indexTrigger(1);

MocapEnd = time(indexTrigger-1+length(LfpGeomDC(:,81))-TmaxEnd*floor(Fs));

display(LfpGeomDC(indexTrigger-1+length(LfpGeomDC(:,81))-TmaxEnd*floor(Fs),81))

%% save Mocap trigger
MocapTrigger = [MocapStart MocapEnd];
save(strcat(strcat(directory,'/MatlabData/'),SplitName{3}(1:end-6),'-MocapTrigger.mat'),'MocapTrigger')

