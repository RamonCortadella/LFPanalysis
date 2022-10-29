
% barplot of brain states prevalence
directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';

Rec1States = load(strcat(directory,'MatlabData/Rec1-PerStates.mat'));
Rec1States = Rec1States.PerStates;%variable same as data
Rec2States = load(strcat(directory,'MatlabData/Rec2-PerStates.mat'));
Rec2States = Rec2States.PerStates;
Rec3States = load(strcat(directory,'MatlabData/Rec3-PerStates.mat'));
Rec3States = Rec3States.PerStates;
Rec4States = load(strcat(directory,'MatlabData/Rec4-PerStates.mat'));
Rec4States = Rec4States.PerStates;
Rec5States = load(strcat(directory,'MatlabData/Rec5-PerStates.mat'));
Rec5States = Rec5States.PerStates;
Rec6States = load(strcat(directory,'MatlabData/Rec6-PerStates.mat'));
Rec6States = Rec6States.PerStates;
Rec9States = load(strcat(directory,'MatlabData/Rec9-PerStates.mat'));
Rec9States = Rec9States.PerStates;

structRecs = struct('Rec1States',Rec1States,'Rec2States',Rec2States,'Rec3States',Rec3States,'Rec4States',Rec4States,'Rec5States',Rec5States,'Rec9States',Rec9States);%strcut('fiedl', value...)

fieldRecs = fieldnames(structRecs);%list of fieldnames
numhours = 100;
PerREM = zeros(numhours);%matrix 100x100 with zeros
PerSWS = zeros(numhours);
PerTHE = zeros(numhours);
PerNThe = zeros(numhours);

TimeStateStruct = struct('PerREM',PerREM,'PerSWS',PerSWS,'PerTHE',PerTHE,'PerNThe',PerNThe);%structures tieh the matrix of 0 from the 4 staes
fieldStates = fieldnames(TimeStateStruct);%list statesnames
indHour = 0; 
matrix = zeros(numhours,length(fieldStates));
EndRecStamp = 0';
maxindHour = zeros(100);
for ii = 1:length(fieldRecs) %iterate over Rec's
    %detHour = 1;
    %indHour = indHour+1;
    ii
    PerRecState = structRecs.(fieldRecs{ii});%local variable ~every rec!structRecs.Rec1State=s=trcutRecs.(fieldRecs{1})
    for iii = 1:length(fieldStates) %iterate over states 4 states PerREM,...
        PerBrainState = PerRecState.(fieldStates{iii});%rec1state.perrem %combiation of recs and states 
        iii
        for iv = 1:length(PerBrainState(:,1)) %iterate over events in the PerState files
            detHour = floor(PerBrainState(iv,1)/3600)+1;
            indHour = EndRecStamp+detHour;
            %if PerBrainState(iv,2) > 3600*detHour
            if PerBrainState(iv,2)>3600*detHour
                duration = 3600*detHour-PerBrainState(iv,1);
                matrix(indHour,iii) = matrix(indHour,iii) + duration;
                %TimeStateStruct.(fieldStates{iii})(indHour) = TimeStateStruct.(fieldStates{iii})(indHour) + duration;
                duration = PerBrainState(iv,2)-3600*detHour;
                indHour = indHour+1;
                matrix(indHour,iii) = matrix(indHour,iii) + duration;
                %maxindHour(iv)= indHour;
                %detHour = detHour+1
                %TimeStateStruct.(fieldStates{iii})(indHour) = TimeStateStruct.(fieldStates{iii})(indHour) + duration;
            else
                duration = PerBrainState(iv,2)-PerBrainState(iv,1);
                matrix(indHour,iii) = matrix(indHour,iii) + duration;
                %TimeStateStruct.(fieldStates{iii})(indHour) = TimeStateStruct.(fieldStates{iii})(indHour) + duration
            end 
            
        end
       
    end
     EndRecStamp = indHour;
    
end
B = matrix(~all(matrix == 0,2),:)
figure()
stacked = B;
bar(stacked,'stacked')


legend("PerREM", "PerSWS", 'PerTHE','PerNThe','PerMicroA');

                