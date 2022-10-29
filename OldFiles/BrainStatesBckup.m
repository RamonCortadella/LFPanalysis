close('all')
clearvars

directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';

rec = 'Rec9';
%combine Motor states with ephys to perform brain state classification

SpecStates = load(strcat(strcat(directory,'MatlabData/'),rec,'-SpecStates.mat'));
SpecStates = SpecStates.SpecStates;
DCspec = SpecStates.WhitenedSpecStates; %this is not DC channel but it is lower time resolution spectrogram compared to the SpecHVS

SpecHVS = load(strcat(strcat(directory,'MatlabData/'),rec,'-SpecHVS.mat'));
SpecHVS = SpecHVS.SpecHVS;
ACspec = SpecHVS.WhitenedSpecHVS;

MocapTrigger = load(strcat(strcat(directory,'MatlabData/'),rec,'-MocapTrigger.mat'));

TriggerStart = MocapTrigger.MocapTrigger(1);
TriggerStop = MocapTrigger.MocapTrigger(2);
TimeTolerance = 1;

LFP_traces = load(strcat(strcat(directory,'MatlabData/'),rec,'-Traces.mat'));
LfpISA = LFP_traces.LFP_traces.LfpDC_Chs;


FsMocap = 180;
FsLFP = 651.04166667;
FsLFPDC = FsLFP/100;

MotorState.MotorStateMap(:,1) = MotorState.MotorStateMap(:,1)+TriggerStart;
MinThetaDeltaRatio = 0.8;
MinTimeRatio = 20;
MinSpindleDeltaRatio = 3.8;
MinDurationHVS_sec = 1;
[StatePeriods, StateMap, StateTitle, bands] = BrainStatesER_R(DCspec.ws',DCspec.wt',DCspec.wf',MotorState.MotorStateMap,[3 5],[6 9],MinThetaDeltaRatio, MinTimeRatio);
[HVSPeriods, HVSMap, StateTitle, bandsHVS] = HVSclass_R(ACspec.ws',ACspec.wt',ACspec.wf',[0 4], [6 30], MinSpindleDeltaRatio, MinDurationHVS_sec);

%% plot motor states

Periods = Periods.Periods;
PerRun = Periods(1);
PerQui = Periods(2);
PerSle = Periods(3);

PerRun = PerRun{1};
PerQui = PerQui{1};
PerSle = PerSle{1};
% 

% merge sleep states separated by short RUN/quite states (classify them as
% microarousals

counter = 1;
lenPerSle = length(PerSle);
PerMicroA = [];

for i= 1:lenPerSle-1
    display(lenPerSle)
    display(i)
    if PerSle(i+1,1)-PerSle(i,2) <= 10
        PerMicroA(counter,1) = PerSle(i,2);
        PerMicroA(counter,2) = PerSle(i+1,1);
        
        PerSle(i,2)=PerSle(i+1,2);
        PerSle(i+1,:)=[];
        counter = counter+1;
        
        
        display(length(PerSle))
        while PerSle(i+1,1)-PerSle(i,2) <= 10

            PerMicroA(counter,1) = PerSle(i,2);
            PerMicroA(counter,2) = PerSle(i+1,1);
            counter = counter+1; 
            if i+2 <= length(PerSle)
                PerSle(i,2)=PerSle(i+1,2);
                PerSle(i+1,:)=[];
            else
                break
            end
        
            
        end
    end
    display(lenPerSle-1-counter)
    if i >= lenPerSle-1-counter
        display('this is the expected end')
        display(i)
        break
    end
    
end

figure()
ax1 = subplot(7,1,1);
hold on
for i = 1:length(PerQui(:,1))
    display(i)
    a(i) = area([PerQui(i,1)+TriggerStart PerQui(i,2)+TriggerStart],[1 1]);
    a(i).FaceColor = [0 1 0];
    a(i).LineStyle = 'none';
    
    a(i).FaceAlpha = 0.2;
end

for i = 1:length(PerRun(:,1))
    display([PerRun(i,1) PerRun(i,2)])
    b(i) = area([PerRun(i,1)+TriggerStart PerRun(i,2)+TriggerStart],[1 1]);
    b(i).FaceColor = [1 0 0];
    b(i).LineStyle = 'none';
    b(i).FaceAlpha = 0.2;
end



for i = 1:length(PerSle(:,1))
    
    display([PerSle(i,1) PerSle(i,2)])
    c(i) = area([PerSle(i,1)+TriggerStart PerSle(i,2)+TriggerStart],[1 1]);
    c(i).FaceColor = [0 0 1];
    c(i).LineStyle = 'none';
    
    c(i).FaceAlpha = 0.2;
end




% legend('Quite')
plot(MotorState.MotorStateMap(:,1), MotorState.MotorStateMap(:,3),'k')
ylabel('v (cm/s)')
xlabel('t (s)')

%% plot brain states
PerTHE = StatePeriods(1);
PerSWS = StatePeriods(2);
PerREM = StatePeriods(3);
PerNThe = StatePeriods(4);
PerHVS = HVSPeriods{1};

PerTHE = PerTHE{1};
PerSWS = PerSWS{1};
PerREM = PerREM{1};
PerNThe = PerNThe{1};
% 
% % merge Theta states separated by less than 20s
counter2 = 1;
lenPerTHE = length(PerTHE(:,1));
for i= 1:lenPerTHE-1

%     if PerTHE(i+1,1)-PerTHE(i,2) <= 20
% 
%         PerTHE(i,2)=PerTHE(i+1,2);
%         PerTHE(i+1,:)=[];
%         counter2 = counter2+1;
%         


    while PerTHE(i+1,1)-PerTHE(i,2) <= 20

        counter2 = counter2+1; 

        if i+2 <= length(PerTHE)
            PerTHE(i,2)=PerTHE(i+1,2);
            PerTHE(i+1,:)=[];
        else
            break
        end
    end
    % remove NThe periods between the merged THE
%         index = find((PerNThe(:,1)>=PerREM(i,1) & PerNThe(:,1)<=PerREM(i,2)));
%         display('these are the indices')
%         display(indexNThe)
%         PerSWS(indexSWS,:) = [];
%         
%        
    indexNThe = find((PerNThe(:,1)>=PerTHE(i,1) & PerNThe(:,1)<=PerTHE(i,2)));
    PerNThe(indexNThe,:) = [];
    display(indexNThe)
    %\

%     display(lenPerTHE-1-counter2)
    if i >= lenPerTHE-1-counter2
        break
    end
        
    
end

% merge REM states separated by less than 20s
counter2 = 1;
lenPerREM = length(PerREM(:,1));
for i= 1:lenPerREM-1

    if PerREM(i+1,1)-PerREM(i,2) <= 20

        PerREM(i,2)=PerREM(i+1,2);
        PerREM(i+1,:)=[];
        counter2 = counter2+1;
        

        
        while PerREM(i+1,1)-PerREM(i,2) <= 20
                        
            counter2 = counter2+1; 
            
            PerREM(i,2)=PerREM(i+1,2);
            PerREM(i+1,:)=[];
           
        end
        % remove REM and NThe periods between the merged SWS
        indexSWS = find((PerSWS(:,1)>=PerREM(i,1) & PerSWS(:,1)<=PerREM(i,2)));
        display('these are the indices')
        display(indexSWS)
        PerSWS(indexSWS,:) = [];
        
       
        indexNThe = find((PerNThe(:,1)>=PerREM(i,1) & PerNThe(:,1)<=PerREM(i,2)));
        PerNThe(indexNThe,:) = [];
        display(indexNThe)
        %\
    end
    
%     display(lenPerREM-1-counter2)
    if i >= lenPerREM-1-counter2
        break
    end
    
end
% merge SWS states separated by less than 20s
counter3 = 1;
lenPerSWS = length(PerSWS(:,1));
for i= 1:lenPerSWS-1

    if PerSWS(i+1,1)-PerSWS(i,2) <= 20

        PerSWS(i,2)=PerSWS(i+1,2);
        PerSWS(i+1,:)=[];
        counter3 = counter3+1;
        


        while PerSWS(i+1,1)-PerSWS(i,2) <= 20
                        
            counter3 = counter3+1; 
            
            PerSWS(i,2)=PerSWS(i+1,2);
            PerSWS(i+1,:)=[];
           
        end
        % remove REM and NThe periods between the merged SWS
        if PerREM
            indexREM = find((PerREM(:,1)>=PerSWS(i,1) & PerREM(:,1)<=PerSWS(i,2)));
            display('these are the indices')
            display(indexREM)
            PerREM(indexREM,:) = [];
        end
        
        if PerNThe
            indexNThe = find((PerNThe(:,1)>=PerSWS(i,1) & PerNThe(:,1)<=PerSWS(i,2)));
            PerNThe(indexNThe,:) = [];
            display(indexNThe)
        end
        %\

    end
    
%     display(lenPerSWS-1-counter3)
    if i >= lenPerSWS-1-counter3
        break
    end
    
end

% 
% % merge NThe states separated by less than 20s
counter2 = 1;
lenPerNThe = length(PerNThe(:,1));
for i= 1:lenPerNThe-1
%     display(i)
%     display(length(PerNThe(:,1)))
%     display('**')
    if PerNThe(i+1,1)-PerNThe(i,2) <= 20

        PerNThe(i,2)=PerNThe(i+1,2);
        PerNThe(i+1,:)=[];
        counter2 = counter2+1;
        
        if i >= lenPerNThe-1-counter2
            break
        end
        
        while PerNThe(i+1,1)-PerNThe(i,2) <= 20
            counter2 = counter2+1; 
        
            PerNThe(i,2)=PerNThe(i+1,2);
            PerNThe(i+1,:)=[];
            
        end

        %\
    
        % remove NThe periods between the merged THE
        indexSWS = find((PerSWS(:,1)>=PerNThe(i,1) & PerSWS(:,1)<=PerNThe(i,2)));
        display('these are the indices')
        display(indexSWS)
        PerSWS(indexSWS,:) = [];
    %         
    %        
        indexTHE = find((PerTHE(:,1)>=PerNThe(i,1) & PerTHE(:,1)<=PerNThe(i,2)));
        PerTHE(indexTHE,:) = [];
        display(indexTHE)
        
        
    end
    if i >= lenPerNThe-1-counter2
        break
    end
%     display(lenPerTHE-1-counter2)

    
end

ax2 = subplot(8,1,2);
hold on
for i = 1:length(PerSWS(:,1))
%     display(i)
    d(i) = area([PerSWS(i,1) PerSWS(i,2)],[1 1]);
    d(i).FaceColor = [0 0 1];
    d(i).LineStyle = 'none';
    
    d(i).FaceAlpha = 0.2;
end

for i = 1:length(PerTHE(:,1))
%     display([PerTHE(i,1) PerTHE(i,2)])
    e(i) = area([PerTHE(i,1) PerTHE(i,2)],[1 1]);
    e(i).FaceColor = [1 0 0];
    e(i).LineStyle = 'none';
    e(i).FaceAlpha = 0.2;
end

for i = 1:length(PerREM(:,1))
    
%     display([PerREM(i,1) PerREM(i,2)])
    f(i) = area([PerREM(i,1) PerREM(i,2)],[1 1]);
    f(i).FaceColor = [0 1 0];
    f(i).LineStyle = 'none';
    
    f(i).FaceAlpha = 0.2;
end

for i = 1:length(PerNThe(:,1))
    
%     display([PerREM(i,1) PerREM(i,2)])
    g(i) = area([PerNThe(i,1) PerNThe(i,2)],[1 1]);
    g(i).FaceColor = [0 0 0];
    g(i).LineStyle = 'none';
    
    g(i).FaceAlpha = 0.2;
end

if length(PerHVS) >= 1
    for i = 1:length(PerHVS(:,1))

        display([PerHVS(i,1) PerHVS(i,2)])
        y(i) = area([PerHVS(i,1) PerHVS(i,2)],[1 1]);
        y(i).FaceColor = [0 0 0];
        y(i).LineStyle = 'none';

        y(i).FaceAlpha = 1;
    end
end


if length(PerMicroA) >= 1
    for i = 1:length(PerMicroA(:,1))

        display([PerMicroA(i,1) PerMicroA(i,2)])
        j(i) = area([PerMicroA(i,1)+TriggerStart PerMicroA(i,2)+TriggerStart],[1 1]);
        j(i).FaceColor = [0 0 0];
        j(i).LineStyle = 'none';

        j(i).FaceAlpha = 1;
    end
end

%%
FinalTime = max(max([PerREM',PerNThe',PerSWS',PerMicroA',PerTHE']));

if (TriggerStop - TimeTolerance <= FinalTime) &  ( FinalTime <= TriggerStop + TimeTolerance)
    display(strcat('Properly synch Mocap with ',num2str(TimeTolerance),'s time tolerance'))
else
    display('NOT SYNCH WITH MOCAP')
end

    
    
%% 

ax3 = subplot(7,1,3);
h= pcolor(DCspec.wt,DCspec.wf,log10(abs(DCspec.ws')));
ylim([1 50])
caxis([min(min(log10(abs(DCspec.ws)')))+2.5 max(max(log10(abs(DCspec.ws'))))-0.5])
set(gca,'YScale','log')
set(gca,'YDir','normal')
colormap('jet')
set(h, 'EdgeColor', 'none')

% ax4 = subplot(4,1,4);
% plot((Sinit+(1:1:length(LfpGeomDC(Sinit:Sfin,1))))/Fs, mean(LfpGeomDC(Sinit:Sfin,[ChInd-1-16:ChInd+1-16, ChInd-1:ChInd+1, ChInd-1+16:ChInd+1+16]),2))
ax4 = subplot(7,1,4);

%%%%%%%%%%%%%%Ratio Theta delta
ratio = bands.ratio;
map = ratio > MinThetaDeltaRatio;
MinTimeRatio = MinTimeRatio;
if MinTimeRatio >= DCspec.wt(2)-DCspec.wt(1)
    
    SampsMapFilter = floor(MinTimeRatio/(DCspec.wt(2)-DCspec.wt(1)));
    map2(1:SampsMapFilter)= map(1:SampsMapFilter);


    for i = SampsMapFilter:length(map)-SampsMapFilter
        if (map(i-SampsMapFilter+1) == 1) && (map(i+SampsMapFilter)==1)
            map2(i-SampsMapFilter+1:i+SampsMapFilter) = 1;
        end
    end
    for i = SampsMapFilter:length(map)-SampsMapFilter
        if (map(i-SampsMapFilter+1) == 0) && (map(i+SampsMapFilter)==0)
            map2(i-SampsMapFilter+1:i+SampsMapFilter) = 0;
        end
    end      

end
%%%%%%%%%%%%%%%%%Ratio spindle delta
ratioHVS = bandsHVS.ratio;

mapmapHVS = ratioHVS > MinSpindleDeltaRatio;

if MinDurationHVS_sec >= ACspec.wt(2)-ACspec.wt(1)
    
    SampsMapFilter = floor(MinDurationHVS_sec/(ACspec.wt(2)-ACspec.wt(1)));
    mapmapHVS2(1:SampsMapFilter)= mapmapHVS(1:SampsMapFilter);


    for i = SampsMapFilter:length(mapmapHVS)-SampsMapFilter
        if (mapmapHVS(i-SampsMapFilter+1) == 1) && (mapmapHVS(i+SampsMapFilter)==1)
            mapmapHVS2(i-SampsMapFilter+1:i+SampsMapFilter) = 1;
        end
    end
    for i = SampsMapFilter:length(mapmapHVS)-SampsMapFilter
        if (mapmapHVS(i-SampsMapFilter+1) == 0) && (mapmapHVS(i+SampsMapFilter)==0)
            mapmapHVS2(i-SampsMapFilter+1:i+SampsMapFilter) = 0;
        end
    end      
    mapHVS = mapmapHVS2';
end
%%%%%%%%%%%%%%%%


time = DCspec.wt;

plot(time, map2,'k')
ylim([-0.5 1.5])
ylabel('v (cm/s)')
xlabel('t (s)')


ax5 = subplot(7,1,5);
plot(time,ratio);
ylim([0 2])

ax6 = subplot(7,1,6);
plot((1:1:length(LfpISA))/FsLFPDC, LfpISA)
ylim([-3e3 3e3])

% ax7 = subplot(8,1,7);
% h= pcolor(ACspec.wt,ACspec.wf,log10(abs(ACspec.ws')));
% caxis([min(min(log10(abs(ACspec.ws)')))+0.5 max(max(log10(abs(ACspec.ws'))))-0.5])
% set(gca,'YScale','log')
% set(gca,'YDir','normal')
% colormap('jet')
% set(h, 'EdgeColor', 'none')

ax7 = subplot(7,1,7);

time = ACspec.wt;
plot(time,ratioHVS);
ylim([0 2])

linkaxes([ax1,ax2,ax3,ax4, ax5, ax6, ax7],'x')

%% plot all HVS events
MinSpindleDeltaRatio = 4.4;%MinSpindleDeltaRatio-0.2; (Rec1 3.9, rec2,3,4,5,6, 4.2, rec 7,9, 4.4
MinDurationHVS_sec = 1;
[HVSPeriods, HVSMap, StateTitle, bandsHVS] = HVSclass_R(ACspec.ws',ACspec.wt',ACspec.wf',[0 4], [6 30], MinSpindleDeltaRatio, MinDurationHVS_sec);
PerHVS = HVSPeriods{1};
side = floor(sqrt(length(PerHVS(:,1))))+1;
figure()

for i = 1:length(PerHVS(:,1))
    ax = subplot(side,side,i); 
    tinit = PerHVS(i,1);
    tfin = PerHVS(i,2);
    index = find((ACspec.wt >= tinit-5) & (ACspec.wt <= tfin+5)) ;
    h= pcolor(ACspec.wt(index),ACspec.wf,log10(abs(ACspec.ws(index,:)')));
    ylim([0.5 150])
    caxis([min(min(log10(abs(ACspec.ws(index,:))')))+0.5 max(max(log10(abs(ACspec.ws(index,:)'))))-0.5])
    set(gca,'YScale','log')
    set(gca,'YDir','normal')
    colormap('parula')
    set(h, 'EdgeColor', 'none')

end

PerStates = struct('PerREM',PerREM,'PerSWS',PerSWS,'PerTHE',PerTHE,'PerNThe',PerNThe,'PerMicroA',PerMicroA, 'PerHVS', PerHVS);

save(strcat(strcat(directory,'MatlabData/'),rec,'-PerStates.mat'),'PerStates')
