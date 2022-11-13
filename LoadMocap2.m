function LoadMocap2(FileName,T,IndDB, InputPath, OutputPath,OutputPathTable, FileNameDB, varargin)
% this function loads and processes (low pass filter and formating)
% mocap data and calls "motorstates" function which classifies the motor
% state
[fMode, vThreshold, duration, FallSleepTime, MinSleepTime] = DefaultArgs(varargin,{'compute',0.05, 5, 20, 40});

%body of the function
switch fMode
    case 'compute'
        %if you give output from the function make it structure :
        nCh = T.NumCh(IndDB);
        
        fn = split(FileName,'-');
        rn = split(fn{3},'.');
        MocapFile = strcat(fn{1},'-',fn{2},'-rec',rn{1}(4:end),'.csv');
        MocapT = readtable(strcat(InputPath,MocapFile)); %load csv file

        %load params
        Fs = T.Fs(IndDB);
        FsMocap = T.FsMocap(IndDB);
        
        
        Kernel = FsMocap;
        var = zeros([length(table2array(MocapT(5:end,2))), 4]);
        cvar = zeros([length(table2array(MocapT(5:end,2)))+Kernel-1, 3]); %init variable for filtered coordinates
        dcvar = zeros([length(table2array(MocapT(5:end,2)))+Kernel-2, 3]);%init variable for derivative of filtered coordinates

        %%assign columns and filter with convolution
        for i = 1:4
            if i>=2
                var(:,i)= str2double(table2array(MocapT(5:end,i+1)))+1; %var(1)=time, var(2) = x, var(3)=y, var(4)=z
                
                cvar(:,i-1)=conv(var(:,i),ones([Kernel,1]));
                dcvar(:,i-1) = (cvar(2:end,i-1)-cvar(1:end-1,i-1))./Kernel;
                
            else
                var(:,i)= str2double(table2array(MocapT(5:end,i+1)));
            end
        end

        ModVel = (dcvar(:,1).^2+dcvar(:,2).^2+dcvar(:,3).^2).^(1/2);
        ModVel(isnan(ModVel))=0;
        [Periods, StateMap, StateTitle] = MotorStates([var(2:end,1), ModVel(Kernel:end)*100],vThreshold,duration, FsMocap, FallSleepTime, MinSleepTime); %min v in cm/s

        MotorStateMap(:,1) = [0:1:length(StateMap)-1]./FsMocap;
        MotorStateMap(:,2) = StateMap;
        %Add broken=1 to database if the velocity length and
        Broken=0;
        if length(MotorStateMap(:,2)) ~= length(ModVel(Kernel:end))
            Broken = 1;
            return
        end
        MotorStateMap(:,3) = ModVel(Kernel:end)*100;

        MocapDuration = var(end,1);

        fn = split(FileName,'-');
        rn = split(fn{3},'.');
        
        save(strcat(OutputPath,fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-MotorState.mat'),'MotorStateMap');
        save(strcat(OutputPath,fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-Periods.mat'),'Periods');
        save(strcat(OutputPath,fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-MocapDuration.mat'),'MocapDuration');
%         save(strcat(OutputPath,fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-MocapDuration.mat'),'MocapDuration');
        
%         T.LostMocapSamples(IndDB) = Broken;
%         out.T = T;
%         writetable(T,strcat(OutputPathTable,FileNameDB,'.xlsx'))
        
    case 'display'
        fn = split(FileName,'-');
        rn = split(fn{3},'.');
        
        MotorStateMap = load(strcat(OutputPath,fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-MotorState.mat'));
        Periods = load(strcat(OutputPath,fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-Periods.mat'));
        MotorStateMap = MotorStateMap.MotorStateMap;
        Periods= Periods.Periods;
        
        figure()
        plot(MotorStateMap(:,1), MotorStateMap(:,3))
        ylabel('v (cm/s)')
        xlabel('t (s)')
        
        %plot motor state classification and export figure
        PerRun = Periods(1);
        PerQui = Periods(2);
        PerSle = Periods(3);

        PerRun = PerRun{1};
        PerQui = PerQui{1};
        PerSle = PerSle{1};
    % 
        figure()
        hold on

    %     PerQui= unique(PerQui);
        sQui = size(PerQui);
        if sQui(1)>1 | sQui(2)>1
            for i = 1:length(PerQui)

                display(PerQui)
                display(i)
                n(i) = area([PerQui(i,1) PerQui(i,2)],[1 1]);
                n(i).FaceColor = [0 1 0];
                n(i).LineStyle = 'none';

                n(i).FaceAlpha = 0.2;
            end
        end
    %     PerRun= unique(PerRun);
        sRun = size(PerRun);
        if sRun(1)>1 | sRun(2)>1
            for i = 1:length(PerRun)
                display([PerRun(i,1) PerRun(i,2)])
                g(i) = area([PerRun(i,1) PerRun(i,2)],[1 1]);
                g(i).FaceColor = [1 0 0];
                g(i).LineStyle = 'none';
                g(i).FaceAlpha = 0.2;
            end
        end

    %     PerSle= unique(PerSle);
        sSle = size(PerSle);
        if sSle(1)>1 | sSle(2)>1
            for i = 1:length(PerSle)
                h(i) = area([PerSle(i,1) PerSle(i,2)],[1 1]);
                h(i).FaceColor = [0 0 1];
                h(i).LineStyle = 'none';

                h(i).FaceAlpha = 0.2;
            end
        end
%         savefig(strcat(OutputPath,'figures/',fn{1},'-',fn{2},'-Rec',rn{1}(4:end),'-FigMotorStates.fig'));
        
    case 'groupstats'
        % load and then compute stats and output 
end
