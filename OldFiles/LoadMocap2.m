function LoadMocap(FileName, InputPath, OutputPath)
    T = readtable(strcat(InputPath,FileName)); %load csv file

    close('all')% filter with kernel=number of samples (for convolution) and plot velocity
    Kernel = 180;
    var = zeros([length(T(5:end,2)), 4]);
    cvar = zeros([length(T(5:end,2))+Kernel, 3]); %init variable for filtered coordinates
    dcvar = zeros([length(T(5:end,2))-1, 3]);%init variable for derivative of filtered coordinates
    
    %%assign columns and filter with convolution
    for i = 1:4
        if i>=2
            var(:,i)= str2double(table2array(T(:,i+1)))+1; %var(1)=time, var(2) = x, var(3)=y, var(4)=z
            cvar(:,i)=conv(var(:,i),ones([kernel,1]));
            dcvar(:,i) = (cvar(2:end,i)-cvar(1:end-1,i))./Kernel;
        else
            var(:,i)= str2double(table2array(T(:,i+1)));
        end
    end
    
    ModVel = (dcvar(:,1).^2+dcvar(:,2).^2+dcvar(:,3).^2).^(1/2);

    figure()
    plot(time(2:end), ModVel(Kernel:end)*100)
    ylabel('v (cm/s)')
    xlabel('t (s)')

    [Periods, StateMap, StateTitle] = MotorStates([time(2:end), ModVel(Kernel:end)*100],0.05,10, 180, 30, 60); %min v in cm/s
    
    MotorStateMap(:,1) = [0:1:length(StateMap)-1]./180;
    MotorStateMap(:,2) = StateMap;
    MotorStateMap(:,3) = ModVel(Kernel:end)*100;

    MocapDuration = time(end);
    
    fn = split(FileName,'-');
    rn = split(fn{3},'.');
    save(strcat(fn{1},'-',fn{2},'-',rn{1},'-MotorState.mat'),'MotorStateMap');
    save(strcat(fn{1},'-',fn{2},'-',rn{1},'-Periods.mat'),'Periods');
    save(strcat(fn{1},'-',fn{2},'-',rn{1},'-MocapDuration.mat'),'MocapDuration');
    
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
    for i = 1:length(PerQui)
        display(i)
        n(i) = area([PerQui(i,1) PerQui(i,2)],[1 1]);
        n(i).FaceColor = [0 1 0];
        n(i).LineStyle = 'none';

        n(i).FaceAlpha = 0.2;
    end

    for i = 1:length(PerRun)
        display([PerRun(i,1) PerRun(i,2)])
        g(i) = area([PerRun(i,1) PerRun(i,2)],[1 1]);
        g(i).FaceColor = [1 0 0];
        g(i).LineStyle = 'none';
        g(i).FaceAlpha = 0.2;
    end

    for i = 1:length(PerSle)

        display([PerSle(i,1) PerSle(i,2)])
        h(i) = area([PerSle(i,1) PerSle(i,2)],[1 1]);
        h(i).FaceColor = [0 0 1];
        h(i).LineStyle = 'none';

        h(i).FaceAlpha = 0.2;
    end
    
    savefig(strcat(fn{1},'-',fn{2},'-',rn{1},'-FigMotorStates.mat');
   
end