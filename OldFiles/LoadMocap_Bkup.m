
T = readtable(FileName);
% rec = FileName(end-7:end-4);
% T = readtable('../Mocap/O14-SH5/Take 2021-12-09 02.16.21 PM-rec1.csv');

timeA = table2array(T(:,2));
time = str2double(timeA(5:end));

xA = table2array(T(:,3));
x = str2double(xA(5:end))+1;


yA =table2array(T(:,4));
y = str2double(yA(5:end))+1;

zA = table2array(T(:,5));
z = str2double(zA(5:end))+1;

% figure()
% plot3(x(1:10000),y(1:10000),z(1:10000))
%% filter and plot velocity
close('all')
Kernel = 180;
cx = conv(x, ones([Kernel,1]));
dx = (cx(2:end)-cx(1:end-1))./Kernel;


cy = conv(y, ones([Kernel,1]));
dy = (cy(2:end)-cy(1:end-1))./Kernel;

cz = conv(z, ones([Kernel,1]));
dz = (cz(2:end)-cz(1:end-1))./Kernel;

% figure()
% plot(cx,cy,cz)


ModVel = (dx.^2+dy.^2+dz.^2).^(1/2);

figure()
plot(time(2:end), ModVel(Kernel:end)*100)
ylabel('v (cm/s)')
xlabel('t (s)')

[Periods, StateMap, StateTitle] = MotorStates([time(2:end), ModVel(Kernel:end)*100],0.05,10, 180, 30, 60); %min v in cm/s

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

legend('Quite')
plot(time(2:end), ModVel(Kernel:end)*100,'k')
ylabel('v (cm/s)')
xlabel('t (s)')

MotorStateMap(:,1) = [0:1:length(StateMap)-1]./180;
MotorStateMap(:,2) = StateMap;
MotorStateMap(:,3) = ModVel(Kernel:end)*100;

MocapDuration = time(end);

save(strcat(strcat(directory,'MatlabData/'),rec,'-MotorState.mat'),'MotorStateMap');
save(strcat(strcat(directory,'MatlabData/'),rec,'-Periods.mat'),'Periods');
save(strcat(strcat(directory,'MatlabData/'),rec,'-MocapDuration.mat'),'MocapDuration');

%% main movie display plots
video=false;

if video == true

    start = 1;
    stop = 180*100;

    skip = 5; %downsampling of video
    
    figure()
    ax1 = plot3(x(1:10),y(1:10),z(1:10));
    xlim([min(x(start:stop)) max(x(start:stop))])
    ylim([min(y(start:stop)) max(y(start:stop))])
    zlim([min(z(start:stop)) max(z(start:stop))])
    
    counter=0;
    for k=start:skip:stop%size(LfpGeom(:,1),1)
        display(k)
        counter = counter+1;
        plot3(x(k:k+10),y(k:k+10),z(k:k+10))
        xlim([min(x(start:stop)) max(x(start:stop))])
        ylim([min(y(start:stop)) max(y(start:stop))])
        zlim([min(z(start:stop)) max(z(start:stop))])
        M(counter) = getframe(gcf);

    end

    %%
    writerObj = VideoWriter('/storage2/ramon/testVideoMocap.avi');
    writerObj.FrameRate = 180/skip;
    % open the video writer
    open(writerObj);

    % write the frames to the video
    for i=1:floor((stop-start)/skip)%length(M)
        % convert the image to a frame
        frame = M(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end