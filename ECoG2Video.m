function ECoG2Video(LfpGeom, LfpGeomDC, outpath, FsAC, DownFact, FminAC, FmaxAC, FminDC, FmaxDC,OrderFiltAC, OrderFiltDC, SlowDown)
% This function generates a video from ECoG data 
%
% outpath : path where the video is saved eg. ='/user/VideoName.avi';
%
% FsAC : sampling frequency of AC channels
% DownFact : downsampling factor between AC and DC channels
%
% Tstart : time of video beginning w.r.t. AC and DC data
% Tstop : time of video end w.r.t. AC and DC data
%
% FminAC : minimum frequency of AC channels video
% FmaxAC : maximum frequency of AC channels video
% FminDC : minimum frequency of DC channels video
% FmaxDC : maximum frequency of DC channels video
% OrderFiltAC : order of bandpass Butter filter applied on AC channels
% OrderFiltDC : order of bandpass Butter filter applied on DC channels 
%
% SlowDown : factor to slow down video
% --------------
% By Ramon Garcia Cortadella 27.03.2022, inspired in code by Anton Sirota
% --------------

% Settings -----

Fs = FsAC;
FsDC = FsAC/DownFact;

Coord(:,1) = reshape(repmat([1:16],16,1),[],1);
Coord(:,2) = repmat([1:16]',16,1);

map = nan(256,1);

%% Filter ephys data

LfpSpindle = ButFilter(LfpGeom,OrderFiltAC,[FminAC FmaxAC]/(Fs)*2,'bandpass');

LfpISA = ButFilter(LfpGeomDC,OrderFiltDC,[FminDC FmaxDC]/(FsDC)*2,'bandpass');

display('data filtered for video generation')

%% main movie display plots
close('all')

start = 1; %start in samples 
stop = size(LfpSpindle,1);   

tax = [1:stop]*1000/Fs;
[bX, bY]= meshgrid([1:16]',[1:16]');
figure(110);clf
ax1 = subplot(4,4,[1:2 5 6]);
imagesc(tax,[],LfpSpindle(1:end,:)'); caxis([-1500 1500]);
%imagescnan(lsegdat',[-1 1]*4000);
tax = [1:size(LfpSpindle,1)]*1000/Fs;
h = Lines(tax(1),[],'k');
ax2 = subplot(4,4,[9 10]);
hold on
% plot(tax,LfpGeom(:,187),'k'); 
plot(tax,LfpSpindle(:,linspace(1,241,16)+linspace(0,15,16)));
% plot(tax,LfpLowGamma(:,187),'r');
%imagescnan(lsegdat',[-1 1]*4000);
h1 = Lines(tax(1),[],'k');
linkaxes([ax1, ax2],'x')
xlim([start stop])

ylim([-1.5*abs(min(min(LfpSpindle(:,linspace(1,241,16)+linspace(0,15,16))))) 1.5*abs(max(max(LfpSpindle(:,linspace(1,241,16)+linspace(0,15,16)))))])

ax3 = subplot(4,4,[13 14]);

taxDC = [1:size(LfpISA,1)]*1000/(FsDC);
plot(taxDC,LfpISA(:,linspace(1,241,16)+linspace(0,15,16)));
% plot(tax,LfpLowGamma(:,187),'r');
%imagescnan(lsegdat',[-1 1]*4000);
h2 = Lines(taxDC(1),[],'k');
xlim([start stop])
ylim([-1.5*abs(min(min(LfpISA(:,linspace(1,241,16)+linspace(0,15,16))))) 1.5*abs(max(max(LfpISA(:,linspace(1,241,16)+linspace(0,15,16)))))])




hold off
%%
counter= 0;
for k=start:stop%size(LfpGeom(:,1),1)
    counter = counter+1;
%  k=1;
 % while 1
%     map = LfpGeom(k,:);
%     
%     F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
%     bF = F(bX,bY);
%     bF = imgaussfilt(bF,1);
%     bF=medfilt2(bF);

%     ax3 = subplot(3,2,2);cla
%     imagesc(bF);

%     pbaspect(ax3,[2,2,1])
%     caxis([-1 1]*prctile(abs(map),99));
%     title(num2str(k));

    % add detected local extrema
%     hold on
%     next = Extr.frame==k & Extr.Sign==-1;
%     plot(Extr.Coord(next,1), Extr.Coord(next,2),'w*','MarkerSize',20);
%     next = Extr.frame==k & Extr.Sign==1;
%     plot(Extr.Coord(next,1), Extr.Coord(next,2),'wo','MarkerSize',20);
%     

    %%
    map = LfpSpindle(k,:);   
    F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
    bF = F(bX,bY);

%     bF = imgaussfilt(bF,1);
    bF = imgaussfilt(bF,0.5);
%     bF=medfilt2(bF);
    ax4 = subplot(4,4,[3:4 7 8]);cla
    imagesc(bF);
    pbaspect(ax4,[2,2,1])
%     caxis([-4000 4000]);
    caxis(prctile(bF(:),[1 99]));
    %%
    kDC = floor(start/(Fs/FsDC))+mod(floor(counter*10/(Fs/FsDC)),floor(stop-start)/(Fs/FsDC))+1;%floor(start/(Fs/FsDC))+mod(k,floor(stop-start)/10)+1;%;%
    if k <= floor(counter*10/(Fs/FsDC))
        break
    end
%     display(floor(start/(Fs/FsDC)))
%     display(mod(floor(counter*10/(Fs/FsDC)),floor(stop-start)/(Fs/FsDC))+1)
%     display(kDC)
    
    kDCh = start+mod(counter*10,floor(stop-start));%mod(k,floor(stop-start)/10)+1;% kDC-floor(start/(Fs/FsDC))+start;
    map = LfpISA(kDC,:);   
    F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
    bF = F(bX,bY);
%     
    bF = imgaussfilt(bF,0.5);
%     bF=medfilt2(bF);
    ax5 = subplot(4,4,[11:12 15 16]);cla
    imagesc(bF);
    pbaspect(ax5,[2,2,1])
%     caxis([-1200 1200]);
    caxis(prctile(bF(:),[1 99]));    



    subplot(4,4,[1:2 5 6]);
    set(h,'XData',[1 1]*(k));
    subplot(4,4,[9 10]);
    set(h1,'XData',[1 1]*(k));
    subplot(4,4,[13 14]);
    set(h2,'XData',[1 1]*(kDCh));


    M(k-start+1) = getframe(gcf);

end

%%
writerObj = VideoWriter(outpath);
writerObj.FrameRate = Fs/(SlowDown);
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(M)
    % convert the image to a frame
    frame = M(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
end