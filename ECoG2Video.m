function ECoG2Video(LfpGeom, LfpGeomDC, outpath, FsAC, DownFact,SpeedDCFact,SpeedACFact, FminAC, FmaxAC, FminDC, FmaxDC,OrderFiltAC, OrderFiltDC, SlowDown, CnstScale, GaussianFact ,nCh,nChCol,nChRow,double)
% This function generates a video from ECoG data 
%
% LfpGeom : data AC channels in [length signal, number of channels], index
% 1 is anterior left channel, index (last) is posterior right.
%
% outpath : path where the video is saved eg. ='/user/VideoName.avi';
%
% FsAC : sampling frequency of AC channels
% DownFact : downsampling factor between AC and DC channels
%
% SpeedACFact : Takes every 'SpeedACFact'th frame to create the top plot
% SpeedDCFact : Takes every 'SpeedDCFact'th frame to create the bottom plot
%               it loops until video in top plot ends
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
% SlowDown : frame rate is Fs/SlowDown
% --------------
% By Ramon Garcia Cortadella 27.03.2022, inspired in code by Anton Sirota
% --------------

% Settings -----

if nargin < 15
    GaussianFact = false;
end
    
Fs = FsAC;
FsDC = FsAC/DownFact;

Coord(:,1) = reshape(repmat([1:nChRow],nChCol,1),[],1);
Coord(:,2) = repmat([1:nChCol]',nChRow,1);

map = nan(nCh,1);

%% Filter ephys data
display(max(max(LfpGeom)))
LfpSpindle = ButFilter(LfpGeom,OrderFiltAC,[FminAC FmaxAC]/(Fs)*2,'bandpass');

display(max(max(LfpSpindle)))

display([FminDC FmaxDC])
LfpISA = ButFilter(LfpGeomDC,OrderFiltDC,[FminDC FmaxDC]/(FsDC)*2,'bandpass');
display(max(max(LfpISA)))

display('data filtered for video generation')

%% main movie display plots
close('all')

start = 1; %start in samples 
stop = size(LfpSpindle,1);   
display('stop')
display(size(LfpSpindle))

tax = [1:stop]*1000/Fs;
[bX, bY]= meshgrid([1:nChRow]',[1:nChCol]');
figure(110);clf
ax1 = subplot(4,4,[1:2 5 6]);
imagesc(tax,[],LfpSpindle(1:end,:)'); caxis([-1000 1000]);
% caxis([-0.6*abs(min(min(LfpSpindle(1:end,:)))) 0.6*abs(max(max(LfpSpindle(1:end,:))))])
%imagescnan(lsegdat',[-1 1]*4000);
tax = [1:size(LfpSpindle,1)]*1000/Fs;
h = Lines(tax(1),[],'k');
if double == true
    ax2 = subplot(4,4,[9 10 13 14]);
else
    ax2 = subplot(4,4,[9 10]);
end

hold on
% plot(tax,LfpGeom(:,187),'k'); 
plot(tax,LfpSpindle(:,linspace(1,241,16)+linspace(0,15,16)));
% plot(tax,LfpLowGamma(:,187),'r');
%imagescnan(lsegdat',[-1 1]*4000);
h1 = Lines(tax(1),[],'k');
linkaxes([ax1, ax2],'x')
xlim([start stop]*1000/Fs)

ylim([-1.5*abs(min(min(LfpSpindle(:,linspace(1,241,16)+linspace(0,15,16))))) 1.5*abs(max(max(LfpSpindle(:,linspace(1,241,16)+linspace(0,15,16)))))])

if double == false
    ax3 = subplot(4,4,[13 14]);
    taxDC = [1:size(LfpISA,1)]*1000/(FsDC);
    plot(taxDC,LfpISA(:,linspace(1,241,16)+linspace(0,15,16)));
    % plot(tax,LfpLowGamma(:,187),'r');
    %imagescnan(lsegdat',[-1 1]*4000);
    h2 = Lines(taxDC(1),[],'k');
    xlim([start stop]*1000/Fs)
    ylim([-1.5*abs(min(min(LfpISA(:,linspace(1,241,16)+linspace(0,15,16))))) 1.5*abs(max(max(LfpISA(:,linspace(1,241,16)+linspace(0,15,16)))))])

end




hold off
%%
counter= 0;
for k=start:SpeedACFact:stop%size(LfpGeom(:,1),1)
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
    display(size(bF))
  
    if GaussianFact
        bF = imgaussfilt(bF,GaussianFact);
    end
    if double == true
        ax4 = subplot(4,4,[3:4 7 8 11 12 15 16]);cla
        imagesc(bF);
        pbaspect(ax4,[2,4,1])  
    else
        ax4 = subplot(4,4,[3:4 7 8]);cla
        imagesc(bF);
        pbaspect(ax4,[2,2,1])   
    end

    if CnstScale == true
%         caxis([-0.3*abs(min(min(LfpSpindle(1:end,:)))) 0.3*abs(min(min(LfpSpindle(1:end,:))))])
        caxis([-100 100]);
    else
        caxis(prctile(bF(:),[1 99])); 
    end
    colorbar
    %%
    kDC = 1+mod(floor(counter*SpeedDCFact/(Fs/FsDC)),floor(stop-start)/(Fs/FsDC))+1;%floor(start/(Fs/FsDC))+mod(k,floor(stop-start)/10)+1;%;%

%     display(floor(start/(Fs/FsDC)))
%     display(mod(floor(counter*10/(Fs/FsDC)),floor(stop-start)/(Fs/FsDC))+1)
%     display(kDC)
    
%     kDCh = start+mod(counter*SpeedDCFact,floor(stop-start));%mod(k,floor(stop-start)/10)+1;% kDC-floor(start/(Fs/FsDC))+start;
    if floor(counter*SpeedDCFact/(Fs/FsDC)) >= floor(stop-start)/(Fs/FsDC)
        break
    end
    map = LfpISA(kDC,:);   
    F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map');
    bF = F(bX,bY);
    if GaussianFact
        bF = imgaussfilt(bF,GaussianFact);
    end
    if double == false
           
        ax5 = subplot(4,4,[11:12 15 16]);cla
        imagesc(bF);
        pbaspect(ax5,[2,2,1])
        if CnstScale == true
%             caxis([-0.2*abs(min(min(LfpISA(1:end,:)))) 0.2*abs(max(max(LfpISA(1:end,:))))])
            caxis([-100 100]);
        else
            caxis(prctile(bF(:),[1 99]));    
        end
        colorbar
    end
    
    
    subplot(4,4,[1:2 5 6]);
    set(h,'XData',[1 1]*(k*1000/Fs));
    
    if double == true
        subplot(4,4,[9 10 13 14]);
    else
        subplot(4,4,[9 10]);
    end
    
    set(h1,'XData',[1 1]*(k*1000/Fs));

    if double == false
        subplot(4,4,[13 14]);
        set(h2,'XData',[1 1]*(kDC*1000/FsDC));
    end
    


    M(k-start+1) = getframe(gcf);

end

%%
writerObj = VideoWriter(outpath);
writerObj.FrameRate = Fs/(SlowDown);
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:SpeedACFact:length(M)
    % convert the image to a frame
    frame = M(i) ;    
    display(frame)
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
end