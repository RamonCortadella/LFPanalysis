function Video_powerbands(ICt,AvIC, AvPowerT,AvPower,outpath, T, indDB, varargin)
% This function generates a video from IC and power bands  
%
% Lfp : data AC channels in [length signal, number of channels], index
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

[SkipSampsFact, SpeedFact, CnstScale] = DefaultArgs(varargin,{1, 0.1, true});


nChRow = T.nRows(indDB);
nChCol = T.nCols(indDB);
Fs = T.Fs(indDB);
FsSpec = 1/(AvPowerT(2)-AvPowerT(1));
%% movie displas
size(ICt')
size(AvIC')
size([ICt',AvIC'])
size(AvPowerT)
AvICint = Interpolate([ICt',AvIC'],AvPowerT,'trim','off');
size(AvICint)

bands = fieldnames(AvPower);
nBands = length(fieldnames(AvPower));


tax = AvICint(:,1);
f=figure(110);clf
f.Position(3:4) = [300,1000];


hold on
ax1 = subplot(nBands+1,1,1);
plot(tax,AvICint(:,2));
display(tax(1))
h1 = Lines(tax(1),[],'k');
%     linkaxes([ax1, ax2],'x')

Min = struct();
% hold off
for iBands = 1:nBands
    Min.(bands{iBands}) = max(max(abs(AvPower.(bands{iBands}))));
end
%%
counter= 0;
for k=1:SkipSampsFact:length(AvICint(:,1))
    counter = counter+1;
    display(k)
    
    for iBands = 1:nBands

%         tax = AvICint(:,1);
     
        ax2 =subplot(nBands+1,1,iBands+1);

        bF = reshape(AvPower.(bands{iBands})(k,:),[nChRow,nChCol]);
  
        imagesc(bF);
        if CnstScale == true
            caxis([-0.7*Min.(bands{iBands}) 0.7*Min.(bands{iBands})]);
        else
            caxis(prctile(bF(:),[1 99])); 
        end
        colorbar
%         hold off

        subplot(nBands+1,1,1);
%         plot([AvICint(k,1) AvICint(k,1)],[-max(AvICint(:,2)) max(AvICint(:,2))],'k')

        set(h1,'XData',[1 1]*AvICint(k,1));%AvICint(k,1)/FsSpec
%         set(gcf, 'PaperSize', [6 2]);
        M(k) = getframe(gcf);
    end

end

%%
writerObj = VideoWriter(outpath);
writerObj.FrameRate = SpeedFact*Fs/SkipSampsFact;
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:SkipSampsFact:length(M)
    % convert the image to a frame
    frame = M(i) ;    
%     display(frame)
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
end