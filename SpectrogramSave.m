function SpectrogramSave(Lfp,Fs,ChIndHVS,ChIndTheta,OutPathSpectrogram,FileName,nRows)

    % Parameters spectrogram
    nFFTtheta = 3072;
    nFFTHVS = 1024;
%     display(Lfp,'lfp')
    
    %% compute spectrogram whitened and mean-of-sveral-channels for AC theta and HVS channels
    
    [wx,~] = WhitenSignal(Lfp,[],[],[],1);
%     for i=1:length(wx(1,:))
%         wx(isnan(wx(:,i)),i)=0;
%     end
    ChIndThetaList= [ChIndTheta-1-nRows:ChIndTheta+1-nRows, ChIndTheta-1:ChIndTheta+1, ChIndTheta-1+nRows:ChIndTheta+1+nRows];
    ChIndHVSList = [ChIndHVS-1-nRows:ChIndHVS+1-nRows, ChIndHVS-1:ChIndHVS+1, ChIndHVS-1+nRows:ChIndHVS+1+nRows];

%     display(Lfp(:,ChIndThetaList))
    
    [ws,wf,wt] = mtcsglong(mean(wx(:,ChIndHVSList),2),nFFTHVS,Fs);%,1024);
    [s,f,t] = mtcsglong(mean(Lfp(:,ChIndHVSList),2),nFFTHVS,Fs);%,1024);
    
    [wst,wft,wtt] = mtcsglong(mean(wx(:,ChIndThetaList),2),nFFTtheta,Fs);%,1024);
    [st,ft,tt] = mtcsglong(mean(Lfp(:,ChIndThetaList),2),nFFTtheta,Fs);%,1024);

    %% save spectrograms 
    fn = split(FileName,'-');
    rn = split(fn{3},'.');
    
    WhitenedSpecHVS = struct('ws',ws,'wf',wf,'wt',wt);
    RawSpecHVS = struct('s',s,'f',f,'t',t);
    SpecHVS = struct('WhitenedSpecHVS',WhitenedSpecHVS,'RawSpecHVS',RawSpecHVS);
    save(strcat(OutPathSpectrogram,fn{1},'-',fn{2},'-',rn{1},'-SpecHVS.mat'),'SpecHVS')


    WhitenedSpecStates = struct('ws',wst,'wf',wft,'wt',wtt);
    RawSpecStates = struct('s',st,'f',ft,'t',tt);
    SpecStates = struct('WhitenedSpecStates',WhitenedSpecStates,'RawSpecStates',RawSpecStates);
    save(strcat(OutPathSpectrogram,fn{1},'-',fn{2},'-',rn{1},'-SpecStates.mat'),'SpecStates')
    %% plot spectrogram
    Fmax = 50;
    NFmaxTheta = floor(Fmax*nFFTtheta/Fs);
    figure()
    h= pcolor(tt,ft(1:NFmaxTheta),log10(abs(st(:,1:NFmaxTheta)')));

    ylim([Fs*6/nFFTtheta Fmax])
    caxis([min(min(log10(abs(st(:,1:NFmaxTheta)'))))+1.9 max(max(log10(abs(st(:,1:NFmaxTheta)'))))-0.3])
    set(gca,'YScale','log')
%     set(gca,'YDir','normal')
%     colormap('jet')
    set(h, 'EdgeColor', 'none')
%     savefig(strcat(OutPathSpectrogram,'figures/',fn{1},'-',fn{2},'-',rn{1},'-Spectrogram.fig'));
   
%     ax2 = subplot(3,1,2); %plot sig AC
%     plot((Sinit+(1:1:length(LfpGeom(Sinit:Sfin,1))))/LFPfs, LfpGeom(Sinit:Sfin,ChIndHVS))
% 
%     ax3 = subplot(3,1,3); %plot sigDC
%     plot((SinitDC+(1:1:length(LfpGeomDC(SinitDC:SfinDC,1))))/FsDC, LfpGeomDC(SinitDC:SfinDC,ChIndHVS))
% 
%     set(h, 'EdgeColor', 'none')
%     linkaxes([ax1,ax2,ax3],'x')
%     xlim([Sinit/Fs Sfin/Fs])
