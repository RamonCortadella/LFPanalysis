function ICA2(FileName,T,indDB,extension,OutputPathStates, varargin)
%This function performs ICA on the LFP file indicated by FileName and using
%the metadata contained in the table T in the position indDB. Extension is
%used to specify the preprocessing stage of the Lfp file used (e.g.
%.interp.Lfp). OutputPathStates is the path to the states folder where the
%output ICA results are saved. Val is the list of files selected for the
%analysis. If val is provided Lfp and brain states are concatenated before
%computing ICA. NumICs indicates the number of ICs to extract. BrainState,
%if provided, limits the Lfp data analysed to the specified brain state (e.g. PerSWS, PerREM...).
%MinDurationStates discards brain state periods shorter that this value. 
%TimeScaleGlue adjust the slope of the sigmoid window boundaries used at
%the glue point between brain states and Lfp recordings. Freq bands is a
%cell array containing the labels of the frequency bands analysed. All the
%fields in the cell array are computed sequentially and the results saved
%together in a structure array. ThRoughness and interactiveDisplay
%('on'/'off') can be used in "groupstatistics" case to filter ICs with
%inhomogeneous coefficients interactively. 

[fMode, NumICs,BrainState,MinDurationStates,TimeScaleGlue,FreqBands,val,ThRoughness,interactiveDisplay] = DefaultArgs(varargin,{'compute',30,[],40,0.2,{'isa','delta'},[],1,'off'});
display(fMode)
display(NumICs)
%body of the function
switch fMode
    case 'compute'
        f.isa = [0.02 0.5];
        f.delta = [0.5 5];
        f.theta = [5 9];
        f.spindle = [9 20];
        f.beta = [20 30];
        f.lowgamma = [30 60];
        f.highgamma = [60 180];
        
        display(['loading ' FileName])        
        display(NumICs,'Number of ICs')
        display(FreqBands)
        
        nCh = T.NumCh(indDB);
        Fs = T.Fs(indDB);      
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB);
        
        if val
            Lfp=[];
            for i = 1:length(val)
                %force use of interpolated data and devices with same
                %number of channels, Fs, nCols and nRows
                % deal with brain states
                if contains(extension,'interp')
                    FileName = ['../../../data/DB/Files/' val{i} '/' val{i} extension];
                    display(FileName, 'loading file for Lfp concatenation across recordings')
                    LfpT = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
                    Lfp = [Lfp; LfpT];
                    
                    if BrainState
                        display(indDB)
                        display(T.RecordingId{indDB})
                        display(strcat(OutputPathStates,T.RecordingId{indDB},'-PerStates.mat'))
                        PerStates = load(strcat(OutputPathStates,T.RecordingId{indDB},'-PerStates.mat'));
                        SelectedPeriods = PerStates.PerStates.(BrainState);
                        DurationStates = SelectedPeriods(:,2)-SelectedPeriods(:,1);
                        SelectedPeriods = SelectedPeriods(find(DurationStates >=MinDurationStates),:); 
            %             DurationStates = DurationStates(find(DurationStates >=MinDurationStates));
                    end
                else
                    display('Extension must contain "interp" processing step for comparable recordings')
                    break
                end
            end
        else
       
            Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
            BadChannels = str2num(T.badChannels{indDB});
            if contains(FileName,'interp')
                BadChannels = [];
            end

            if BadChannels
                Lfp(:,BadChannels) = 0*Lfp(:,BadChannels);
            end

            if BrainState
                display(indDB)
                display(T.RecordingId{indDB})
                display(strcat(OutputPathStates,T.RecordingId{indDB},'-PerStates.mat'))
                PerStates = load(strcat(OutputPathStates,T.RecordingId{indDB},'-PerStates.mat'));
                SelectedPeriods = PerStates.PerStates.(BrainState);
                DurationStates = SelectedPeriods(:,2)-SelectedPeriods(:,1);
                SelectedPeriods = SelectedPeriods(find(DurationStates >=MinDurationStates),:); 
    %             DurationStates = DurationStates(find(DurationStates >=MinDurationStates));
            end
        end
        for iFreq = 1:length(FreqBands)
            LfpFilt = ButFilter(Lfp,2,f.(FreqBands{iFreq})/(Fs)*2,'bandpass');
            
            display(FreqBands{iFreq})
            
            LfpState = [];
            %concatenate Lfp for selected periods using a zeroing
            for i = 1:length(SelectedPeriods(:,1))
                LfpStateTemp = LfpFilt(floor(Fs*SelectedPeriods(i,1)):floor(Fs*SelectedPeriods(i,2)),:); 

                WindowTemp = ones([length(LfpStateTemp(:,1)),1]);
                Window1 = sigmoid([1:length(LfpStateTemp(:,1))]/Fs-10*TimeScaleGlue,TimeScaleGlue);
                Window= WindowTemp.*Window1'.*flip(Window1');
                Window = repmat(Window,1,nCh);
                display(size(LfpStateTemp))
                display(size(Window))
                LfpStateTemp = LfpStateTemp.*Window;


                LfpState = [LfpState; LfpStateTemp];
            end
            [ICAsig,A,W] = fastica(LfpState','numOfIC',NumICs);%, 'interactivePCA','off'); % [ICAsig,A,W]  A = mixing matrix (X = A.S), W = demixing matrix (A^-1) (S = WX)


    %       save identified sources and loadings
            out.(FreqBands{iFreq}).ICAsig= ICAsig;
            out.(FreqBands{iFreq}).A = A;
            out.(FreqBands{iFreq}).W = W;
            out.(FreqBands{iFreq}).Freq = f.(FreqBands{iFreq});
            
            
        end
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end

        save(strcat(OutputPathStates,T.FileName{indDB}(1:end-4),Coupling,'.ICA.mat'),'out')
        
    case 'display'
        
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end
        
        load(strcat(OutputPathStates,T.FileName{indDB}(1:end-4),Coupling,extension(1:end-4),'.mat'),'out');
        for iFreq = 1:length(FreqBands)
            
            A = out.(FreqBands{iFreq}).A;
            display(length(A(1,:)),'Index of IC in the displayed computation: CHECK FOR INCONSISTENCIES!!!')
            
            nRows = T.nRows(indDB);
            nCols = T.nCols(indDB);
            Fs = T.Fs(indDB);

            %display components map
            [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
            v = reshape(bX,[],1);
            w = reshape(bY,[],1);
            figure()
            title(FreqBands{iFreq})
            for i = 1:NumICs
                subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
                F=  scatteredInterpolant(v,w,A(:,i));
                bF = F(bX,bY);
                pcolor(bF)
    %             caxis([-0.3 0.3])
                colorbar

                title(['index=' num2str(i)])
            end
            
            
            
            figure()
            title(FreqBands{iFreq})
            hold on
            nfft =2^14;
            nperov = 2^12;
            for indIC = 1:length(out.(FreqBands{iFreq}).ICAsig(:,1))
                subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,indIC)

                [ppx, f] = pwelch(out.(FreqBands{iFreq}).ICAsig(indIC,:),nfft,nperov,nfft,Fs);
                plot(f,ppx)
                set(gca,'XScale','log')
                set(gca,'YScale','log')
            end
            
            

    %         tRatio = SpecStates.RawSpecStates.t;
    %         tICA = (1:length(out.ICAsig(1,:)))/Fs;

    %         figure()
    %         hold on
    %         for indIC = 1:length(out.ICAsig(:,1))
    %             subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,indIC)
    %             
    %             display(indIC,'computing coherence for IC with index')
    %             TandIC = [tICA',out.ICAsig(indIC,:)'];
    %             ICAsig2 = Interpolate(TandIC,tRatio,'trim','off');
    %             ICAsig2(isnan(ICAsig2))=0;
    %             
    %             [cxy,f] = mscohere(ICAsig2(:,2),PerStates.bands.ratio,hamming(100),80,100,1/(tRatio(2)-tRatio(1)));
    %             plot(f,cxy)
    %             set(gca,'XScale','log')
    %             set(gca,'YScale','linear')
    %             set(gca,'Ylim',[0 0.3])
    %             set(gca,'Xlim',[0.005 0.5])
    %         end

    %         figure()
    %         ax1 = subplot(3,1,1);
    %         TandIC = [tICA',out.ICAsig(33,:)'];
    %         ICAsig2 = Interpolate(TandIC,tRatio,'trim','off');
    %         ICAsig2(isnan(ICAsig2))=0;
    %         [cxy,f] = mscohere(ICAsig2(:,2),PerStates.bands.ratio,hamming(100),80,100,1/(tRatio(2)-tRatio(1)));
    %         plot(f,cxy)
    %         set(gca,'XScale','log')
    %         set(gca,'YScale','linear')
    %         set(gca,'Ylim',[0 0.3])
    %         set(gca,'Xlim',[0.005 0.5])
    %         
    %         ax2 = subplot(3,1,2);
    %         plot(tRatio,ICAsig2(:,2),'b')
    %         ax3 = subplot(3,1,3);
    %         plot(tRatio,PerStates.bands.ratio,'r')
    %         
    %         linkaxes([ax1, ax2, ax3],'x');

            %dim A = (nCh rows x NumICs cols)
            %plot coefficients (loadings)
            
            SelectICs = input('indices of ICs to further display');
            [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
            v = reshape(bX,[],1);
            w = reshape(bY,[],1);


            for i = SelectICs
                figure()
                subplot(1,3,1)
                F=  scatteredInterpolant(v,w,A(:,i));
                bF = F(bX,bY);
                imagesc(bF)
                colorbar
                subplot(1,3,[2,3])
                tICA = (1:length(out.(FreqBands{iFreq}).ICAsig(i,:)))/Fs;
                plot(tICA,out.(FreqBands{iFreq}).ICAsig(i,:))
                title(FreqBands{iFreq})
            end
        end
       
    case 'groupstats'
        for fN = val
            load(strcat(OutputPathStates,T.FileName{indDB}(end-4:end),'-ICA.mat'),'out');
            load(strcat(OutputPathStates,T.RecordingId{indDB},'-PerStates.mat'),'PerStates');
            load(strcat(OutputPathStates,T.RecordingId{indDB},'-SpecStates.mat'),'SpecStates');
        end
%             tRatio = SpecStates.RawSpecStates.t;
%             tICA = (1:length(out.ICAsig(1,:)))/T.Fs;
%             
%             for indIC = 1:length(ICAsig(:,1))
%                 TandIC = [tICA';ICAsig(indIC,:)'];
%                 ICAsig2 = Interpolate(TandIC,tRatio,'trim','off');
%                 [cxy,f] = mscohere(ICAsig2(:,2),PerStates.bands.ratio,[],[],[],1/(tRatio[2]-tRatio[1]));
%             strcat([InputPath, fN(1:end-4),extension])
%         end
        % load and then compute stats and output 
        
        
        
%         if strcmp(interactiveDisplay,'on')
%                 Ac = [];
% 
%             while true
%                 clf('reset')
%                 title(FreqBands{iFreq})
% 
%                 loop = input('(re)-compute? (1-yes/0-no)');
%                 if loop ==1
%                     for i = 1:NumICs %find good components in terms of smoothness
% 
%                         SmInd = abs(sum(diff(A(:,i)))/mean(abs(A(:,i))));
%                         display(SmInd)
%                         if SmInd>= ThRoughness
%                             Ac = [Ac, A(:,i)];
%                         subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
%                         F=  scatteredInterpolant(v,w,A(:,i));
%                         bF = F(bX,bY);
%                         pcolor(bF)
%                         colorbar
%                         end
%                     end
%                 end
%                 ThRoughness = input('Change in roughness index threshold for IC discard');
% 
%             end
%         end
end

