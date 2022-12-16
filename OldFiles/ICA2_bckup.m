function ICA2(FileName,T,indDB,OutputPathStates,extension, varargin)
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
%inhomogeneous coefficients interactively. PhAmpKlustersAnimal indicates
%the index of the animal for which the Ph-Amp coupling is evaluated, if
%only the significant klusters are taken into account

[fMode, NumICs,BrainState,val,MinDurationStates,TimeScaleGlue,FreqBands,ThRoughness,interactiveDisplay,PhAmpKlustersAnimal] = DefaultArgs(varargin,{'compute',30,[],[],40,0.2,{'isa','delta','theta','spindle','beta','lowgamma','highgamma'},1,'off',3}); %'isa','delta','theta','spindle','beta','lowgamma','highgamma'
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
        f.lowgamma = [30 50];
        f.highgamma = [50 100];
           
        display(FreqBands)
        
        nCh = T.NumCh(indDB);
        Fs = T.Fs(indDB);      
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB);
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end
        
        if strcmp(Coupling, '.AC') | strcmp(Coupling, '')
            if T.depth(indDB) == 1 & T.SingleShank(indDB) == 0
                DownFact = 6;
            elseif T.depth(indDB)==0 
                DownFact = 3;
            end
        end
        Fs = Fs/DownFact;
        display(Fs)
        szVal = size(val);
        if szVal(2)>=1
            InputPathState = '/storage2/ramon/data/DB/Files/';
            Lfp=[];
            
            SelectedPeriods = [];
            for i = 1:length(val)
                %force use of interpolated data and devices with same
                %number of channels, Fs, nCols and nRows
                % deal with brain states
                Tend = 0;
                if contains(extension,'interp')
                    dn = split(val{1},'-');
                    Device = [dn{1} '-' dn{2}];
                    if contains(T.FileName(indDB),'ECoG')
                        FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.ECoG' Coupling extension];
                    elseif  contains(T.FileName(indDB),'Depth')
                        FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} '.Depth' Coupling extension];
                    else
                        FileName = ['../../../data/DB/Files/' Device '/' val{i} '/' val{i} Coupling extension];
                    end
                    
                    display(FileName, 'loading file for Lfp concatenation across recordings')
                    LfpT = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
                    LfpT = ButFilter(LfpT,2,1/DownFact,'low');
                    
                    LfpT = LfpT(1:DownFact:end,:);
                    szLfp = size(Lfp);
                    Tend = szLfp(1)/Fs;
                    
                    Lfp = [Lfp; LfpT];
                    
                    if BrainState
                        display(indDB)
                        display(T.RecordingId{indDB})
                        
                        display(strcat(InputPathState,Device,'/',val{i},'/States/',val{i},'.PerStates.mat'))
                        
                        PerStates = load(strcat(InputPathState,Device,'/',val{i},'/States/',val{i} ,'.PerStates.mat'));
                        SelectedPeriodsTemp = PerStates.PerStates.(BrainState);
                        DurationStates = SelectedPeriodsTemp(:,2)-SelectedPeriodsTemp(:,1);
                        SelectedPeriodsTemp = SelectedPeriodsTemp(find(DurationStates >=MinDurationStates),:); 
                        SelectedPeriodsTemp = SelectedPeriodsTemp + Tend;
                        
                        SelectedPeriods = [SelectedPeriods; SelectedPeriodsTemp];
            %             DurationStates = DurationStates(find(DurationStates >=MinDurationStates));
                    end
                else
                    display('ATTENTION! : Extension must contain "interp" processing step for comparable recordings!!!!! This process will crash in 100s')
                    pause(100)
                    break
                end
            end
        else
            display(FileName, 'Loading file name')
            Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
            Lfp = ButFilter(Lfp,2,1/DownFact,'low');
                    
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
                display(strcat(OutputPathStates,T.RecordingId{indDB},'.PerStates.mat'))
                PerStates = load(strcat(OutputPathStates,T.RecordingId{indDB},'.PerStates.mat'));
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
                display(size(LfpFilt),'size Lfp')
                display(floor(Fs*SelectedPeriods(i,:)), 'min and max samp Period')
                LfpStateTemp = LfpFilt(floor(Fs*SelectedPeriods(i,1)+1):floor(Fs*SelectedPeriods(i,2))-1,:); 

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
            out.Freqs = FreqBands;
            
        end
        
        if szVal(2)>=1
            save(strcat(OutputPathStates,'Group',Coupling,'.',BrainState,'.ICA.mat'),'out', '-v7.3')
        else
            save(strcat(OutputPathStates,T.RecordingId{indDB},Coupling,'.',BrainState,'.ICA.mat'),'out', '-v7.3')
        end
        
    case 'display'
        
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end
        
        szVal = size(val);
        if szVal(2)>=1
            load(strcat(OutputPathStates,'Group',Coupling,'.',BrainState,'.ICA.mat'),'out')
        else
            load(strcat(OutputPathStates,T.RecordingId{indDB},Coupling,'.',BrainState,'.ICA.mat'),'out')
        end
        
        for iFreq = 1:length(FreqBands)
            
            A = out.(FreqBands{iFreq}).A;
            display(length(A(1,:)),'Index of IC in the displayed computation: CHECK FOR INCONSISTENCIES!!!')
            
            nRows = T.nRows(indDB);
            nCols = T.nCols(indDB);
            Fs = T.Fs(indDB);
            
            if strcmp(Coupling, '.AC') | strcmp(Coupling, '')
                if T.depth(indDB) == 1 & T.SingleShank(indDB) == 0
                    DownFact = 6;
                elseif T.depth(indDB)==0 
                    DownFact = 3;
                end
            end
            Fs = Fs/DownFact;
            
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
                imagesc(bF)
    %             caxis([-0.3 0.3])
                colorbar

                title(['index=' num2str(i)])
            end
            
            
            
%             figure()
%             title(FreqBands{iFreq})
%             hold on
%             nfft =2^14;
%             nperov = 2^12;
%             for indIC = 1:length(out.(FreqBands{iFreq}).ICAsig(:,1))
%                 subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,indIC)
% 
%                 [ppx, f] = pwelch(out.(FreqBands{iFreq}).ICAsig(indIC,:),nfft,nperov,nfft,Fs);
%                 plot(f,ppx)
%                 set(gca,'XScale','log')
%                 set(gca,'YScale','log')
%             end
            
            

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
     

    case 'Significance'
        
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end    
        
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB); 
        nCh = T.NumCh(indDB);    
        
        dev = split(val{1},'-');
        dev = [dev{1} '-' dev{2}];
        WithinGroup = [];
        for i = 1:length(val)
            WithinGroup = [WithinGroup , contains(val{i},dev)];
        end
        
        if strcmp(Coupling,'')
            CouplingCell = {''};
        else
            CouplingCell = {'.DC','.AC'};
        end
        dev0 = '';
        counter = 0;
        
        for i = 1:length(val)
            dev = split(val{i},'-');
            dev = [dev{1} '-' dev{2}];
            if length(find(WithinGroup == 1)) == length(val) %across sessions in one animal
                
                for Coupling = CouplingCell
                    load(strcat(OutputPathStates,dev,'/',val{i},'/States/',val{i},Coupling,'.',BrainState,'.ICA.mat'),'out')  
                    if strcmp(Coupling,'.DC')
                        out0 = out;
                    else
                        fn = fieldnames(out);
                        for iF = 1:length(fieldnames(out))
                            out0.(fn{iF}) = out.(fn{iF});
                        end
                    end
                end
            else
                counter = counter+1;
                if strcmp(dev0,dev)
                    continue
                else
                    for iC = 1:length(CouplingCell)
                        Coupling = CouplingCell{iC};
                        display(strcat(OutputPathStates,dev,'/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'))
                        load(strcat(OutputPathStates,dev,'/GroupsAnalysis/','Group',Coupling,'.',BrainState,'.ICA.mat'),'out')
                        if strcmp(Coupling,'.DC')
                            out0 = out;
                        else
                            fn = fieldnames(out);
                            for iF = 1:length(fieldnames(out))
                                out0.(fn{iF}) = out.(fn{iF});
                            end
                        end
                    end
                end
            end
            display(fieldnames(out0),'fieldnames out0 **********')
            
            dev0 = dev;
            for iFreq = 1:length(FreqBands)
                display(size(out0.(FreqBands{iFreq}).A))
                if i == 1
                    groupA.(FreqBands{iFreq}) = zeros([nCh,NumICs]); %size predefinition and indexing below needed in those cases where the number of IC is below target because of convergence issues
                    groupA.(FreqBands{iFreq})(:,1:length(out0.(FreqBands{iFreq}).A(1,:))) = out0.(FreqBands{iFreq}).A; 
                else
                    A = zeros([nCh,NumICs]);
                    A(:,1:length(out0.(FreqBands{iFreq}).A(1,:))) = out0.(FreqBands{iFreq}).A;
                    groupA.(FreqBands{iFreq}) = cat(3,groupA.(FreqBands{iFreq}),A);
                end
            end
        end
        
        for iFreq = 1:length(FreqBands)
            [out2.clustering.(FreqBands{iFreq}),~,out2.linkpvalues.(FreqBands{iFreq})] = isctest(groupA.(FreqBands{iFreq}),0.05,0.05,'mixing');
        end
        
        [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
        v = reshape(bX,[],1);
        w = reshape(bY,[],1);
        
        kl = out2.clustering.(FreqBands{4});
        nkl = length(kl(:,1,1));
        for ikl = 1:nkl
            figure()
            nk=length(kl(ikl,:,1));
            for k = 1:nk
                if kl(ikl,k,1) ~= 0

                    subplot(floor(sqrt(nk))+1,floor(sqrt(nk))+1,k)
                    F=  scatteredInterpolant(v,w,groupA.(FreqBands{4})(:,kl(ikl,k,1),k));
                    bF = F(bX,bY);
                    imagesc(bF)
                    colorbar
                    title([num2str(kl(ikl,k,1)) 'th component of' num2str(k) 'th subject in kluster ' num2str(ikl) 'th'])
                end
            end
        end
        
        if length(find(WithinGroup == 1)) == length(val)
            save(strcat(OutputPathStates,dev,'/GroupsAnalysis/','GroupKlusters','.ICA.mat'),'out2')
        else
            display(fieldnames(out2.clustering),'fieldnames **')
            save(strcat(OutputPathStates,'/GroupsAnalysis/','GroupKlusters','.ICA.mat'),'out2')
        end
%           
  
    case 'PhAmpCoupling'
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end
        
        if strcmp(Coupling,'')
            CouplingCell = {''};
        else
            CouplingCell = {'.DC','.AC'};
        end
        for iC = 1:length(CouplingCell)
            Coupling = CouplingCell{iC};
            load(strcat(OutputPathStates,'Group',Coupling,'.',BrainState,'.ICA.mat'),'out'); 
            if strcmp(Coupling,'.DC')
                out0 = out;
            else
                fn = fieldnames(out);
                for iF = 1:length(fieldnames(out))
                    out0.(fn{iF}) = out.(fn{iF});
                end
            end
        end
        
        Fs = T.Fs(indDB);
        Fs1 = 65.10416667;
        if strcmp(Coupling, '.AC') | strcmp(Coupling, '')
            if T.depth(indDB) == 1 & T.SingleShank(indDB) == 0
                DownFact = 6;
            elseif T.depth(indDB)==0 
                DownFact = 3;
            end
        end
        
        Fs2 = 651.0416667/DownFact;
%         Fs = Fs/DownFact;
        
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB);        
        
        PhAmpC_amp=[];
        PhAmpC_angle=[];
        
        if PhAmpKlustersAnimal
            display(strcat('../../../data/DB/Files/GroupsAnalysis/','GroupKlusters','.ICA.mat'))
            load(strcat('../../../data/DB/Files/GroupsAnalysis/','GroupKlusters','.ICA.mat'),'out2')
            out1 = out2;
            display(fieldnames(out1.clustering))
            for iFreq = 1:length(fieldnames(out1.clustering))
                counter = 0;
                for ik = 1:length(out1.clustering.(FreqBands{iFreq})(:,1))
                    indT = find(out1.clustering.(FreqBands{iFreq})(ik,:)~=0);
                    if length(indT)== length(out1.clustering.(FreqBands{iFreq})(1,:))
                        counter = counter+1;
                        if counter ==1
                            goodKList.(FreqBands{iFreq}) = ik;
                        else
                            goodKList.(FreqBands{iFreq}) = [goodKList.(FreqBands{iFreq}), ik];
                        end
                    end
                end
            end
        end
        
        if length(fieldnames(out0)) >=2
            for IndLoopF = 1:length(fieldnames(out0))-1
                for iFreq = 1:length(fieldnames(out0))-IndLoopF
                    display(iFreq,'frequency index Ph-Amp')
                    display(IndLoopF,'Inter frequency step')
                    for i = 1:length(out0.(FreqBands{iFreq}).ICAsig(:,1))
                        display(i,'IC index of low freqs')
                        
                        for i2 = 1:length(out0.(FreqBands{iFreq+IndLoopF}).ICAsig(:,1))
% --------------------- implementation 1
%                             pow = abs(hilbert(out.(FreqBands{iFreq+IndLoopF}).ICAsig(i2,:)));
%                             ph = angle(hilbert(out.(FreqBands{iFreq}).ICAsig(i,:)));
%                             powsum = sum(pow);
%                             R = sum(exp(1i*ph).*pow)./powsum;
%                             PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = abs(R);
%                             PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = angle(R);
                            
% --------------------- implementation 2
                            
                            if PhAmpKlustersAnimal
                                if find(out1.clustering.(FreqBands{iFreq})(goodKList.(FreqBands{iFreq}),PhAmpKlustersAnimal) == i) & find(out1.clustering.(FreqBands{iFreq+IndLoopF})(goodKList.(FreqBands{iFreq+IndLoopF}),PhAmpKlustersAnimal)==i2)
                                    if (iFreq == 1 | iFreq == 2) & iFreq+IndLoopF >= 3 
                                        
                                        FNi = Fs2/2;
                                        fPow = [out0.(FreqBands{iFreq+IndLoopF}).Freq/FNi; (out0.(FreqBands{iFreq+IndLoopF}).Freq+1)/FNi]; %the frequency band in the second position is only provided to temporarily avoid a bug//it is not used
                                        fPh = [out0.(FreqBands{iFreq}).Freq/FNi; (out0.(FreqBands{iFreq}).Freq+1)/FNi]; 
                                        display(out0.(FreqBands{iFreq+IndLoopF}).Freq,'*')
                                        display(out0.(FreqBands{iFreq}).Freq,'*')
                                        
                                        size([[1:length(out0.(FreqBands{iFreq}).ICAsig(i,:))]'/Fs1 out0.(FreqBands{iFreq}).ICAsig(i,:)'])
                                        
                                        ICAsigInter = Interpolate([[1:length(out0.(FreqBands{iFreq}).ICAsig(i,:))]'/Fs1 out0.(FreqBands{iFreq}).ICAsig(i,:)'],[1:length(out0.(FreqBands{iFreq+IndLoopF}).ICAsig(i,:))]/Fs2,'trim','off');
                                        ICAsigInter = ICAsigInter(:,2)';
                                        
                                    elseif iFreq == 1 & iFreq+IndLoopF ==2
                                        FNi = Fs1/2;
                                        fPow = [out0.(FreqBands{iFreq+IndLoopF}).Freq/FNi; (out0.(FreqBands{iFreq+IndLoopF}).Freq+1)/FNi]; %the frequency band in the second position is only provided to temporarily avoid a bug//it is not used
                                        fPh = [out0.(FreqBands{iFreq}).Freq/FNi; (out0.(FreqBands{iFreq}).Freq+1)/FNi]; 
                                        display(out0.(FreqBands{iFreq+IndLoopF}).Freq,'**')
                                        display(out0.(FreqBands{iFreq}).Freq,'**')
                                        
                                        ICAsigInter = out0.(FreqBands{iFreq}).ICAsig(i,:);
                                        
                                    elseif iFreq >= 3 & iFreq+IndLoopF >=3
                                        FNi = Fs2/2;
                                        fPow = [out0.(FreqBands{iFreq+IndLoopF}).Freq/FNi; (out0.(FreqBands{iFreq+IndLoopF}).Freq+1)/FNi]; %the frequency band in the second position is only provided to temporarily avoid a bug//it is not used
                                        fPh = [out0.(FreqBands{iFreq}).Freq/FNi; (out0.(FreqBands{iFreq}).Freq+1)/FNi]; 
                                        
                                        display(out0.(FreqBands{iFreq+IndLoopF}).Freq,'**')
                                        display(out0.(FreqBands{iFreq}).Freq,'**')
                                        
                                        
                                        ICAsigInter = out0.(FreqBands{iFreq}).ICAsig(i,:);
                                    end
                                    display(size(ICAsigInter))
                                    display(size(out0.(FreqBands{iFreq+IndLoopF}).ICAsig(i2,:)'))
                                    display(fPow*FNi)
                                    display(fPh*FNi)
                                    [out2, ~, ~] = PowerPhasePairsR(ICAsigInter',  fPh, out0.(FreqBands{iFreq+IndLoopF}).ICAsig(i2,:)' , fPow, 1,'but',@PowerModulation);
                                    PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = out2.Ramp(1,1);
                                    PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = out2.phbins(find(sq(out2.pow_dens(:,1,1))==max(sq(out2.pow_dens(:,1,1)))),1,1);
                                else
                                    continue
                                end
                            else
                                [out2, ~, ~] = PowerPhasePairsR(out0.(FreqBands{iFreq}).ICAsig(i,:)',  fPh, out0.(FreqBands{iFreq+IndLoopF}).ICAsig(i2,:)' , fPow, 1,'but',@PowerModulation);
                                PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = out2.Ramp(1,1);
                                PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = out2.phbins(find(sq(out2.pow_dens(:,1,1))==max(sq(out2.pow_dens(:,1,1)))),1,1);
                            end
                        end
                    end

                    
                    figure()
                    imagesc(PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}]))
                    title([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])
                    colormap('parula')
                    colorbar()
                    figure()
                    imagesc(PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}]))
                    title([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])
                    colormap('hsv')
                    colorbar()
                end
            end
        else
            display('Phase-amplitude coupling not computed because only one frequency provided')
        end
        while true
            AnglularPlotIndices = input('IN BRAKETS! [ind of the ph-freq, steps between ph and amp freq, IC of ph, IC of amplitude');

            [out2, ~, ~] = PowerPhasePairsR(out.(FreqBands{AnglularPlotIndices(1)}).ICAsig(AnglularPlotIndices(3),:)',  fPh, out.(FreqBands{AnglularPlotIndices(1)+AnglularPlotIndices(2)}).ICAsig(AnglularPlotIndices(4),:)' , fPow, 1,'but',@PowerModulation);
            
            
            [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
            v = reshape(bX,[],1);
            w = reshape(bY,[],1);
            figure()
            subplot(1,3,1)
            F=  scatteredInterpolant(v,w,out.(FreqBands{AnglularPlotIndices(1)}).A(:,AnglularPlotIndices(3)));
            bF = F(bX,bY);
            imagesc(bF)
            colorbar
            title('IC of ph-frequency')
            subplot(1,3,2)
            F=  scatteredInterpolant(v,w,out.(FreqBands{AnglularPlotIndices(1)+AnglularPlotIndices(2)}).A(:,AnglularPlotIndices(4)));
            bF = F(bX,bY);
            imagesc(bF)
            colorbar
            title('IC of pow-frequency')
            subplot(1,3,3)
            polarplot(out2.phbins(:,1,1),sq(out2.pow_dens(:,1,1))'); axis tight; box off
        end
end


