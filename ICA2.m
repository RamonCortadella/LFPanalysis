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

[fMode, NumICs,BrainState,val,FreqBands,DevType,FreqBandsHigh,CouplingLow,CouplingHigh,MinDurationStates,TimeScaleGlue,ThRoughness,interactiveDisplay,PhAmpKlustersAnimal] = DefaultArgs(varargin,{'compute',30,[],[],{'isa','delta','theta','spindle','beta','lowgamma','highgamma'},'',{'delta'},'.DC','.AC',40,0.2,1,'off',[]}); %'isa','delta','theta','spindle','beta','lowgamma','highgamma'
display(fMode)
display(NumICs)
display(FreqBands)
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
        
        % Gets AC or DC coupling from a single index of the files
        % selection. Thus it works on groups of specific DC/AC coupled data
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end
        
        % if AC coupled data it is downsampled 
        DownFact = [];
        if strcmp(Coupling, '.AC') | strcmp(Coupling, '')
            if T.depth(indDB) == 1 & T.SingleShank(indDB) == 0
                DownFact = 6;
            elseif T.depth(indDB)==0 
                DownFact = 3;
            end
            
            Fs = Fs/DownFact;
        end
        
        display(Fs)
        
        
        szVal = size(val);
        if szVal(2)>=1 %if val is provided multiple LFP files are concatenated, if not, only one is taken specified by FileName
            InputPathState = '/storage2/ramon/data/DB/Files/';
            Lfp=[];
            
            SelectedPeriods = [];
            for i = 1:length(val)
                
                Tend = 0;
                if contains(extension,'interp') %force use of interpolated files
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
                    if DownFact 
                        LfpT = ButFilter(LfpT,2,1/DownFact,'low');
                        LfpT = LfpT(1:DownFact:end,:);
                    end
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
            out.val = val;
            out.MinDurationStates = MinDurationStates;
            
        end
        
        if szVal(2)>=1
            save(strcat(OutputPathStates,'Group',Coupling,'.',BrainState,'.ICA.mat'),'out', '-v7.3')
        else
            save(strcat(OutputPathStates,T.RecordingId{indDB},Coupling,'.',BrainState,'.ICA.mat'),'out', '-v7.3')
        end
    

    case 'display'
        % the val must be generated accordingly to the data because indDB
        % is used to find the DC/AC coupling
        if T.CoupledAC(indDB) == 1 & T.CoupledDC(indDB) ==0
            Coupling ='.AC';
        elseif T.CoupledAC(indDB) == 0 & T.CoupledDC(indDB) ==1
            Coupling ='.DC';
        else
            Coupling ='';
        end
        
        szVal = size(val);
        if szVal(2)>=1
            load(strcat(OutputPathStates,'Group',Coupling,'.',BrainState,DevType,'.ICA.mat'),'out')
        else
            load(strcat(OutputPathStates,T.RecordingId{indDB},Coupling,'.',BrainState,DevType,'.ICA.mat'),'out')
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
                Fs = Fs/DownFact;
            
            end
            
            %display components map
            [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
            v = reshape(bX,[],1);
            w = reshape(bY,[],1);
            figure()
            title(FreqBands{iFreq})
            for i = 1:NumICs
                subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
                if T.depth(indDB)==1 & T.SingleShank(indDB)==1
                    imagesc(A(:,i))
                else
                    F=  scatteredInterpolant(v,w,A(:,i));
                    bF = F(bX,bY);
                    imagesc(bF)
                end
    %             caxis([-0.3 0.3])
                colorbar

                title(['index=' num2str(i)])
            end
            
            
            SelectICs = input('indices of ICs to further display');
            [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
            v = reshape(bX,[],1);
            w = reshape(bY,[],1);


            for i = SelectICs
                figure()
                subplot(1,3,1)
                if T.depth(indDB)==1 & T.SingleShank(indDB)==1
                    imagesc(A(:,i))
                else
                    F=  scatteredInterpolant(v,w,A(:,i));
                    bF = F(bX,bY);
                    imagesc(bF)
                end
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
        
        %if val contains more than one device the statistical test is run
        %across devices (i.e. animals). If not it is computed across
        %recordings in the indicated device/animal
        
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
                    load(strcat(OutputPathStates,dev,'/',val{i},'/States/',val{i},Coupling,'.',BrainState,DevType,'.ICA.mat'),'out')  
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
                        display(strcat(OutputPathStates,dev,'/GroupsAnalysis/','Group',Coupling,'.',BrainState,DevType,'.ICA.mat'))
                        load(strcat(OutputPathStates,dev,'/GroupsAnalysis/','Group',Coupling,'.',BrainState,DevType,'.ICA.mat'),'out')
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
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB); 
        nCh = T.NumCh(indDB);    
        
        
        iFreq = 1; %for the Ph-amp coupling only the first band in FreqBands and the first in FreqBandsHigh are taken
        
        load(strcat(OutputPathStates,'Group',CouplingLow,'.',BrainState,DevType,'.ICA.mat'),'out'); 
        outL = out;
        AL = outL.(FreqBands{1}).A;
        load(strcat(OutputPathStates,'Group',CouplingHigh,'.',BrainState,DevType,'.ICA.mat'),'out'); 
        outH = out;
        AH = outH.(FreqBandsHigh{1}).A;
        
        NumICsL = length(outL.(FreqBands{1}).ICAsig(:,1));
        NumICsH = length(outH.(FreqBandsHigh{1}).ICAsig(:,1));
        %display components map
        [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
        v = reshape(bX,[],1);
        w = reshape(bY,[],1);
        figure()
        title(FreqBands{iFreq})
        for i = 1:NumICsL
            subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
            F=  scatteredInterpolant(v,w,AL(:,i));
            bF = F(bX,bY);
            imagesc(bF)
%             caxis([-0.3 0.3])
            colorbar

            title(['index=' num2str(i)])
        end
        figure()
        title(FreqBandsHigh{iFreq})
        for i = 1:NumICsH
            subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs))+1,i)
            F=  scatteredInterpolant(v,w,AH(:,i));
            bF = F(bX,bY);
            imagesc(bF)
%             caxis([-0.3 0.3])
            colorbar

            title(['index=' num2str(i)])
        end
        
        indicesIC = input('indices of low frequency ICs [in brackets]');        
        indicesIC2 = input('indices of high frequency ICs [in brackets]');
        
        Fs = T.Fs(indDB);
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB);        
        
        PhAmpC_amp=[];
        PhAmpC_angle=[];
        PhAmpC_amp_PowSc=[];
        PhAmpC_angle_PowSc=[];
        
        if strcmp(CouplingLow,'.DC') & strcmp(CouplingHigh,'.AC')
            Fs = T.Fs(indDB-1);
            FsLow = T.Fs(indDB);
%             if strcmp(Coupling, '.AC') | strcmp(Coupling, '')
            if T.depth(indDB) == 1 & T.SingleShank(indDB) == 0
                DownFact = 6;
            elseif T.depth(indDB)==0 
                DownFact = 3;
            end
            Fs = Fs/DownFact;
%             end
        end
        
        display(iFreq,'frequency index Ph-Amp')
        outL.(FreqBands{iFreq}).ICAsigInt=[];
        for i = indicesIC %1:length(out.(FreqBands{iFreq}).ICAsig(:,1))
            display(i,'IC index of low freqs')

            for i2 = indicesIC2 %1:length(out.(FreqBands{iFreq+IndLoopF}).ICAsig(:,1))
% --------------------- implementation 1
%                             pow = abs(hilbert(out.(FreqBands{iFreq+IndLoopF}).ICAsig(i2,:)));
%                             ph = angle(hilbert(out.(FreqBands{iFreq}).ICAsig(i,:)));
%                             powsum = sum(pow);
%                             R = sum(exp(1i*ph).*pow)./powsum;
%                             PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = abs(R);
%                             PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBands{iFreq+IndLoopF}])(i,i2) = angle(R);

% --------------------- implementation 2
                FNi = Fs/2;
                fPow = [outH.(FreqBandsHigh{iFreq}).Freq/FNi; (outH.(FreqBandsHigh{iFreq}).Freq+1)/FNi]; %the frequency band in the second position is only provided to temporarily avoid a bug//it is not used
                fPh = [outL.(FreqBands{iFreq}).Freq/FNi; (outL.(FreqBands{iFreq}).Freq+1)/FNi]; 
                display(fPow)
                display(fPh)
                % interpolate DC coupled ICs
                if strcmp(CouplingLow,'.DC') & strcmp(CouplingHigh,'.AC')
                    size([[1:length(outL.(FreqBands{iFreq}).ICAsig(i,:))]'/FsLow, outL.(FreqBands{iFreq}).ICAsig(i,:)'])
                    size([1:length(outH.(FreqBandsHigh{iFreq}).ICAsig(1,:))]'/Fs)
                    sig = Interpolate([[1:length(outL.(FreqBands{iFreq}).ICAsig(i,:))]'/FsLow, outL.(FreqBands{iFreq}).ICAsig(i,:)'],[1:length(outH.(FreqBandsHigh{iFreq}).ICAsig(1,:))]'/Fs,'trim','off');
                    size(sig)
                    size(sig(:,2)')
                    outL.(FreqBands{iFreq}).ICAsigInt(i,:) = sig(:,2)';
                    size(outL.(FreqBands{iFreq}).ICAsigInt(i,:)')
                    size(outH.(FreqBandsHigh{iFreq}).ICAsig(i2,:)')
                    size(outL.(FreqBands{iFreq}).ICAsig(i,:)')
                    [out2, ~, ~] = PowerPhasePairsR(outL.(FreqBands{iFreq}).ICAsigInt(i,10:end-1000)',  fPh, outH.(FreqBandsHigh{iFreq}).ICAsig(i2,10:end-1000)' , fPow, 1,'but',@PowerModulation);
                    [out3, ~, ~] = PowerPhasePairsR2(outL.(FreqBands{iFreq}).ICAsigInt(i,10:end-1000)',  fPh, outH.(FreqBandsHigh{iFreq}).ICAsig(i2,10:end-1000)' , fPow, 1,'but',@PowerModulation_PowerScaled);
                else
               
                    [out2, ~, ~] = PowerPhasePairsR(outL.(FreqBands{iFreq}).ICAsig(i,:)',  fPh, outH.(FreqBandsHigh{iFreq}).ICAsig(i2,:)' , fPow, 1,'but',@PowerModulation);
                    [out3, ~, ~] = PowerPhasePairsR2(outL.(FreqBands{iFreq}).ICAsig(i,:)',  fPh, outH.(FreqBandsHigh{iFreq}).ICAsig(i2,:)' , fPow, 1,'but',@PowerModulation_PowerScaled);
                end
                PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}])(i,i2) = out2.Ramp(1,1);
                PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}])(i,i2) = out2.phbins(find(sq(out2.pow_dens(:,1,1))==max(sq(out2.pow_dens(:,1,1)))),1,1);
                PhAmpC_amp_PowSc.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}])(i,i2) = out3.Ramp(1,1);
                PhAmpC_angle_PowSc.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}])(i,i2) = out3.phbins(find(sq(out3.pow_dens(:,1,1))==max(sq(out3.pow_dens(:,1,1)))),1,1);

            end
        end


        figure()
        title([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}])
        subplot(1,2,1)
        imagesc(PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]))
        colormap('parula')
        colorbar()
        subplot(1,2,2)
        a = PhAmpC_amp.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]);
        imagesc(reshape(a(find(a~=0)),[length(indicesIC),length(indicesIC2)]))
        colormap('parula')
        colorbar()
        
        
        figure()
        subplot(1,2,1)
        imagesc(PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]))
        colormap('hsv')
        colorbar()
        subplot(1,2,2)
        a = PhAmpC_angle.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]);
        imagesc(reshape(a(find(a~=0)),[length(indicesIC),length(indicesIC2)]))
        colormap('hsv')
        colorbar()

        figure()
        subplot(1,2,1)
        imagesc(PhAmpC_amp_PowSc.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]))
        colormap('parula')
        colorbar()
        subplot(1,2,2)
        a = PhAmpC_amp_PowSc.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]);
        imagesc(reshape(a(find(a~=0)),[length(indicesIC),length(indicesIC2)]))
        colormap('parula')
        colorbar()
        title([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq} '-power weigthed'])
        
        
        figure()
        subplot(1,2,1)
        imagesc(PhAmpC_angle_PowSc.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]))
        colormap('hsv')
        colorbar()
        subplot(1,2,2)
        a = PhAmpC_angle_PowSc.([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq}]);
        imagesc(reshape(a(find(a~=0)),[length(indicesIC),length(indicesIC2)]))
        colormap('hsv')
        colorbar()
        title([FreqBands{iFreq} 'to' FreqBandsHigh{iFreq} '-power weigthed'])
          
        while true
            out2=struct();
            AngularPlotIndices = input('IN BRAKETS! [IC of ph, IC of amplitude]');
            display(AngularPlotIndices)

            [out2, ~, ~] = PowerPhasePairsR(outL.(FreqBands{iFreq}).ICAsig(AngularPlotIndices(1),10:end-1000)',  fPh, outH.(FreqBandsHigh{iFreq}).ICAsig(AngularPlotIndices(2),10:end-1000)' , fPow, 1,'but',@PowerModulation);
            
            display(continues)
            
            [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
            v = reshape(bX,[],1);
            w = reshape(bY,[],1);
            figure()
            subplot(1,3,1)
            F=  scatteredInterpolant(v,w,outL.(FreqBands{1}).A(:,AngularPlotIndices(1)));
            bF = F(bX,bY);
            imagesc(bF)
            colorbar
            title('IC of ph-frequency')
            subplot(1,3,2)
            F=  scatteredInterpolant(v,w,outH.(FreqBandsHigh{1}).A(:,AngularPlotIndices(2)));
            bF = F(bX,bY);
            imagesc(bF)
            colorbar
            title('IC of pow-frequency')
            subplot(1,3,3)
            polarplot(out2.phbins(:,1,1),sq(out2.pow_dens(:,1,1))'); axis tight; box off
            
        end
        
%         figure()
%         ax1 = subplot(3,1,1);
%         tRatio = SpecStates.RawSpecStates.t;
%         tICA = (1:length(out.(FreqBands{iFreq}).ICAsig(i,:)))/Fs;
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

end
