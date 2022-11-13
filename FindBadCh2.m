function out = FindBadCh2(FileName , InputPath, T, IndDB, OutputPathTable,FileNameDB, varargin)
%this function takes finds the channels with too high and too low power and
%returns the array of indices corresponding to these "bad channels"

%- Lfp is an (1-dimensional - n x m) x time array
%- Fs (sampling frequency)
%- nRows and nCols are the number of rows (top to bottom) and columns (left
%  to right)

% take 10 random segments of 50 sec for
% nearest/next-nearest channels and central channel

[fMode, factor, interactive] = DefaultArgs(varargin,{'compute',2.5,'off'});

%body of the function
switch fMode
    case 'compute'
        display(['loading ' FileName])
        %if you give output from the function make it structure :
        nCh = T.NumCh(IndDB);
        Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
        %load params
        nRows = T.nRows(IndDB);
        nCols = T.nCols(IndDB);
        Chs = [1:nRows*nCols];
        Fs = T.Fs(IndDB);

        SegmentSize = floor(length(Lfp(:,1))/11);%in units of samples
        SegmentsInit = [1:10]'*SegmentSize;%in units of samples
        Samps=[];
        % LfpCorr = zeros(size(Lfp));
        if length(Lfp(:,1))/Fs <= 600
            display('¡¡¡¡¡¡ RECORDING TOO SHORT TO COMPUTE RMS RELIABLY !!!!!')
            return
        end

        for i = 1:10
            Samps = [Samps , SegmentsInit(i): SegmentsInit(i)+Fs*50];
        end
        %generate coordinates arrays M1 and M2
        M1 = repmat([1:nRows]',1,nCols);
        M2 = repmat([1:nCols],nRows,1);

        LfpFilt = ButFilter(Lfp,4,[1, 10]/(Fs/2),'bandpass');
        if T.CoupledAC(IndDB)
            LfpFiltH = ButFilter(Lfp,4,[60, 120]/(Fs/2),'bandpass');
        end
        dLfp = Lfp(2:end,:)-Lfp(1:end-1,:);

        RMSarray = sqrt(mean((LfpFilt(Samps,:)-mean(LfpFilt(Samps,:))).^2));
        counter = -1;
        while true
            if strcmp(interactive, 'on')
                loop = input('(re)-compute? (1-yes/0-no)');
            else
                loop=1;
                counter = counter+1;
            end
            if loop ==1 & counter <= 0
                display('continued')
                
                BadChannels=find(RMSarray==0);
                for i = 1:length(Lfp(1,:))

                    lenZeros = length(find(Lfp(2:end,i)==0 & dLfp(:,i)==0));
                    if lenZeros>=1000000
                        BadChannels = [BadChannels,i];
                    end
                end
                BadChannels=unique(sort(BadChannels));
        %         RMSarray(find(RMSarray==0))=[];

                for i=1:nCols*nRows %compute noise ratios for each channel

                    frameA = repmat([-2*nRows -nRows 0 nRows 2*nRows],5,1);
                    frameB = repmat([-2 -1 0 1 2]',1,5);
                    frame = frameA + frameB;
                    display(i,'channel')

                    indices=0;
                    counter=0;
                    if nCols>1
                        f = frame+i;
                        for ii= f(1:end)
                            if ii >= 1 & ii <= nCols*nRows
                                if abs(M1(ii)-M1(i))<=2 & abs(M2(ii)-M2(i))<=2
                                    if length(find(BadChannels==ii))==0
                                        counter=counter+1;
                                        indices(counter) = M1(ii)+(M2(ii)-1)*nRows; %select only neighbour/neares-neighbour channels (edge effects taken into account)
                                    end
                                end
                            end
                        end
                    else
                        f = [-2,-1,1,2]+i;

                        for ii =f(1:end)
                            if ii>= 1 & ii <= nCols*nRows
                                counter=counter+1;
                                indices(counter)=ii;
                            end
                        end
                    end

                % Compute correction ratio and correction sign
                    if T.CoupledAC(IndDB)
         
                        RMScenter = sqrt(mean((LfpFilt(Samps,i)-mean(LfpFilt(Samps,i))).^2));
                        RMScenterH = sqrt(mean((LfpFiltH(Samps,i)-mean(LfpFiltH(Samps,i))).^2));

                        if indices == 0
                            RMSneighb = 1e-20 ; 
                            RMSneighbH = 1e-20;
                        else
                            RMSneighb = sqrt(mean((LfpFilt(Samps,indices)-mean(LfpFilt(Samps,indices))).^2));
                            RMSneighbH = sqrt(mean((LfpFiltH(Samps,indices)-mean(LfpFiltH(Samps,indices))).^2));

                        end


                        if log10(RMScenter) >= median(log10(RMSneighb))+log10(factor) | log10(RMScenter) <= median(log10(RMSneighb))-log10(factor)
                            BadChannels = [BadChannels i];
                        elseif log10(RMScenterH) >= median(log10(RMSneighbH))+log10(factor) | log10(RMScenterH) <= median(log10(RMSneighbH))-log10(factor)
                            BadChannels = [BadChannels i];
                        end
                    else
                        RMScenter = sqrt(mean((LfpFilt(Samps,i)-mean(LfpFilt(Samps,i))).^2));
                        if indices == 0
                            RMSneighb = 1e-20 ; 
                        else
                            RMSneighb = sqrt(mean((LfpFilt(Samps,indices)-mean(LfpFilt(Samps,indices))).^2));
                        end

                        if log10(RMScenter) >= median(log10(RMSneighb))+log10(factor) | log10(RMScenter) <= median(log10(RMSneighb))-log10(factor)
                            BadChannels = [BadChannels i];
                        end
                    end
                end
                BadChannels = unique(sort(BadChannels));

                out.BadChannels = BadChannels;
                out.RMSarray =RMSarray;
                out.LowThr = median(log10(RMSneighb/factor));
                out.HighThr = median(log10(RMSneighb*factor));
                
                T.badChannels(IndDB)= {mat2str(BadChannels)};
                writetable(T,strcat(OutputPathTable,FileNameDB,'.xlsx'))
                
                Mask = zeros([nRows*nCols,1]);
                Mask(BadChannels)=1;
                Mask = reshape(Mask, nRows,nCols);
                
                figure()
                imshow(Mask,[],'InitialMagnification',10000)
                colormap('parula')
                set(gca,'Units','Normalized','OuterPosition',[0 0 1 1])
                if strcmp(interactive, 'on')
                    deltaThreshold = input('Change in threshold for bad channel detection (type number and/or press INTRO)');
                else 
                    deltaThreshold=[];
                end
                if deltaThreshold
                    factor = factor + deltaThreshold;
                end
                display(factor,'current factor')
                
                
        %         save([FileBase '.' mfilename '.mat'],'out');
            else
%                 path = strcat(InputPath,'/States/',T.RecordingId{IndDB},'-BadChannels.mat');
%                 save(path,'out');
                T.FactorBadChannels(IndDB) = factor;
                out.T =T;
%                 writetable(T,strcat(OutputPathTable,FileNameDB,'.xlsx'))
                display('jumped to next')
                break
            end
        end
    case 'display'
        load(strcat(InputPath,'/States/',T.RecordingId{IndDB},'-BadChannels.mat'));
        figure()
        hold on
        nbins=round((max(log10(out.RMSarray))-min(log10(out.RMSarray)))/0.02);
        display(nbins)
        histogram(log10(out.RMSarray),nbins)

        histogram(log10(out.RMSarray(out.BadChannels)),nbins)
%         xline(median(log10(out.RMSarray)),'r')
        xlim([2.2 3.8])
        
    case 'groupstats'
        % load and then compute stats and output 
end






