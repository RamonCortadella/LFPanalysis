function LfpCorr = GmSmoothening(Lfp, band, Fs, nRows, nCols, Chs, iterations, BadChannels)
%this function takes an array of signals (1-dimensional - n x m) x time and reduces the
%gain noise by calculating the ratio between the RMS of nearest/next-nearest
%neighbours and the central channel. Possible phase inversion (i.e. wrong 
%sign of gain is detected by calculating signal diferences and sums and
%comparing their RMS value

%- Lfp is an (1-dimensional - n x m) x time array
%- band defines the frequency in Hz of the band taken to calculate gain
%  differences (e.g. [1,10])
%- Fs (sampling frequency)
%- nRows and nCols are the number of rows (top to bottom) and columns (left
%  to right)
%- Chs is an array of channel indices taken into account (e.g. [1:256], or [257:288])
%- Iterations defines the number of smoothenings of Gm

% take 20 random segments of 5 sec (exclude last 10 s) for
% nearest/next-nearest channels and central channel


SegmentSize = floor(length(Lfp(:,1))/21);%in units of samples
SegmentsInit = [1:20]'*SegmentSize;%in units of samples
Samps=[];
LfpCorr = zeros(size(Lfp));

for i = 1:20
    Samps = [Samps , SegmentsInit(i): SegmentsInit(i)+Fs*5];
end
%generate coordinates arrays M1 and M2
M1 = repmat([1:nRows]',1,nCols);
M2 = repmat([1:nCols],nRows,1);

for it = 1:iterations
    LfpFilt = ButFilter(Lfp,4,band/(Fs/2),'bandpass');
    for i=Chs %smoothen gain for ith channel
        display(i,'Channel')
        frame = [i-nRows-1, i-nRows, i-nRows+1, i-1, i+1, i+nRows-1, i+nRows, i+nRows+1];
        indices=0;
        counter=0;
        for ii=frame
            if ii >= 1 & ii <= length(Chs)
                if abs(M1(ii)-M1(i))<=1 & abs(M2(ii)-M2(i))<=1
                    if length(find(BadChannels==ii))==0
                        counter=counter+1;
                        coords = vertcat([M1(ii),M2(ii)]);
                        indices(counter) = M1(ii)+(M2(ii)-1)*nRows; %select only neighbour/neares-neighbour channels (edge effects taken into account)
                    end
                end
            end
        end

        if counter==0
            CorrRatio = 0;
            CorrSign = 1;
        else
            RMSneighb = sqrt(mean((LfpFilt(Samps,indices)-mean(LfpFilt(Samps,indices))).^2));
            Diff = sqrt(mean(((LfpFilt(Samps,indices)-mean(LfpFilt(Samps,indices)))-(LfpFilt(Samps,i)-mean(LfpFilt(Samps,indices)))).^2)) - sqrt(mean(((LfpFilt(Samps,indices)-mean(LfpFilt(Samps,indices)))+(LfpFilt(Samps,i)-mean(LfpFilt(Samps,indices)))).^2));      
    % Compute correction ratio and correction sign
            RMScenter = sqrt(mean((LfpFilt(Samps,i)-mean(LfpFilt(Samps,i))).^2));
            CorrRatio = median(RMSneighb/RMScenter);
            CorrSign = (Diff < 0);
            CorrSign(find(CorrSign==0)) = -1;
            CorrSign = median(CorrSign);
        end
        LfpCorr(:,i) = Lfp(:,i)*CorrSign*CorrRatio;
    end
    Lfp = LfpCorr;
end

    
    
    
    
    