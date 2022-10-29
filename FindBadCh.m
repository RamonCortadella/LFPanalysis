function [BadChannels,RMSarray] = FindBadCh(Lfp, Fs, nRows, nCols, factor)
%this function takes finds the channels with too high and too low power and
%returns the array of indices corresponding to these "bad channels"

%- Lfp is an (1-dimensional - n x m) x time array
%- Fs (sampling frequency)
%- nRows and nCols are the number of rows (top to bottom) and columns (left
%  to right)

% take 10 random segments of 50 sec for
% nearest/next-nearest channels and central channel


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
dLfp = LfpFilt(2:end,:)-LfpFilt(1:end-1,:);

RMSarray = sqrt(mean((LfpFilt(Samps,:)-mean(LfpFilt(Samps,:))).^2));
BadChannels=find(RMSarray==0);
for i = 1:length(Lfp(1,:))
    
    lenZeros = find(Lfp(2:end,i)==0 & dLfp(:,i)==0);
    if lenZeros>=1
        BadChannels = [BadChannels,i];
    end
end
BadChannels=unique(sort(BadChannels));
RMSarray(find(RMSarray==0))=[];

for i=1:nCols*nRows %compute noise ratios for each channel
   
    frameA = repmat([-2*nRows -nRows 0 nRows 2*nRows],5,1);
    frameB = repmat([-2 -1 0 1 2]',1,5);
    frame = frameA + frameB;

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
%         display(f)
        for ii =f(1:end)
            if ii>= 1 & ii <= nCols*nRows
%                 display(ii)
                counter=counter+1;
                indices(counter)=ii;
            end
        end
    end
%     display(indices,'Indices')

% Compute correction ratio and correction sign
    RMScenter = sqrt(mean((LfpFilt(Samps,i)-mean(LfpFilt(Samps,i))).^2));
    RMSneighb = sqrt(mean((LfpFilt(Samps,indices)-mean(LfpFilt(Samps,indices))).^2));

%     RMSneighbSTD = std(log10(RMSneighb));
    if log10(RMScenter) >= mean(log10(RMSneighb))+log10(factor) | log10(RMScenter) <= mean(log10(RMSneighb))-log10(factor)
        BadChannels = [BadChannels i];
    end

end
BadChannels = unique(sort(BadChannels));