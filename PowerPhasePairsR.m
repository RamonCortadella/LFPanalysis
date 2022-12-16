function [out PhMat PowMat] = PowerPhasePairs(xPh, FrBinsPh, xPow,FrBinsPow,ResampleNum,FilterType, varargin)
%function [out PhMat PowMat] = PowerPhasePairs(xPh, FrBinsPh, xPow,FrBinsPow,Resample,FilterType, FunHandle, Args)
% FrBinsXXX - edges of frequency bins for phase/power calculation in norm.
% units Hz/F_Nuiqwist
% xPh is 1d vector from which the phase will be computed
% xPow could be a matrix (many channels)
%   or structure with FileName, nChannels,
UniPhase = 0;
nBinsPow = size(FrBinsPow,1);
nBinsPh = size(FrBinsPh,1);
n = size(xPh,1);
if isempty(FilterType)
    filttype    = 'but';
else
    filttype = FilterType;
end
if size(xPh,1)~=n error('xPow and xPh shoudl be of same size'); end
%if size(xPh,2)>1  error('so far works with 1d matrix of phase - one reference phase'); end
nPhChan = size(xPh,2);
%if any(FrBinsPow>0.5) |any(FrBinsPh>0.5)
%    error('phases have to be from 1 ot 1/2');
%end

% if nargin>6 & ~isa(varargin{1},'function_handle')
%     error('there should be a function handle argument');
% end

% allocate memory
nr=size(resample(rand(n,1),1,ResampleNum),1);
if nBinsPh>1
    PowMat = zeros(nr, nBinsPow);
end
PhMat = zeros(nr, nPhChan, nBinsPh);
if FrBinsPh==0
    ph = unwrap(xPh);
    if ResampleNum>1
        ph = ph(1:ResampleNum:end);
        ph = ph(1:nr);
    end
    PhMat = mod(ph+pi,2*pi)-pi;
else
    
    %ER: -------Original code which was moved outside the loop ------------------------
    %disp('A bug was fixed here by Evgeny on 17.10.12!!')
    if ResampleNum>1
        xPh = resample(xPh,1,ResampleNum);
    end
    %-----------------------------------------------------------------------

    for kPh = 1:nBinsPh
        if 0
            ph = angle(hilbert(ButFilter(xPh, 4, FrBinsPh(kPh,:), 'bandpass')));
            if ResampleNum>1

                ph = unwrap(ph);
                ph = resample(ph,1,ResampleNum);
                ph = mod(ph+pi,2*pi)-pi;
            end
        else
            %ER: -------Original code which was moved outside the loop, because xPh must be downsampled only ones! ------------------------
%             if ResampleNum>1
%                 xPh = resample(xPh,1,ResampleNum);
%             end
            %---------------------------------------------
%             display(FrBinsPh(kPh,:)*ResampleNum)
%             display(FrBinsPh(kPh,:))
            ph = angle(hilbert(ButFilter(xPh, 2, FrBinsPh(kPh,:)*ResampleNum, 'bandpass')));
            if UniPhase
                for kPhCh=1:nPhChan
                    ph(:,kPhCh) =MakeUniformDistr(ph(:,kPhCh),-pi, pi);
                end
            end

        end
        if size(ph,2)==1
            PhMat(:,1,kPh) = ph;
        else
            PhMat(:,:,kPh) = ph;
        end
    end
end
% if nargin<5
%     %just return the power and phase matrices
%     out = 'no_function_handle';
%     if nPow>1 warning('returning only one column power values for the sake of memory');end
% end
if ~isstruct(xPow)
    nPow = size(xPow,2);
else
    nPow = xPow.nPow;
end

out = struct();
for knPow = 1:nPow
    if ~isstruct(xPow)
        myxPow = xPow(:,knPow);
    else
        myxPow = xPow.FunHandle(knPow);
    end
    if nBinsPh>1
        for kPow=1:nBinsPow
            PowMat(:,kPow) = FilterPower(kPow);
        end
        for kPhCh = 1:nPhChan
            for kPh = 1:nBinsPh
                for kPow=1:nBinsPow
                    pm = feval(varargin{1},PhMat(:,kPhCh,kPh),PowMat(:,kPow),varargin{2:end});
                    out = SetStructFields(out,pm,{kPh,kPow,knPow,kPhCh});
                end
            end
        end
    else
        tic;
        for kPow=1:nBinsPow
            tic;
            pow = FilterPower(kPow);
            pm = feval(varargin{1},PhMat,pow,varargin{2:end});
            out = SetStructFields(out,pm,{kPow,knPow});
           % fprintf('Freq. bin %d  took %f sec\n',kPow,toc);
        end
        %fprintf('Channel %d  took %f sec\n',knPow,toc);
    end
    %keyboard
end
out.fPh = mean(FrBinsPh,2);
out.fPow = mean(FrBinsPow,2);

    function pow = FilterPower(kPow)

        switch filttype
            case 'but'
                feeg = ButFilter(myxPow, 2, FrBinsPow(kPow,:), 'bandpass');
            case 'fir'
                feeg = FirFilter(myxPow, 0, FrBinsPow(kPow,:), 'bandpass');
            case 'mt'
                feeg = MTFilter(myxPow, FrBinsPow(kPow,:)*Fs/2, Fs,0);
        end
        pow = abs(hilbert(feeg));
        if ResampleNum>1
            pow = resample(pow,1,ResampleNum);
        end

    end
end
