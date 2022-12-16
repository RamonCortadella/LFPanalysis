function out = PowerModulation(th,powL,pow,varargin)
%function out = PowerModulation(th,pow,Shuffle, UniformTrans)
%first compute the regular resultant legnth, MI and the significance 
% using parameters of Shuffle structure: Type (random, shift), nShuffle,
% MaxShift
[Shuffle, UniformTrans ] = DefaultArgs(varargin,{[],0});
N = length(th);
powsum = sum(pow);
powLsum = sum(powL);
if UniformTrans
    th = mud(th);
end
R = sum(exp(1i*th).*pow.*powL)./(powsum*powLsum);
out.Ramp = abs(R);
out.Rth = angle(R);


% now compute the MI via the enropy
nbins =18;
%th = mod(th+pi,2*pi)-pi; % to have it from -pi to pi always
phedges = linspace(-pi,pi,nbins+1); % 20 deg binsx
out.phbins = (phedges(1:end-1)+phedges(2:end))/2;
%get the density-like normalization of pow
[ph_cnt phind] = histcI(th,phedges);
pow_cnt = accumarray(phind,pow,[nbins 1],@sum);
pow_dens = pow_cnt./powsum;
Hpow = -sum(pow_dens.*log(pow_dens));
Hmax = log(nbins);
out.MI = (Hmax-Hpow)/Hmax;
out.pow_dens = pow_dens';
if ~isempty(Shuffle) & isstruct(Shuffle) 
    %then compute the significance
    %MaxShift = 125;
    
    Rsh = []; MIsh = [];
    for k=1:Shuffle.nShuffle
        %rand('state',sum(100*clock));
        switch Shuffle.Type
            case 'random'
                % via shuffling the pow values randomly
                pow_sh = randsample(pow,N);
            case 'shift'
                % via shifting by > 1 sec left/right
                shift = Shuffle.MaxShift+round((rand(1,1))*Shuffle.MaxShift);
                ind = mod([1:length(pow)]+shift,length(pow));
                ind(ind==0) =length(pow);
                pow_sh = pow(ind);
        end
        Rsh(k) = abs(sum(exp(1i*th).*pow_sh))./powsum;

        pow_cntsh = accumarray(phind,pow_sh,[nbins 1],@sum);
        pow_denssh = pow_cntsh./powsum;
        Hpowsh = -sum(pow_denssh.*log(pow_denssh));
        MIsh(k) = (Hmax-Hpowsh)/Hmax;

    end
    %stat_pm = CatStruct(stat_pm);
    if 0
        out.stat(shtype).Rpval = sum(Rsh>out.Ramp)/Shuffle.nShuffle;
        out.stat(shtype).MIpval = sum(MIsh>out.MI)/Shuffle.nShuffle;
        out.stat(shtype).Rsh_mean = mean(Rsh);
        out.stat(shtype).Rsh_std = std(Rsh);

        out.stat(shtype).MIsh_mean = mean(MIsh);
        out.stat(shtype).MIsh_std = std(MIsh);

        out.stat(shtype).Rzscore = (out.Ramp-out.stat(shtype).Rsh_mean)/out.stat(shtype).Rsh_std;
        out.stat(shtype).MIzscore = (out.MI-out.stat(shtype).MIsh_mean)/out.stat(shtype).MIsh_std;
    else
        out.Rpval = sum(Rsh>out.Ramp)/Shuffle.nShuffle;
        out.MIpval = sum(MIsh>out.MI)/Shuffle.nShuffle;
        Rsh_mean = mean(Rsh);
        Rsh_std = std(Rsh);

        MIsh_mean = mean(MIsh);
        MIsh_std = std(MIsh);

        out.Rzscore = (out.Ramp-Rsh_mean)/Rsh_std;
        out.MIzscore = (out.MI-MIsh_mean)/MIsh_std;

    end


end