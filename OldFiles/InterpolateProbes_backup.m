function LfpInterp = InterpolateProbes(Lfp, badindex,nChRow,nChCol,depth, SingleShank)
%this function uses scatter interpolant to interpolate bad channels in an
%array. Lfp(:,1) contains all samples for a channel in position 1. Channels
%are sorted geometrically (i.e. 1,2,3,... along a column, and then next
%column). For a linear probe, a clone column is created for the scatter
%interpolant to work.

    LfpInterp = zeros([length(Lfp(:,1))-length(badindex), nChRow*nChCol]);


    if depth ==1 & SingleShank == 1 % a clone column is created for a linear probe
        nChCol = 2;
    end
%     if nChCol == 32
%         nChRow = 32;
%     end
    [bX, bY]= meshgrid([1:nChRow]',[1:nChCol]');

    if depth==0
        bX(:,floor(nChCol/2)+1:end) = bX(:,floor(nChCol/2)+1:end) + 10; %add spacing between Cols 1 to 8 and 9 to 16 to avoid interpolation across hemispheres   
    end
    if nChCol == 32
        bY(
    end
    
    v = bX;
    w = bY;
    if depth ==1 & SingleShank == 1
        v(:,badindex)=[];
        badindex(badindex==0)=[];
        Nbad = length(badindex);
        W(:,1:Nbad) = [];
    end
    v = reshape(v,[],1);
    w = reshape(w,[],1);
    if depth==0 | SingleShank==0
        v(badindex)=[];
        w(badindex)=[];
    end


    if depth ==1 & SingleShank == 1
        for k =1:length(Lfp(:,1))
            display(k,'sample number')

            map = Lfp(k,:);
            map(badindex)=[];%eliminate bad channels

            map = [map;map];
            map=reshape(map,[],1);

            F=  scatteredInterpolant(v,w,map');
            bF = F(bX,bY);
            LfpInterp(k,:) = reshape(bF,[],1);

        end
    else
        for k =1:length(Lfp(:,1))
            display(k,'sample number')
            map = Lfp(k,:);
            map(badindex)=[];%eliminate bad channels
            F=  scatteredInterpolant(v,w,map');
            bF = F(bX,bY);
            LfpInterp(k,:) = reshape(bF,[],1);
        end
    end