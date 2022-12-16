function LfpInterp = InterpolateProbes2(FileName, T,indDB,OutputPath)
%this function uses scatter interpolant to interpolate bad channels in an
%array. Lfp(:,1) contains all samples for a channel in position 1. Channels
%are sorted geometrically (i.e. 1,2,3,... along a column, and then next
%column). For a linear probe, a clone column is created for the scatter
%interpolant to work.

% [fMode] = DefaultArgs(varargin,{'compute'});
    
    nCh = T.NumCh(indDB);
    Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
    display('Loading Lfp')
    display('Lfp loaded')
    display(size(Lfp),'sizeLFP')

    badindex = str2num(T.badChannels{indDB});
    nChRow = T.nRows(indDB);
    nChCol = T.nCols(indDB);
    depth = T.depth(indDB);
    SingleShank = T.SingleShank(indDB);

    if depth ==1 & SingleShank == 1 % a clone column is created for a linear probe
        nChCol = 2;
    end
    
%     if nChCol == 32
%         nChRow = 32;
%     end
    [bX, bY]= meshgrid([1:nChCol]',[1:nChRow]');
    
    if depth==0
        bX(:,floor(nChCol/2)+1:end) = bX(:,floor(nChCol/2)+1:end) + 11; %add spacing between Cols 1 to 8 and 9 to 16 to avoid interpolation across hemispheres   
    end
    %adjust geometrically the 32x16 arrays
    if nChCol == 32
        for i=1:floor(nChRow/2)
            if rem(i,2)==0
                bY(:,i) = bY(:,i)*2-1;
                bY(:,i+16) = bY(:,i)*2;
            else
                bY(:,i) = bY(:,i)*2;
                bY(:,i+16) = bY(:,i)*2-1;
            end
        end
    end
    %adjust geometrically the 32x8 depth arrays
    if depth == 1 & nChCol == 8
        for i = 1:nChCol
            bX(:,i)= (i-1)*80+1;
        end
        for i=1:nChRow
            bY(i,:) = (i-1)*13+1; 
        end
    end
    
    v = bX;
    w = bY;

    if depth ==1 & SingleShank == 1
        v(badindex,:)=[];
        badindex(badindex==0)=[];
        Nbad = length(badindex);
        w(1:Nbad,:) = [];
    end
    v = reshape(v,[],1);
    w = reshape(w,[],1);
    if SingleShank==0
        v(badindex)=[];
        w(badindex)=[];
    end

    if nChCol == 32
        [bX, bY]= meshgrid([1:nChCol]',[1:32]');
        LfpInterp = zeros([length(Lfp(:,1)), 32*nChCol]);
    else
        LfpInterp = zeros([length(Lfp(:,1)), nChRow*nChCol]);  
    end

    if depth ==1 & SingleShank == 1
        [bX, bY]= meshgrid([1:1]',[1:nChRow]');
        Lfp = [Lfp,Lfp];
        badindex = [badindex,badindex+32];
        LfpInterp = zeros([length(Lfp(:,1)), nChRow*1]);  
    end
%     display(size(Lfp),'sizeLfp')
    for k =1:length(Lfp(:,1))
%         display(k,'sample number')
        map = Lfp(k,:);
        map(badindex)=[];%eliminate bad channels
        F=  scatteredInterpolant(v,w,map');
        bF = F(bX,bY);
        LfpInterp(k,:) = reshape(bF,[],1);
    end
    
    SaveBinary(strcat(FileName(1:end-4),'.interp.lfp'), LfpInterp);
end