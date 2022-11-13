function T = CorrectBadChannels(T,indDB)
    nCh = T.NumCh(indDB);
    nRows = T.nRows(indDB);
    nCols = T.nCols(indDB);
    
    BadChannels = str2num(T.badChannels{indDB});
    
    if nCh~= 256 | nRows~=16 | nCols~=16
        display('Wrong file!!!!')
    else
      
        BadChannelsCorr(mod(BadChannels,2)==0) = BadChannels(mod(BadChannels,2)==0)-1;
        BadChannelsCorr(mod(BadChannels,2)~=0) = BadChannels(mod(BadChannels,2)~=0)+1;
        
        
        T.badChannels(indDB)= {mat2str(BadChannelsCorr)};
        display('returned Bad Channels corrected')
    end