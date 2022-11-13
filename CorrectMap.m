function CorrectMap(FileName,T,indDB)
    nCh = T.NumCh(indDB);
    nRows = T.nRows(indDB);
    nCols = T.nCols(indDB);
    
    if nCh~= 256 | nRows~=16 | nCols~=16
        display('Wrong file!!!!')
    else
        
        Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
        oddind = [1:2:256];
        evenind = [2:2:256];
        LfpCorr(:,oddind) = Lfp(:,evenind);
        LfpCorr(:,evenind) = Lfp(:,oddind);

        SaveBinary([FileName(1:end-4) '.lfp'], LfpCorr);
        display('saved corrected file')
    end