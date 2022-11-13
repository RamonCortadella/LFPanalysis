function LfpReduced = PCA_ICA2(FileName,T,indDB,OutputPath, varargin)
%function out = FunctionName(FileBase ,fMode, Arg1, Arg2)
%here is the help message :
[fMode, NumPCs, NumICs] = DefaultArgs(varargin,{'compute',100,50});
display(fMode)
display(NumPCs)
%body of the function
switch fMode
    case 'compute'
        display(['loading ' FileName])        
        
        nCh = T.NumCh(indDB);
        Lfp = LoadBinaryDAT(FileName, [0:nCh-1], nCh,1)';
        
        Fs = T.Fs(indDB);      
        nRows = T.nRows(indDB);
        nCols = T.nCols(indDB);

        [coeff, Data_PCA,~,~,~,~] = pca(Lfp,'NumComponents', NumPCs);

        %plot coefficients (loadings)
        [bX, bY]= meshgrid([1:nCols]',[1:nRows]');
        v = reshape(bX,[],1);
        w = reshape(bY,[],1);
        figure()
        for i = 1:NumPCs
            subplot(floor(sqrt(NumPCs))+1,floor(sqrt(NumPCs)),i)
            F=  scatteredInterpolant(v,w,coeff(:,i));
            bF = F(bX,bY);
            pcolor(bF)
            caxis([-0.3 0.3])
            colorbar
        end

        % Get Lfp with reduced dimensionality and apply ICA
        LfpReduced= Data_PCA*coeff';

        SaveBinary(strcat(FileName(1:end-4),'.PCA','.lfp'),LfpReduced);   
 
        
%         [ICAsig,A,W] = fastica(LfpReduced,'numOfIC',NumICs); %A = mixing matrix (X = A.S), W = demixing matrix (A^-1) (S = WX)
%         
%         for i = 1:NumICs
%             subplot(floor(sqrt(NumICs))+1,floor(sqrt(NumICs)),i)
%             F=  scatteredInterpolant(v,w,coeff(:,i));
%             bF = F(bX,bY);
%             pcolor(bF)
%             caxis([-0.3 0.3])
%             colorbar
%         end
        
        
    case 'display'
        load([FileBase '.' mfilename '.mat'],'out');
        %plot stuff ..
        
    case 'groupstats'
        % load and then compute stats and output 
end


