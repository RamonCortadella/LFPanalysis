function PCA_ICA(Lfp, Fs, )
q = 30;


Fs = 651.04166667;
DownSampDC = 20;
FsDC = Fs/DownSampDC;

SpeedDCFact=100;
SpeedACFact = 1;
OrderFiltAC = 4;
OrderFiltDC = 2;

FminAC = 2;
FmaxAC = 10;
FminDC = 0.05;
FmaxDC = 0.2;
SlowDown = 30;
CnstScale = false;
GaussianFact= 0.5;

[coeff, Data_PCA,latent,tsquared,explained,mu] = pca(LfpGeom, 'NumComponents', q);


[bX, bY]= meshgrid([1:nCols]',[1:nRows]');
v = reshape(bX,[],1);
w = reshape(bY,[],1);
for i = 1:q
    subplot(floor(sqrt(q))+1,floor(sqrt(q)),i)
    F=  scatteredInterpolant(v,w,coeff(:,i));
    bF = F(bX,bY);
    pcolor(bF)
    caxis([-0.3 0.3])
    colorbar
end
%%
ArtifactPCs= [2,5];

Data_PCA_Artifact = Data_PCA(:,ArtifactPCs);
Data_PCA_ArtifactDC = LfpGeomDC*coeff(:,ArtifactPCs);

counter = 0;
LfpNoArtifact = LfpGeom;
LfpNoArtifactDC = LfpGeomDC;
for i = ArtifactPCs
    counter = counter +1;
    LfpNoArtifact = LfpNoArtifact-Data_PCA_Artifact(:,counter)*coeff(:,i)';
    LfpNoArtifactDC = LfpNoArtifactDC-Data_PCA_ArtifactDC(:,counter)*coeff(:,i)';
end
