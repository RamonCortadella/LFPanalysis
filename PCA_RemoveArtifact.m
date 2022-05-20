directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/';
Par =  LoadXml(strcat(directory,'DatData/ClippedMapped/B13289O14-DH1-Rec3interAC.xml')); %Rec1, Rec2, Rec3,Rec4, Rec5, Rec6, Rec7,Rec8
videoOutPath = strcat(directory,'MatlabData/Videos/delta-ISA-NoArifact-gaussian0p5.avi');
videoOutPathArtif = strcat(directory,'MatlabData/Videos/delta-ISA-Arifact-gaussian0p5.avi');
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

for i = [3]
    d = dir(strcat(directory,'DatData/ClippedMapped/B13289O14-DH1-Rec',int2str(i),'inter*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');


    ACLfp = [];
    DCLfp = [];
    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:255], Par.nChannels,1)';

        SplitName = split(d(fn).name,'-');

        if SplitName{3}(end-5:end-4)=='AC'
            LfpGeom = Lfp;
        end
        if SplitName{3}(end-5:end-4)=='DC'
            LfpGeomDC = Lfp;
        end
    end
end

[coeff, Data_PCA,latent,tsquared,explained,mu] = pca(LfpGeom, 'NumComponents', q);

[bX, bY]= meshgrid([1:16]',[1:16]');
Coord(:,1) = reshape(repmat([1:16],16,1),[],1);
Coord(:,2)= repmat([1:16]',16,1);


close('all')
for i = 1:30
    subplot(6,6,i)
    map = coeff(:,i);
    F=  scatteredInterpolant(Coord(:,1),Coord(:,2),map);
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
ECoG2Video(LfpNoArtifact(floor(length(LfpNoArtifact(:,1))/10):2*floor(length(LfpNoArtifact(:,1))/10),:), LfpNoArtifactDC(1:floor(length(LfpNoArtifactDC(:,1))),:), videoOutPath, Fs, DownSampDC, SpeedDCFact,SpeedACFact, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,CnstScale, GaussianFact)
ECoG2Video(LfpGeom(floor(length(LfpGeom(:,1))/10):2*floor(length(LfpGeom(:,1))/10),:), LfpGeomDC(1:floor(length(LfpGeomDC(:,1))),:), videoOutPathArtif, Fs, DownSampDC, SpeedDCFact,SpeedACFact, FminAC, FmaxAC, FminDC, FmaxDC, OrderFiltAC, OrderFiltDC, SlowDown,CnstScale, GaussianFact)
