close('all')

Directories = ["../../../data/ASIC512/B13907W21-T1-rat01601/DatData/SIG_B13907W21-T1_InVivo-13_F-AC";];
% a(0)  
for iD = 1:length(Directories)
    Directory = Directories(iD);
%     Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/DatData/ClippedMapped/B13289O14-DH1-Rec3' 'AC.xml']);
    Par =  LoadXml(['../../../data/ASIC512/B13907W21-T1-rat01601/DatData/SIG_B13907W21-T1_InVivo-13_F-AC' '.xml']);

    d = dir(strcat(Directory, '*.dat'));%DC-LP30Hz-Notch50-100Hz.dat');

    LFPFs = 500;
    DownFact = 100;
    Tstab = 0;
    NCh = 512;
    nChRow = 16;
    nChCol = 32;
    
    %% get ephys data
    ACLfp = [];
    DCLfp = [];
    ClipLfp = [];

    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:NCh-1], Par.nChannels,1)';

        SplitName = split(d(fn).name,'-');
        
        if SplitName{4}(end-5:end-4)=='AC'
    %     myrLfp = resample(Lfp,1,10); % RESAMPLE TO 100 HZ !!!
            ACLfp = Lfp;
%         elseif  SplitName{3}(end-5:end-4)=='DC' 
%             DCLfp = Lfp;
%         elseif SplitName{3}(end-7:end-4)=='Clip'
%             ClipLfp = Lfp;
        end
    end

    %% rearrange and interpolate spatiall
    LfpGeom = zeros([size(ACLfp)]);
%     LfpGeomDC = zeros([size(DCLfp)]);
    LfpGeomClip2 =  zeros([size(ACLfp)]);

    for i = 1:length(Par.AnatGrps)
        for ii = 1:length(Par.AnatGrps(i).Channels)
            LfpGeom(:,ii+(i-1)*nChCol) = ACLfp(:,Par.AnatGrps(i).Channels(ii)+1);
%             LfpGeomDC(:,ii+(i-1)*16) = DCLfp(:,Par.AnatGrps(i).Channels(ii)+1);
            LfpGeomClip2(:,ii+(i-1)*nChCol) = logical(ones([length(LfpGeomClip2(:,1)),1])*Par.AnatGrps(i).Skip(ii));
            display(Par.AnatGrps(i).Skip(ii))
        end
    end
    

    figure()

    LfpStruct = struct('AC', LfpGeom);
    fn = fieldnames(LfpStruct);

    for k=1:numel(fn)
        Lfp = LfpStruct.(fn{k});
        if fn{k} == 'DC'
            LfpGeomClip= LfpGeomClip2(1:DownFact:end,:);
            LFPfs = LFPFs/DownFact;
        else
            LfpGeomClip= LfpGeomClip2;
            LFPfs = LFPFs;
        end
        LfpClipInter = [];
        LfpClipInterSum = [];
        
        for i = 1:length(Lfp(1,:))
            ind0 = find(LfpGeomClip(:,i)==1);
            if length(ind0)>= 1
                ind1 = find(ind0 >= floor(Tstab*LFPfs));
                
                if length(ind1)>= 1
                    % fill index of clip for transitions of less than 1 s long
                    deltaInd = ind1(2:end)-ind1(1:end-1);
                    ind2 = find(deltaInd<=LFPfs & deltaInd>= 2);
                    ind3 = ind0;
                    
                    if length(ind2)>=1
                        for ii = ind2
                            ind3 = [ind3 linspace(ind1(ii),ind1(ii)+floor(LFPfs),floor(LFPfs)+1)];
                            ind3 = unique(ind3);
                        end
                    end
                    LfpGeomClip(ind3,i) = 1;
                end
            end
        end
        
        
        
        
        for i =1: length(Lfp(1,:))  
            display(i)

            %add only if the indices are valid
            if (mod(i,nChCol) >= 2) 
                mask1 = ((1-LfpGeomClip(:,i-1)).*(1-LfpGeomClip(:,i+1)));%this is the mask where the two indices are valid
                display(max(LfpGeomClip(:,i-1)))
                display(max(LfpGeomClip(:,i+1)))
                ind0 = find(mask1 == 0);
                if length(ind0)>= 1
                    ind1 = find(ind0 >= floor(Tstab*LFPfs));
                    
                    if length(ind1)>= 1
                        % fill index of clip for transitions of less than 1 s long
                        deltaInd = ind1(2:end)-ind1(1:end-1);
                        ind2 = find(deltaInd<=LFPfs & deltaInd>= 2);
                        ind3 = ind0;
                        
                        if length(ind2)>=1
                            for ii = ind2
                                ind3 = [ind3 linspace(ind1(ii),ind1(ii)+floor(LFPfs),floor(LFPfs)+1)];
                                ind3 = unique(ind3);
                            end
                        end
                        mask1(ind3) = 0;
                    end
                end
                int1 = mask1.*(Lfp(:,i-1)+Lfp(:,i+1))/2; %this is the mean where the two indices are valid
            else
                display('in')
%                 display(min(mask1))
                mask1 = LfpGeomClip(:,19)*0;
                int1 = mask1;
            end
            if (floor(i/nChCol) >=1) & (floor(i/nChCol) <= nChRow-2) & (mod(i,nChCol)~=0) 
                mask2 = ((1-LfpGeomClip(:,i-nChCol)).*(1-LfpGeomClip(:,i+nChCol)));
                ind0 = find(mask2 == 0);
                if length(ind0)>= 1
                    ind1 = find(ind0 >= floor(Tstab*LFPfs));
                    if length(ind1)>= 1
                        % fill index of clip for transitions of less than 1 s long
                        deltaInd = ind1(2:end)-ind1(1:end-1);
                        ind2 = find(deltaInd<=LFPfs & deltaInd>= 2);
                        ind3 = ind0;
                        
                        if length(ind2)>=1
                            for ii = ind2
                                ind3 = [ind3 linspace(ind1(ii),ind1(ii)+floor(LFPfs),floor(LFPfs)+1)];
                                ind3 = unique(ind3);
                            end
                        end
                        mask2(ind3) = 0;
                    end
                end
                int2 = mask2.*(Lfp(:,i-nChCol)+Lfp(:,i+nChCol))/2;
            else
                display('in')
%                 display(min(mask2))
                mask2 = LfpGeomClip(:,19)*0;
                int2 = mask2;
            end
            if (floor(i/nChCol) >=1) & (floor(i/nChCol) <= nChRow-2) & (mod(i,nChCol) >= 2) 
                mask3 = ((1-LfpGeomClip(:,i-nChCol-1)).*(1-LfpGeomClip(:,i+nChCol+1)));
                ind0 = find(mask3 == 0);
                if length(ind0)>= 1
                    ind1 = find(ind0 >= floor(Tstab*LFPfs));
                    if length(ind1)>= 1
                        % fill index of clip for transitions of less than 1 s long
                        deltaInd = ind1(2:end)-ind1(1:end-1);
                        ind2 = find(deltaInd<=LFPfs & deltaInd>= 2);
                        ind3 = ind0;
                        
                        if length(ind2)>=1
                            for ii = ind2
                                ind3 = [ind3 linspace(ind1(ii),ind1(ii)+floor(LFPfs),floor(LFPfs)+1)];
                                ind3 = unique(ind3);
                            end
                        end
                        mask3(ind3) = 0;
                    end
                end
                int3 = mask3.*(Lfp(:,i-nChCol-1)+Lfp(:,i+nChCol+1))/2;
            else
                display('in')
%                 display(min(mask3))
                mask3 = LfpGeomClip(:,19)*0;
                int3 = mask3;
            end
            if (floor(i/nChCol) >=1) & (floor(i/nChCol) <= nChRow-2) & (mod(i,nChCol) >= 2)
                mask4 = ((1-LfpGeomClip(:,i+nChCol-1)).*(1-LfpGeomClip(:,i-nChCol+1)));
                ind0 = find(mask4 == 0);
                if length(ind0)>= 1
                    ind1 = find(ind0 >= floor(Tstab*LFPfs));
                    if length(ind1)>= 1
                        % fill index of clip for transitions of less than 1 s long
                        deltaInd = ind1(2:end)-ind1(1:end-1);
                        ind2 = find(deltaInd<=LFPfs & deltaInd>= 2);
                        ind3 = ind0;
                        if length(ind2)>=1
                            for ii = ind2
                                ind3 = [ind3 linspace(ind1(ii),ind1(ii)+floor(LFPfs),floor(LFPfs)+1)];
                                ind3 = unique(ind3);
                            end
                        end
                        mask4(ind3) = 0;
                    end
                end
                int4 = mask4.*(Lfp(:,i+nChCol-1)+Lfp(:,i-nChCol+1))/2;
            else
                display('in')
%                 display(min(mask4))
                mask4 = LfpGeomClip(:,19)*0;
                int4 = mask4;
            end
            Interp = (int1+int2+int3+int4);
            MaskSum = (mask1+mask2+mask3+mask4);
            MaskT = mask1|mask2|mask3|mask4;
            index0 = find(MaskSum == 0);
            if length(index0)>= 1
                MaskSum(index0) = 1000000;
            end


            Interp = LfpGeomClip(:,i).* Interp .* MaskT ./ MaskSum;

            display(min(Lfp(:,i)))
            display(max(Lfp(:,i)))
            display(max(abs(MaskSum)))

            Lfp(:,i) = (1-LfpGeomClip(:,i)).*Lfp(:,i) + Interp;
            LfpClipInter(:,i) = LfpGeomClip(:,i) & (1-mask1) & (1-mask2) & (1-mask3) & (1-mask4);
        end
        LfpClipInterSum = sum(LfpClipInter,2);

        plot(LfpClipInterSum)
        LfpStruct.(fn{k}) = Lfp;

    end

    %%
    SaveBinary(strcat(Directory,'interAC.lfp'), LfpStruct.AC)
%     SaveBinary(strcat(Directory,'interDC.lfp'), LfpStruct.DC)
%     SaveBinary(strcat(Directory,'interClip.lfp'), LfpClipInter)
end