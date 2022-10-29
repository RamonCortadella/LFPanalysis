close('all')
close('all')

Directories = ["../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec1";
               "../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec2";
               "../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec3";
               "../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec4";
               "../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec11";
               "../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec12";];
%                "../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec7";
%                "../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec8";
%                "../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec9";
%                "../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec10";
%                "../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec11";
%                "../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec13"];
            
for iD = 1:length(Directories)
    Directory = Directories(iD);
    
    if iD==1
        Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec1' 'AC.xml']);
    end
    if iD==2
        Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec2' 'AC.xml']);
    end
    if iD==3
        Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec3' 'AC.xml']);
    end
    if iD==4
        Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec4' 'AC.xml']);
    end
    if iD==5
        Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec11' 'AC.xml']);
    end
    if iD==6
        Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH2-01551/DatData/CorrectedTriggerFastMedian/B13289O14-DH2-Rec12' 'AC.xml']);
    end
%     if iD==7
%         Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec7' 'AC.xml']);
%     end
%     if iD==8
%         Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec8' 'AC.xml']);
%     end
%     if iD==9
%         Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec9' 'AC.xml']);
%     end
%     if iD==10
%         Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec10' 'AC.xml']);
%     end
%     if iD==11
%         Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec11' 'AC.xml']);
%     end
%     if iD==13
%         Par =  LoadXml(['../../../data/LargeScale/B13289O14-DH1-01463/DatData/CorrectedTriggerFastMedian/B13289O14-DH1-Rec13' 'AC.xml']);
%     end
%     Par =  LoadXml(['../../../data/ASIC1024/B14062W18-T1-rat01601/Rec_map_SIG_B14062W18-T1_54_S' '.xml']);

    d = dir(strcat(Directory, '*.dat'));%DC-LP30Hz-Notch50-100Hz.dat');

    LFPFs = 651.04166667;
    DownFact = 10;
    Tstab = 60;
    NCh = 256;
    
    %% get ephys data
    ACLfp = [];
    DCLfp = [];
    ClipLfp = [];

    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:NCh-1], Par.nChannels,1)';

        SplitName = split(d(fn).name,'-');
        
        if SplitName{3}(end-5:end-4)=='AC'
    %     myrLfp = resample(Lfp,1,10); % RESAMPLE TO 100 HZ !!!
            ACLfp = Lfp;
        elseif  SplitName{3}(end-5:end-4)=='DC' 
            DCLfp = Lfp;
        elseif SplitName{3}(end-7:end-4)=='Clip'
            ClipLfp = Lfp;
        end
    end

    %% rearrange and interpolate spatiall
    LfpGeom = zeros([size(ACLfp)]);
    LfpGeomDC = zeros([size(DCLfp)]);
    LfpGeomClip2 =  zeros([size(ClipLfp)]);

    for i = 1:length(Par.AnatGrps)
        for ii = 1:length(Par.AnatGrps(i).Channels)
            LfpGeom(:,ii+(i-1)*16) = ACLfp(:,Par.AnatGrps(i).Channels(ii)+1);
            LfpGeomDC(:,ii+(i-1)*16) = DCLfp(:,Par.AnatGrps(i).Channels(ii)+1);
            LfpGeomClip2(:,ii+(i-1)*16) = logical(ClipLfp(:,Par.AnatGrps(i).Channels(ii)+1)+Par.AnatGrps(i).Skip(ii));
        end
    end
    

    figure()

    LfpStruct = struct('DC',LfpGeomDC, 'AC', LfpGeom);
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
            if (mod(i,16) >= 2) 
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
            if (floor(i/16) >=1) & (floor(i/16) <= 14) & (mod(i,16)~=0) & (i<=16*7 | i>=16*9+1)
                mask2 = ((1-LfpGeomClip(:,i-16)).*(1-LfpGeomClip(:,i+16)));
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
                int2 = mask2.*(Lfp(:,i-16)+Lfp(:,i+16))/2;
            else
                display('in')
%                 display(min(mask2))
                mask2 = LfpGeomClip(:,19)*0;
                int2 = mask2;
            end
            if (floor(i/16) >=1) & (floor(i/16) <= 14) & (mod(i,16) >= 2) & (i<=16*7 | i>=16*9+1)
                mask3 = ((1-LfpGeomClip(:,i-16-1)).*(1-LfpGeomClip(:,i+16+1)));
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
                int3 = mask3.*(Lfp(:,i-16-1)+Lfp(:,i+16+1))/2;
            else
                display('in')
%                 display(min(mask3))
                mask3 = LfpGeomClip(:,19)*0;
                int3 = mask3;
            end
            if (floor(i/16) >=1) & (floor(i/16) <= 14) & (mod(i,16) >= 2) & (i<=16*7 | i>=16*9+1)
                mask4 = ((1-LfpGeomClip(:,i+16-1)).*(1-LfpGeomClip(:,i-16+1)));
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
                int4 = mask4.*(Lfp(:,i+16-1)+Lfp(:,i-16+1))/2;
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
    SaveBinary(strcat(Directory,'interDC.lfp'), LfpStruct.DC)
    SaveBinary(strcat(Directory,'interClip.lfp'), LfpClipInter)
end