close('all')


Directories = [
               "../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O24-DH1SL7-Rec3AC";
               "../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O24-DH1SL7-Rec4AC";
               "../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O24-DH1SL7-Rec5AC";
               "../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O24-DH1SL7-Rec1AC";
                ];
for iD = 1:length(Directories)
    Directory = Directories(iD);
    if iD==1
        Par =  LoadXml(['../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O14-DH1SL7-Rec3' 'AC.xml']);
    end
    if iD==2
        Par =  LoadXml(['../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O14-DH1SL7-Rec4' 'AC.xml']);
    end
    if iD==3
        Par =  LoadXml(['../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O14-DH1SL7-Rec5' 'AC.xml']);
    end
    if iD==4
        Par =  LoadXml(['../../../../data/LargeScale/B13289O24-DH1-01604/DatData/CorrectedTriggerFastMedian/B13289O14-DH1SL7-Rec1' 'AC.xml']);
    end
    d = dir(strcat(Directory, '*.dat'));%DC-LP30Hz-Notch50-100Hz.dat');

    LFPFs = 651.04166667;
    DownFact = 20;
    Tstab = 60;
    %% get ephys data
    ACLfp = [];
    DCLfp = [];
    ClipLfp = [];

    for fn = 1:length(d)
        FileName = [d(fn).folder ,'/',d(fn).name];
        Lfp = LoadBinaryDAT(FileName, [0:255], Par.nChannels,1)';

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

    %% rearrange and interpolate spatially
    LfpGeom = zeros([size(ACLfp)]);
    LfpGeomDC = zeros([size(DCLfp)]);
    LfpGeomClip2 =  zeros([size(ClipLfp)]);
    LfpGeomDepthAC = zeros([length(ACLfp) 32]);
    LfpGeomDepthDC = zeros([length(DCLfp) 32]);
    LfpGeomDepthClip = zeros([length(DCLfp) 32]);
    
    counter = 0;
    PositionDepth = [7 6 5 8 3 4 1 2];
    for i = 1:length(Par.AnatGrps)
        if length(nonzeros(i == [9 11 13 15]))>=1
            counter = counter+1;
        end
        for ii = 1:length(Par.AnatGrps(i).Channels)
            LfpGeom(:,ii+(i-1)*16) = ACLfp(:,Par.AnatGrps(i).Channels(ii)+1);
            LfpGeomDC(:,ii+(i-1)*16) = DCLfp(:,Par.AnatGrps(i).Channels(ii)+1);
            if length(nonzeros(i == [9 11 13 15]))>=1
                LfpGeomClip2(:,ii+(i-1)*16) = 1;
                LfpGeom(:,ii+(i-1)*16) = 0;
                LfpGeomDC(:,ii+(i-1)*16) = 0;
                if ii>=9
                    continue
                end
                LfpGeomDepthAC(:,counter+(ii-1)*4) = ACLfp(:,Par.AnatGrps(i).Channels(PositionDepth(ii))+1);
                LfpGeomDepthDC(:,counter+(ii-1)*4) = DCLfp(:,Par.AnatGrps(i).Channels(PositionDepth(ii))+1);
                LfpGeomDepthClip(:,counter+(ii-1)*4) = logical(ClipLfp(:,Par.AnatGrps(i).Channels(ii)+1)+Par.AnatGrps(i).Skip(ii));
            else         
                LfpGeomClip2(:,ii+(i-1)*16) = logical(ClipLfp(:,Par.AnatGrps(i).Channels(ii)+1)+Par.AnatGrps(i).Skip(ii));

            end
        end
    end
%     interpolate depth probe
    for i=1:32
        if (i>=1 & i<=31)
            if LfpGeomDepthClip(1,i)==1
                LfpGeomDepthAC(:,i) = (LfpGeomDepthAC(:,i-1) + LfpGeomDepthAC(:,i+1))/2;
                LfpGeomDepthDC(:,i) = (LfpGeomDepthDC(:,i-1) + LfpGeomDepthDC(:,i+1))/2;
            end
        end
    end
    
    figure()

    LfpStruct = struct('DC',LfpGeomDC, 'AC', LfpGeom);
    fn = fieldnames(LfpStruct);

    for k=1:numel(fn)
        if fn{k}=='AC'
            continue
        end
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

%         for i = 1:length(Lfp(1,:))
%             ind0 = find(LfpGeomClip(:,i)==1);
%             if length(ind0)>= 1
%                 ind1 = find(ind0 >= floor(60*LFPfs));
%                 if length(ind1)>= 1
%                     LfpGeomClip(ind0(ind1(1)):end,i) = 1;
%                 end
%             end
%         end
%         

        
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
%                 display(max(LfpGeomClip(:,i-1)))
%                 display(max(LfpGeomClip(:,i+1)))
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
%                 display('in')
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
%                 display('in')
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
%                 display('in')
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
%                 display('in')
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

%             display(min(Lfp(:,i)))
%             display(max(Lfp(:,i)))
%             display(max(abs(MaskSum)))

            Lfp(:,i) = (1-LfpGeomClip(:,i)).*Lfp(:,i) + Interp;
            LfpClipInter(:,i) = LfpGeomClip(:,i) & (1-mask1) & (1-mask2) & (1-mask3) & (1-mask4);
        end
        LfpClipInterSum = sum(LfpClipInter,2);

        plot(LfpClipInterSum)
        if fn{k} == 'DC'
            Lfp = [Lfp,LfpGeomDepthDC];
        else
            Lfp = [Lfp,LfpGeomDepthAC];
        end
        LfpStruct.(fn{k}) = Lfp;

    end
    %%
    clearvars -except LfpStruct Directory LfpClipInter
    
    SaveBinary(strcat(Directory,'interAC.lfp'), LfpStruct.AC)
    SaveBinary(strcat(Directory,'interDC.lfp'), LfpStruct.DC)
    SaveBinary(strcat(Directory,'interClip.lfp'), LfpClipInter)
end