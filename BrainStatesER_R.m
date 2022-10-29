function [StatePeriods, StateMap, StateTitle, bands] = BrainStatesER_R(s, t, f, varargin)
%BrainStatesER is a function which determines brain state (THETA, SWS. REM) using LFP and movement.
%This function is used by DetectBrainState.m which runs over a session data.
%
%USAGE: [StatePeriods, StateMap, StateTitle, bands] = BrainStatesER(s, t, f, <MotorStateMap>, <FreqRangeDelta>, <FreqRangeTheta>, <MinThetaDeltaRatio>, <RatioSmoothWindowLength>)
%
%INPUT:
%  s              is a NxM matrix with a spectrogram
%  t              is a 1xN vector with time bins for the spectrogram (sec)
%  f              is a 1xM vector with frequency bins for the spectrogram (Hz)
% <MotorStateMap> is a 2xN matrix, where the first colum contains time bins (can differ from the spectrogram time bins)
%                 and the second column contains labels of motor states at those time bin.
%                 The matrix must be created in advance by MotorStates.m (1-RUN, 2-QUIET, 3-SLEEP)
%                 Default = [].
% <FreqRangeDelta>     is a 1x2 vector with a delta frequency band (Hz). Default=[0 4] Hz.
% <FreqRangeTheta>     is a 1x2 vector with a theta frequency band (Hz). Default=[5 12] Hz.
% <MinThetaDeltaRatio> is the minimum theta/delta ratio to be classified as THETA state. Default = 5.
% <RatioSmoothWindowLength> is a length (spec time bins) of a smoothing window for the theta-delta ratio. Default = 2 bins.                           
% 
%                   
%OUTPUT: 
%  StatePeriods    is a 1x3 cell vector, where each cell contains a list of [start stop] pairs 
%                  for periods of THETA, SWS and REM brain states.
%  StateMap        is a 1xN vector with labeles of brain states at each time bin,
%                  (1-THETA, 2-SWS, 3-REM).
%
%EXAMPLE:   [StatePeriods, StateMap, StateTitle, bands] = BrainStatesER(spec.y, spec.t, spec.f, [t_spd MotorStateMap])
%
%DEPENDENCIES: 
% Labbox: DefaultArgs, contigous, MergeRangesShortGaps
% FMAToolbox: Interpolate, SpectrogramBands
% 
%
% The function was inspired by BrainStates from the FMAtoolbox (M. Zugaro)
% Evgeny Resnik
% version 02.12.2016 
% 
%
%NOTES:
%1) I decided to apply a fixed threshold for theta/delta ratio instead of clustering for the following reasons.
%   Clustering gives proper segregation between SWS and REM only if both states are present during SLEEP motor state. 
%   My tests have revealed, that when data have only SWS state, clustering may just further segregate SWS periods into substates.
%   Another issue, when clustering is done separately on SLEEP and RUN motor states, as expected, it gives very different threshold 
%   values for theta/delta ratio to distinguish SWS vs. REM and THETA vs. NONTHETA states.
%   Hence, i might come back to clustering as soon as we have EMG signal.
%2) Theta/theta_out ratio turned out to be less effective in state segregation as a conventional theta/delta ratio because
%   the band of 12-15 Hz is often contaminated by the first theta harmonic.
%3) In the case when animal motor data (speed or accelerometer) are not available, THETA and SWS states are detected purely based on 
%   theta/delta ratio without relation to actual animal motor activity. REM state can not be defined without motor data and will be empty.
%
%TO DISCUSS:
%1) Should we add additional validations of the detected REM periods. For example
%   - the maximum time interval between the preceeding SWS period and the REM period is <10sec; 
%   - there must be at least 120 sec of SLEEP motor state before the onset of the REM period;
%2) How should we name theta states: .THETA or .RUN? An option: RUN  is a simlink to THETA.
%   OR: motor states are saved into motor_run, motor_quiet, motor_sleep
%
%TO DO:
% - add HHM segmentation as in CheckEegStates




if nargin<1
    error(['USAGE:  [StatePeriods, StateMap, StateTitle, bands] = BrainStatesER(s,t,f, <MotorStateMap>, <FreqRangeDelta>, <FreqRangeTheta>, <MinThetaDeltaRatio>, <RatioSmoothWindowLength>)'])
end


%Parse input parameters
[MotorStateMap, FreqRangeDelta, FreqRangeTheta, MinThetaDeltaRatio, MinTimeRatio, RatioSmoothWindowLength] = DefaultArgs(varargin,{ [], [1 4], [5 12], 2, 20, 2 });


% %--------------------------------------------------------------------------------------%
% %-----------------------  Constant parameters -----------------------------------------%
% % %Frequency bands for elta and theta rhythms
% % %Example of values used in FMAToolbox (0-4 Hz / 7-10 Hz) and CheckEegState (1-5 Hz / 5-12 Hz)
% % FreqRangeDelta = [1 4];
% % FreqRangeTheta = [5 12];
% % %Minimum theta/delta ratio to be classified as THETA state (otherwise NONTHETA)
% % MinThetaDeltaRatio = 5;
% 
% 
% %Minimum duration of a single period of the particular brain state
% %(must be optimized so that even short but reasonable periods are kept).
% %(can differ for different brain states).
% MinDurationTHETA_sec = 0; 
% MinDurationREM_sec   = 0; 
% MinDurationSWS_sec   = 0;  
% 
% %Maximum duration of a brief change in the brain state
% %(for now i decided to use a common theshold for all the brain states).
% MaxBriefChangeDuration_sec = 0; 
% 
% %Minimum peak value of theta/delta ratio during REM
% MinPeakRatioREM = 0;


%--------------------------------------------------------------------------------------%
%Flag whether motor data available
if ~isempty(MotorStateMap)
    IfMotorDataAvailable = 1;
else
    IfMotorDataAvailable = 0;
end



if IfMotorDataAvailable==1
    
    %Interpolate motor state map at spectrogram timestamps
    
    MotorStateMap2 = Interpolate(MotorStateMap,t,'trim','off');
    MotorStateMap2(:,2) = round(MotorStateMap2(:,2));

 
    %Convert the common motor state map into individual state-specific binary maps
    run = zeros(size(t));
    run(MotorStateMap2(:,2)==1) = 1;
    run = logical(run);
    
    quiet = zeros(size(t));
    quiet(MotorStateMap2(:,2)==2) = 1;
    quiet = logical(quiet);
    
    sleep = zeros(size(t));
    sleep(MotorStateMap2(:,2)==3) = 1;
    sleep = logical(sleep);
    
end %if IfMotorDataAvailable==1


%Compute instantenious power and power ratios in physiological bands

bands = SpectrogramBands(s,f, 'delta',FreqRangeDelta, 'theta',FreqRangeTheta, 'smooth', RatioSmoothWindowLength);
%Only theta/delta ratio kept in bands.ratio will be used.
ratio = bands.ratio;


%Initialize binary maps for brain states
sws      = logical(zeros(size(t)));
rem      = logical(zeros(size(t)));
theta    = logical(zeros(size(t)));
nontheta = logical(zeros(size(t))); %will not be outputed for now
AwNTh = logical(zeros(size(t)));



%-----------------------------------------------------------------------------%
%       Brain state segmentation using thresholding of delta/theta ratio   
%-----------------------------------------------------------------------------%


%---------- Detect periods of suprathreshold theta/delta ratio ---------------%

%Create a binary map of the state with suprathreshold theta/delta ratio 

map = ratio > MinThetaDeltaRatio;

% merge closer than 20
% counter2 = 1;
% lenMap = length(map);
% for i= 1:lenMap-1
% 
%     if map(i+1,1)-map(i,2) <= 20
% 
%         map(i,2)=map(i+1,2);
%         map(i+1,:)=[];
%         counter2 = counter2+1;
% 
%         while map(i+1,1)-map(i,2) <= 20
%             counter2 = counter2+1; 
%         
%             map(i,2)=map(i+1,2);
%             map(i+1,:)=[];
%             
%         end
% 
%     end
% %     display(lenPerTHE-1-counter2)
%     if i >= lenMap-1-counter2
%         break
%     end
%     
% end

% % Merge map separated by less than MinTimeRatio and then eliminate the
% map shorter than MinTimeRatio
MinTimeRatio = MinTimeRatio;
if MinTimeRatio >= t(2)-t(1)
    
    SampsMapFilter = floor(MinTimeRatio/(t(2)-t(1)));
    map2(1:SampsMapFilter)= map(1:SampsMapFilter);


    for i = SampsMapFilter:length(map)-SampsMapFilter
        if (map(i-SampsMapFilter+1) == 1) && (map(i+SampsMapFilter)==1)
            map2(i-SampsMapFilter+1:i+SampsMapFilter) = 1;
        end
    end
    for i = SampsMapFilter:length(map)-SampsMapFilter
        if (map(i-SampsMapFilter+1) == 0) && (map(i+SampsMapFilter)==0)
            map2(i-SampsMapFilter+1:i+SampsMapFilter) = 0;
        end
    end      
    map = map2';
end

if length(run')>length(map)
    map= [map; zeros([length(run')-length(map), 1])];
else
    run= [run, zeros([length(map)-length(run'), 1])];
end
%Create binary maps of the brain states
if IfMotorDataAvailable==1
    display(length(run'))
    display(length(map))
    %THETA state periods: high theta/delta ratio periods during RUN
    
    theta = run' & map;

    
    %SWS state periods: low theta/delta ratio periods during SLEEP
    sws = sleep' & ~map;

    %REM state periods: high theta/delta ratio periods during SLEEP
    rem = sleep' & map;
    
    AwNTh = (~sleep' & ~map) | (quiet' & map);
else
    
    %if motor data not available, keep THETA and SWS detected pure based on theta/delta ratio
    %REM can not be defined without motor data.
    theta = map;
    sws   = ~map;
    rem   = [];
    AwNTh = [];
    
end %if IfMotorDataAvailable==1


%Convert binary maps into periods (start/stop timestamps in sec)
if ~isempty(theta)
    out = contiguous(theta, 1);
    theta_per = t(out{1,2});
    %convert a column to row when only one value is found
    if isequal(size(theta_per), [2 1])
        theta_per = theta_per';
    end
else
    theta_per = [];
end

if ~isempty(sws)
    out = contiguous(sws, 1);
    sws_per = t(out{1,2});
    %convert a column to row when only one value is found
    if isequal(size(sws_per), [2 1])
        sws_per = sws_per';
    end
else
    sws_per = [];
end

if ~isempty(rem)
    out = contiguous(rem, 1);
    rem_per = t(out{1,2});
    %convert a column to row when only one value is found
    if isequal(size(rem_per), [2 1])
        rem_per = rem_per';
    end
    
else
    rem_per = [];
end

if ~isempty(AwNTh)
    out = contiguous(AwNTh, 1);
    aw_Ntheta = t(out{1,2});
    %convert a column to row when only one value is found
    if isequal(size(aw_Ntheta), [2 1])
        aw_Ntheta = aw_Ntheta';
    end
    
else
    aw_Ntheta = [];
end
clear map out theta sws rem




% %------------------- Discard VERY short (<0.5s) state periods -------------------------%
% %Just to avoid errors later due to periods of 1 time bin duration.
% %Main sorting out short periods will be done below.
% MinDurationTOTAL_sec = 0;
% BadPer = diff(theta_per,1,2) < MinDurationTOTAL_sec;
% theta_per(BadPer,:) = [];
% BadPer = diff(sws_per,1,2) < MinDurationTOTAL_sec; 
% sws_per(BadPer,:) = [];
% BadPer = diff(rem_per,1,2) < MinDurationTOTAL_sec; 
% rem_per(BadPer,:) = [];
% clear BadPer MinDurationTOTAL_sec



% %------------------ Drop out brief changes in brain states -------------------%
% %NOTE: 
% % Droping out means that two state periods separated by the brief gap will be merged.
% %NOTE: 
% % brief changes in the brain state are droped out only if the motor state did not change during that change!
% %IMPORTANT: 
% % Droping out brief gaps between long periods of a particular state contaminates initially detected states.
% % Whether to do it or not depends on what is more critical: to get periods of THETA state or periods of 
% % high-power theta oscillations.
% 
% %Loop across states: apply the same period selection criteria to all states
% for st=1:3
%     
%     %pick the periods from a given state
%     switch st
%         case 1
%             state_per = theta_per;
%         case 2
%             state_per = sws_per;
%         case 3
%             state_per = rem_per;
%     end
%     
%     if isempty(state_per)
%         continue
%     end
%                
%     
%     if IfMotorDataAvailable==1
%         %Drop out brief changes in the brain state only if the motor state did not change
% 
%         %Create a matrix with indices of period pairs which should NOT be merged
%         %despite the short gap between them, because the motor state changed
%         Periods2Keep = [];
%         for k=1:size(state_per,1)-1
%             %first period
%             ind = find(  MotorStateMap2(:,1) > state_per(k,1)  &  MotorStateMap2(:,1) < state_per(k,2) );            
%             PerMotorState = unique(MotorStateMap2(ind,2));            
%             
%             %Deal with the case when the given state period state_per(k,:) is so short, 
%             %that not a single motor state samples lies within this period.
%             %Do not keep such short brain state period.
%             if isempty(ind)                
%                 continue
%             end
%             
%             %gap between first and second periods
%             ind = find(  MotorStateMap2(:,1) > state_per(k,2)  &  MotorStateMap2(:,1) < state_per(k+1,1) );
%             GapMotorState = unique(MotorStateMap2(ind,2));
%             %State periods by definition must have a unique motor state
%             if length(PerMotorState)>1
%                 error('Brain state period with multiple motor states found! It should not happen! Check for bugs! ')
%             end
%             %Do not merge two periods of the given state if the motor state changed in the gap between them
%             if length(GapMotorState)>1 || GapMotorState~=PerMotorState
%                 Periods2Keep = [Periods2Keep; k k+1];
%             end
%             clear ind PerMotorState GapMotorState
%         end %loop across state periods
%         
%         
%     elseif IfMotorDataAvailable==0 
%         %Drop out brief changes WITHOUT taking into account the motor state
%         Periods2Keep = [];
%         
%     end %if IfMotorDataAvailable==1
%     
%     
%     %Drop out brief changes in the brain state
%     state_per = MergeRangesShortGaps(state_per, MaxBriefChangeDuration_sec, Periods2Keep);
%     
%     
%     %Keep the outcome in the state-corresponding variable
%     switch st
%         case 1
%             theta_per = state_per;
%         case 2
%             sws_per = state_per;
%         case 3
%             rem_per = state_per;
%     end    
%         
%     clear state_per Periods2Keep   
% end %loop across states
% 
% 
% 
% %------------------- Discard too short state periods -------------------------%
% BadPer = diff(theta_per,1,2) < MinDurationTHETA_sec;
% theta_per(BadPer,:) = [];
% BadPer = diff(sws_per,1,2) < MinDurationSWS_sec; 
% sws_per(BadPer,:) = [];
% BadPer = diff(rem_per,1,2) < MinDurationREM_sec; 
% rem_per(BadPer,:) = [];
% clear BadPer
% 
% 
% 
% %----- Discard REM periods with too low PEAK value of theta/delta ratio ------%
% if ~isempty(rem_per)
%     
%     %Compute the peak theta/delta ratio for the detected REM periods
%     for k=1:size(rem_per,1)
%         [~, ind] = ismember([rem_per(k,1)  rem_per(k,2)], t);
%         PeakRatioREM(k) = max( ratio(ind(1):ind(2)) );
%     end
%     
%     %discard 'bad' REM periods
%     BadPer = PeakRatioREM < MinPeakRatioREM;
%     rem_per(BadPer,:) = [];
%     clear BadPer k ind
%     
% end %if ~isempty(rem_per)
% 
% 
% %------------ Discard REM periods surrounded by THETA periods ----------------%
% %ASSUMPTION: REM periods can occur only in between epochs of SWS with short time gaps between them.
% %Hence, 
% % 1) the period right before ***and the one right after(not sure we shoul keep it???)*** a given REM period must be either REM or SWS.
% % 2) Time interval between these surronding SWS/REM and the given REM period must be short enough (say, <10 sec).
% %    This second claim is questinable, because it may happen, that a continuous SWS epoch will be interrupted 
% %    by a brief state change, so that the second piece of the SWS epoch turn out to be too short to be pass 
% %    the selection. As a result there will be a long time gap between the SWS epoch and the given REM period.
% 
% if ~isempty(rem_per)
%     
%     %Combine all periods keeping their identity (just for convenience)
%     allper = [theta_per; sws_per; rem_per];
%     allper_state = [ 1*ones(size(theta_per,1),1);  2*ones(size(sws_per,1),1);  3*ones(size(rem_per,1),1) ];
%     
%     %Sort combined periods by the period start in an asceding order
%     [~, ind] = sort(allper(:,1));
%     allper = allper(ind,:);
%     allper_state = allper_state(ind);
%     clear ind
%     
%     %Loop across REM periods
%     BadPer = [];
%     clear TimeGapToREM
%     for k=1:size(rem_per,1)
%         %find the period closest to the given REM period on the left side
%         indL = find( allper(:,2) <= rem_per(k,1)  );
%         indL = max(indL);
%         %find the period closest to the given REM period on the right side
%         indR = find( allper(:,1) >= rem_per(k,2)  );
%         indR = min(indR);
%         
%         %NEW: Time interval between the given REM period and preeceding SWS/REM period
%         %TimeGapToREM(k) = rem_per(k,1) - allper(indL,2);
%         
%         %accumulate a list of REM periods to be discarded
%         if  allper_state(indL)==1 |  allper_state(indR)==1
%             %fprintf('DEBUGGING: bad REM period-%d: %d-REM-%d, (t=%1.0f)\n', k,  allper_state(indL),  allper_state(indR), rem_per(k,1)  )
%             BadPer = [BadPer;  k];
%             display(k)
%         end
%         clear indL indR
%     end %loop across REM periods
%     
%     
%     % %NEW: Maximum time interval between the onset/offset of the given REM period and the preceeding SWS/REM period
%     % MaxTimeGapToREM_sec = 30;
%     % %NEW: find REM periods which are too apart from the preceeding SWS/REM period
%     % BadPer2 = find( TimeGapToREM > MaxTimeGapToREM_sec );
%     
%     %discard identified 'bad' REM periods
%     rem_per(BadPer,:) = [];
%     
%     clear allper allper_state
%     
% end %if ~isempty(rem_per)


%------------------------- Create output variables ----------------------------%
StateTitle = {'THETA','SWS','REM'};
StatePeriods{1} = theta_per;
StatePeriods{2} = sws_per;
StatePeriods{3} = rem_per;
StatePeriods{4} = aw_Ntheta;
%Create binary maps out of periods
StateMap = zeros(size(t,1),1);
for k=1:size(theta_per,1)
    [~, ind] = ismember([theta_per(k,1)  theta_per(k,2)], t);
    StateMap( ind(1):ind(2) ) = 1;
end
for k=1:size(sws_per,1)
    [~, ind] = ismember([sws_per(k,1)  sws_per(k,2)], t);
    StateMap( ind(1):ind(2) ) = 2;
end
for k=1:size(rem_per,1)
    [~, ind] = ismember([rem_per(k,1)  rem_per(k,2)], t);
    StateMap( ind(1):ind(2) ) = 3;
end

for k=1:size(aw_Ntheta,1)
    [~, ind] = ismember([aw_Ntheta(k,1)  aw_Ntheta(k,2)], t);
    StateMap( ind(1):ind(2) ) = 4;
end

clear ind
 


return






% %-----------------------------------------------------------------------------%
% %          OPTIONAL Figure:   Brain state segmentation FIXED THRESHOLD
% %-----------------------------------------------------------------------------%
% figure
% 
% % delete(h)
% h = tight_subplot(4,1, [.04 .01], [0.1 0.1], [0.04 0.02]);
% 
% %Plot spectrogram with detected states
% for st=1:length(StateTitle)
%     
%     axes(h(st)); cla; hold on
%     imagesc(t, f, log(s));
%     axis tight; axis xy;
%     set(gca, 'ylim', [1 20], 'ytick', [5:5:20], 'yticklabel', [5:5:20])
%     title(sprintf('%s periods', StateTitle{st}))
%     
%     %for color scaling use median+/-std across all freq. bins
%     ry = reshape(log(s),[],1);
%     ry = ry(~isinf(ry) & ~isnan(ry));
%     med = median(ry);
%     dev = std(ry);
%     if isnan(dev) dev=med/3; end
%     caxis(med+[-3 3]*dev);
%     clear ry med dev
%     
%     %label delta and theta frequency band
%     if ~isempty(FreqRangeDelta)
%         plot(xlim, FreqRangeDelta(2)*[1 1], 'k:')
%     end
%     if ~isempty(FreqRangeTheta)
%         plot(xlim, FreqRangeTheta(1)*[1 1], 'k:')
%         plot(xlim, FreqRangeTheta(2)*[1 1], 'k:')
%     end    
%     
%     %label motor states with colors: 1-RUN, 2-QUIET, 3-SLEEP
%     if exist('MotorStateMap0', 'var')
%         out = contiguous(MotorStateMap2(:,2), [1 2 3]);
%         for k=1:3
%             MotorPeriods{k} = t(out{k,2});
%             if isequal(size(MotorPeriods{k}), [2 1])
%                 MotorPeriods{k} = MotorPeriods{k}';
%             end
%         end
%         yLim = ylim;
%         MotorStateColor = {'r','g','b'};
%         for ms=1:length(MotorPeriods)
%             for n=1:size(MotorPeriods{ms},1)
%                
%                 if ~isempty(MotorPeriods{ms})
%                     p(ms) = plot(MotorPeriods{ms}(n,:), 0.8*yLim(2)*[1 1], 'Color',MotorStateColor{ms}, 'LineWidth', 2);                    
%                     PresentState(ms) = logical(1);  
%                 else
%                     PresentState(ms) = logical(0);
%                 end                           
%                 
%             end %loop n
%         end  %loop ms 
%         
%         %Create legend labels in accordance with the present states
%         LegendStr = {'RUN','QUIET','SLEEP'};
%         LegendStr(~PresentState) = [];
%         p(~PresentState) = [];
%         clear PresentState
%         
%         if st==1
%             legend(p, LegendStr)
%         end
%         clear out k yLim ms
%     end %if exist('MotorStateMap0', 'var')
%     
%     
%     %Label brain states with period borders
%     if ~isempty(StatePeriods{st})
%         PlotIntervals(StatePeriods{st}, 'bars');
%     end
%     
%     %label only on the upper subplot
%     if st==1
%         ylabel('Frequency, (Hz)');
%     end    
%     
%     linkx    
% end %loop across brain states
% 
% 
% 
% %plot theta/delta ratio
% axes(h(length(StatePeriods)+1)); cla; hold on
% plot(t, bands.ratio, 'k')
% axis tight; 
% yLim = [0 round(max(bands.ratio))];
% ylbl = 0:10:yLim(2);
% set(gca, 'ylim', yLim, 'ytick', ylbl, 'yticklabel', ylbl)
% xlbl = round(linspace(0, t(end), 10));
% set(gca,'xlim',t([1 end]), 'xtick',xlbl, 'xticklabel',xlbl)
% ylabel('Theta/delta ratio');
% xlabel('Time, (s)')
% linkx
% %label thtreshold
% plot(xlim, MinThetaDeltaRatio*[1 1], 'k--')
% %Label brain states with lines
% % PlotIntervals(theta_per, 'bars'); 
% % PlotIntervals(sws_per, 'bars'); 
% % PlotIntervals(rem_per, 'bars'); 
% clear xlbl yLim ylbl
% 
% 
% 
% % %---------------------------- Suptitle ---------------------------------------------%
% % titlestr = sprintf('Identification of anatomical layers: %s, %s, %s, %s, RipChan=%d, ThChan=%d, Ch=%d-%d', ...
% %     FileBase, SignalType, PeriodTitle, SleepPeriodTitle, RipChan, ThChan, Channels([1 end]) );
% % suptitle2(titlestr, .98, 1)
% % clear titlestr
% % 
% % 
% % %Save figure into a file
% % FileOut = sprintf(['%s.%s.%s.%s.%s.%d-%d.%d'], FileBase, mfilename, SignalType, PeriodTitle, SleepPeriodTitle, Channels([1 end]), FigIndex);
% % fprintf(['Saving figure  into a file %s ...'], [FileOut '.jpg'])
% % AxisPropForIllustrator
% % FigPathCommon = [pwd '/figures/'];
% % mkdir2([ FigPathCommon mfilename ])
% % SaveForIllustrator( [FigPathCommon mfilename '/' FileOut],'jpg')
% % % SaveForIllustrator( [FigPathCommon mfilename '/' FileOut],'eps')
% % fprintf('DONE\n')
% 











% %-----------------------------------------------------------------------------%
% %                             NOT IN USE
% %-----------------------------------------------------------------------------%
% %For now i decided no to use clustering, because it gives less predictable outcome.
% 
% %-----------------------------------------------------------------------------%
% %          WAY-1:   Brain state segmentation by clustering 
% %-----------------------------------------------------------------------------%
% %Determine features for automatic data clustering (from the whole session) 
% switch(method)
%     case 'pca',
%         % Compute PCA on spectrogram and reduce dimensionality
%         S = s';
%         S(:,f>30) = 0;
%         [eigenvectors,projected,lambda] = princomp(S,'econ');
%         if nComponents == 0,
%             nComponents = find(cumsum(lambda)/sum(lambda)>0.85);
%             nComponents = nComponents(1);
%         end
%         eigenvectors = eigenvectors(:,1:nComponents);
%         lambda = lambda(1:nComponents);
%         projected = projected(:,1:nComponents);
%         features = projected;
%     case 'ratios',
%         % Use heuristic ratios
%         features = [bands.ratio1  bands.ratio2];        
%     case 'direct',
%         % Use theta/delta ratio
%         features = bands.ratio;        
%     case 'direct_extended',
%         % Use theta/delta_out ratio (like in CheckEegState)
%         %Turned out to be less efficient because of theta harmonic contamination.
%         features = bands.ratio3;
% end %switch(method)
% 
% 
% %Add EMG to feature list
% if ~isempty(emg0),
% 	features = [features emg0];
% end
% 
% 
% %Normalize columns (we need this because K-means clustering uses Euclidean distances)
% features = zscore(features);
% 
% 
% 
% %--------------- 1. Cluster SWS/REM states during SLEEP  -----------------%
% %Keep features only from SLEEP 
% features0 = features(sleep,:);
% 
% %Cluster (K-means)
% cluster = kmeans(features0, nClusters);
% 
% %Use either theta/delta or theta/theta_out ratio to distinguish SWS/REM
% switch(method)
%     case 'direct',
%         %theta/delta ratio
%         r = bands.ratio(sleep);
%     case 'direct_extended',
%         %theta/delta_out ratio
%         r = bands.ratio3(sleep);
% end
% 
% %Compute the mean theta/delta ratio for both detected states (clusters)
% for i = 1:nClusters,
% 	ratio(i) = mean(r(cluster==i));
% end
% 
% %Apply some heuristic analysis
% %1) REM: highest theta/delta ratio during SLEEP
% [~,i] = sort(ratio(:),1,'descend');
% h = i(1);
% rem(sleep) = cluster==h;
% highest = ratio(h);
% 
% %2) SWS: lowest theta/delta ratio during SLEEP
% [~,i] = sort(ratio(:),1,'ascend');
% l = i(1);
% sws(sleep) = cluster==l;
% lowest = ratio(l);
% 
% %3) REM must have a theta/delta ratio at least twice as high as SWS
% %It deals with cases when only one state (SWS or REM) is present and discard, for example, 
% %detected REM periods if their mean theta/delta ratio is not distinct enough relative to the one in SWS.
% %The threshold value here is taken as in the FMAToolbox.
% MinStateRatioDifference = 1; % this was 2 by default
% if highest < MinStateRatioDifference * lowest,
% 	sws(sleep) = 1;
% 	rem(sleep) = 0;
% end
% 
% 
% % %OPTIONAL VARIABLE: just for debugging and tunning purposes
% % %Transition value of theta/delta ratio produced by clustering.
% % %Periods with theta/delta ratio below this value have been assigned to SWS state.
% % %Periods with theta/delta ratio above this value have been assigned to REM state.
% % TransitionValue_SWSREM = mean([ min(r(cluster==h)) max(r(cluster==l)) ]);
% 
% 
% clear features0 clusters lowest highest r h l 
% 
% 
% 
% 
% %--------------- 2. Cluster THETA/NONTHETA states during RUN  -----------------%
% %Keep features only from SLEEP 
% features0 = features(run,:);
% 
% %Cluster (K-means)
% cluster = kmeans(features0, nClusters);
% 
% %Use either theta/delta or theta/theta_out ratio to distinguish THETA/NONTHETA
% switch(method)
%     case 'direct',
%         %theta/delta ratio
%         r = bands.ratio(run);
%     case 'direct_extended',
%         %theta/delta_out ratio
%         r = bands.ratio3(run);
% end
% 
% 
% %Compute the mean theta/delta ratio for both detected states (clusters)
% for i = 1:nClusters,
% 	ratio(i) = mean(r(cluster==i));
% end
% 
% %Apply some heuristic analysis
% %1) THETA: highest theta/delta ratio during RUN
% [~,i] = sort(ratio(:),1,'descend');
% h = i(1);
% theta(run) = cluster==h;
% 
% %3) NONTHETA: lowest theta/delta ratio during RUN
% [~,i] = sort(ratio(:),1,'ascend');
% l = i(1);
% nontheta(run) = cluster==l;
% % lowest = ratio(l);
% 
% 
% % %OPTIONAL VARIABLE: just for debugging and tunning purposes
% % %Transition value of theta/delta ratio produced by clustering.
% % %Periods with theta/delta ratio below this value have been assigned to SWS state.
% % %Periods with theta/delta ratio above this value have been assigned to REM state.
% % TransitionValue_THETA = mean([ min(r(cluster==h)) max(r(cluster==l)) ]);
% 
% 
% clear clusters features0 r h l lowest highest



%-------------------------  Create output variables ----------------------------%
%Create output variables out of binary state maps
StateTitle = {'THETA','SWS','REM','AwNonTheta'};

out = contiguous(theta,1);
StatePeriods{1} = t(out{1,2});
out = contiguous(sws,1);
StatePeriods{2} = t(out{1,2});
out = contiguous(rem,1);
StatePeriods{3} = t(out{1,2});
out = contiguous(AwNtheta,1);
StatePeriods{4} = t(out{1,2});
clear out 

StateMap = zeros(size(t,1),1);
StateMap(theta) = 1;
StateMap(sws)   = 2;
StateMap(rem)   = 3;
StateMap(AwNtheta) = 4;
%-----------------  END OF WAY-1: clustering ---------------------------------------%
%-----------------------------------------------------------------------------------%










%-----------------------------------------------------------------------------%
%           Supplementary functions to stay independent from the labbox
%-----------------------------------------------------------------------------%

% function runs = contiguous(A,varargin)
% %   RUNS = CONTIGUOUS(A,NUM) returns the start and stop indices for contiguous 
% %   runs of the elements NUM within vector A.  A and NUM can be vectors of 
% %   integers or characters.  Output RUNS is a 2-column cell array where the ith 
% %   row of the first column contains the ith value from vector NUM and the ith 
% %   row of the second column contains a matrix of start and stop indices for runs 
% %   of the ith value from vector NUM.    These matrices have the following form:
% %  
% %   [startRun1  stopRun1]
% %   [startRun2  stopRun2]
% %   [   ...        ...  ]
% %   [startRunN  stopRunN]
% %
% %   Example:  Find the runs of '0' and '2' in vector A, where
% %             A = [0 0 0 1 1 2 2 2 0 2 2 1 0 0];  
% %    
% %   runs = contiguous(A,[0 2])
% %   runs = 
% %           [0]    [3x2 double]
% %           [2]    [2x2 double]
% %
% %   The start/stop indices for the runs of '0' are given by runs{1,2}:
% %
% %           1     3
% %           9     9
% %          13    14
% %
% %   RUNS = CONTIGUOUS(A) with only one input returns the start and stop
% %   indices for runs of all unique elements contained in A.
% %
% %   CONTIGUOUS is intended for use with vectors of integers or characters, and 
% %   is probably not appropriate for floating point values.  You decide.  
% %
% 
% if prod(size(A)) ~= length(A),
%     error('A must be a vector.')
% end
% 
% if isempty(varargin),
%     num = unique(A);
% else
%     num = varargin{1};
%     if prod(size(num)) ~= length(num),
%         error('NUM must be a scalar or vector.')
%     end
% end
% 
% for numCount = 1:length(num),
%     
%     indexVect = find(A(:) == num(numCount));
%     shiftVect = [indexVect(2:end);indexVect(end)];
%     diffVect = shiftVect - indexVect;
%     
%     % The location of a non-one is the last element of the run:
%     transitions = (find(diffVect ~= 1));
%     
%     runEnd = indexVect(transitions);
%     runStart = [indexVect(1);indexVect(transitions(1:end-1)+1)];
%     
%     runs{numCount,1} = num(numCount);
%     runs{numCount,2} = [runStart runEnd];
%     
% end





%     %-------- If motor data are available --------%
%     %Compute inter-period time intervals between successive periods of the given state
%     IPI = state_per(2:end,1) - state_per(1:end-1,2);
%     
%     %Find 'brief' inter-period time intervals (=gaps)
%     BadIPI = find( IPI < MaxBriefChangeDuration_sec );
%     
%     %Loop across the detected gaps
%     for k=1:length(BadIPI)        
%         
%         %Motor states the two brain state periods and the gap between them occured at
%         %NOTE: Both brain state periods must have same motor state!
%         %NOTE: GapMotorState may have a more than one motor state!
%         [~, ind] = ismember([state_per(BadIPI(k),1)  state_per(BadIPI(k),2)], t);
%         PerMotorState = unique( MotorStateMap0( ind(1):ind(2) ,2) );
%         [~, ind] = ismember([state_per(BadIPI(k),2)  state_per(BadIPI(k)+1,1)], t);
%         GapMotorState = unique( MotorStateMap0( ind(1):ind(2) ,2) );
%         
%         if length(PerMotorState)>1
%             error('Brain state period with multiple motor states found! It should not happen! Check for bugs! ')
%         end
%         
%         %Do not merge two periods of the given state if the motor state changed in the gap between them
%         if length(GapMotorState)>1 | GapMotorState~=PerMotorState
%             %fprintf('DEBUGGING: skipped short gap-%d (t=%1.0f)\n', BadIPI(k), state_per(BadIPI(k),2)  )            
%             continue
%         end
%         
%         %create a new merged period by dropping out the gap between them
%         newper = [state_per(BadIPI(k),1)  state_per(BadIPI(k)+1,2)];
%         
%         %replace both original periods with the same new one
%         %(it is necessary for proper processing of gaps, repeated periods will be removed later)
%         state_per(BadIPI(k),:)   = newper;
%         state_per(BadIPI(k)+1,:) = newper;
%         
%         clear ind PerMotorState GapMotorState newper
%         
%     end %loop across 'brief' inter-period time intervals
%     
%     
%     %Remove repeated state periods, which appeared after merging
%     state_per = unique(state_per, 'rows');
%     
%        
%     %As a result of merging state periods with the same start but different stop can appear.
%     %If it happens, keep only the longest one.
%     %compute duration of the state periods
%     dt = diff(state_per,1,2);
%     %list of unique period start without repetitions
%     ulist = unique(state_per(:,1));
%     %loop across unique period starts
%     BadPer=[];
%     for k=1:length(ulist)
%         %indicies of state periods with the same start
%         ind = find( state_per(:,1) == ulist(k) );
%         %if more than one state periods with the common start found, keep the longest one (must always be the last one)
%         if length(ind)>1
%             [~, LongestPerInd] = max( dt(ind) );
%             OtherPerInd = setdiff( 1:length(ind) , LongestPerInd );
%             %accumulate a list of repeated periods to be removed
%             BadPer = [BadPer;  ind(OtherPerInd)];
%         end
%         clear ind LongestPerInd OtherPerInd
%     end %loop across unique period starts
%     
%     %Remove repeated periods with the same starts
%     state_per(BadPer,:) = [];
% 
% %Just to doube-check
% if size(state_per,1) ~= length( unique(state_per(:,2)) )
%     error('Brain state periods with same stop found! It should not happen! Check for bugs! ')
% end



