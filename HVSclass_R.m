function [StatePeriods, StateMap, StateTitle, bands] = HVSclass_R(s, t, f, varargin)
%HVSclass_R is a function which determines HVS using LFP.
%This function is used by BrainStateClass_HVS.m which runs over a session data.
%
%USAGE: [HVSPeriods, HVSMap] = HVSclass_R(s, t, f, <FreqRangeDelta>, <FreqRangeSpindle>, <MinSpindleDeltaRatio>, <RatioSmoothWindowLength>)
%
%INPUT:
%  s              is a NxM matrix with a spectrogram
%  t              is a 1xN vector with time bins for the spectrogram (sec)
%  f              is a 1xM vector with frequency bins for the spectrogram (Hz)

% <FreqRangeDelta>     is a 1x2 vector with a delta frequency band (Hz). Default=[0 4] Hz.
% <FreqRangeSpindle>     is a 1x2 vector with a theta frequency band (Hz). Default=[6 30] Hz.
% <MinSpindleDeltaRatio> is the minimum spindle/delta ratio to be classified as HVS event. Default = 5.
% <RatioSmoothWindowLength> is a length (spec time bins) of a smoothing window for the spindle-delta ratio. Default = 2 bins.                           
% 
%                   
%OUTPUT: 
%  StatePeriods    is a 1x3 cell vector, where each cell contains a list of [start stop] pairs 
%                  for periods of THETA, SWS and REM brain states.
%  StateMap        is a 1xN vector with labeles of brain states at each time bin,
%                  (1-HVS 0-else).
%
%EXAMPLE:   [StatePeriods, StateMap] = HVSclass_R(spec.y, spec.t, spec.f, [t_spd MotorStateMap])
%
%DEPENDENCIES: 
% Labbox: DefaultArgs, contigous, MergeRangesShortGaps
% FMAToolbox: Interpolate, SpectrogramBands
% 
%
% The function was inspired by BrainStatesER from Evgeny Resnik
% Ramon Garcia 
% version 23.01.2022 
% 



if nargin<1
    error(['USAGE:  [StatePeriods, StateMap, StateTitle, bands] = BrainStatesER(s,t,f, <MotorStateMap>, <FreqRangeDelta>, <FreqRangeTheta>, <MinThetaDeltaRatio>, <RatioSmoothWindowLength>)'])
end


%Parse input parameters
[ FreqRangeDelta, FreqRangeSpindle, MinSpindleDeltaRatio, MinDurationHVS_sec, RatioSmoothWindowLength] = DefaultArgs(varargin,{ [0 4], [6 30], 1, 1, 2 });


% %--------------------------------------------------------------------------------------%
% %-----------------------  Constant parameters -----------------------------------------%
% % %Frequency bands for elta and theta rhythms
% % %Example of values used in FMAToolbox (0-4 Hz / 7-10 Hz) and CheckEegState (1-5 Hz / 5-12 Hz)
% % FreqRangeDelta = [1 4];
% % FreqRangeSpindle = [6 30];
% % %Minimum spindle/delta ratio to be classified as HVS event (otherwise none)
% % MinThetaDeltaRatio = 5;
% 
% 
% %Minimum duration of a HVS period 

% MinDurationHVS_sec = 1; 
% 
% %Maximum duration of a brief change in the brain state
% %(for now i decided to use a common theshold for all the brain states).
% MaxBriefChangeDuration_sec = 0; 
% 
% %Minimum peak value of theta/delta ratio during REM
% MinPeakRatioREM = 0;



%Compute instantenious power and power ratios in physiological bands

bands = SpectrogramBands(s,f, 'delta',FreqRangeDelta, 'spindles',FreqRangeSpindle, 'smooth', RatioSmoothWindowLength);
%Only theta/delta ratio kept in bands.ratio will be used.
ratio = bands.ratio;


%Initialize binary maps for brain states
HVS      = logical(zeros(size(t)));



%-----------------------------------------------------------------------------%
%       Brain state segmentation using thresholding of delta/theta ratio   
%-----------------------------------------------------------------------------%


%---------- Detect periods of suprathreshold theta/delta ratio ---------------%

%Create a binary map of the state with suprathreshold theta/delta ratio 

map = ratio > MinSpindleDeltaRatio;

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
if MinDurationHVS_sec >= t(2)-t(1)
    
    SampsMapFilter = floor(MinDurationHVS_sec/(t(2)-t(1)));
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




%Convert binary maps into periods (start/stop timestamps in sec)
if ~isempty(map)
    out = contiguous(map, 1);
    map_per = t(out{1,2});
    %convert a column to row when only one value is found
    if isequal(size(map_per), [2 1])
        map_per = map_per';
    end
else
    map_per = [];
end




%------------------------- Create output variables ----------------------------%
StateTitle = {'HVS'};
StatePeriods{1} = map_per;

%Create binary maps out of periods
StateMap = zeros(size(t,1),1);
for k=1:size(map_per,1)
    [~, ind] = ismember([map_per(k,1)  map_per(k,2)], t);
    StateMap( ind(1):ind(2) ) = 1;
end


clear ind
 


return



%-------------------------  Create output variables ----------------------------%
%Create output variables out of binary state maps
% StateTitle = {'THETA','SWS','REM','AwNonTheta'};
% 
% out = contiguous(theta,1);
% StatePeriods{1} = t(out{1,2});
% out = contiguous(sws,1);
% StatePeriods{2} = t(out{1,2});
% out = contiguous(rem,1);
% StatePeriods{3} = t(out{1,2});
% out = contiguous(AwNtheta,1);
% StatePeriods{4} = t(out{1,2});
% clear out 
% 
% StateMap = zeros(size(t,1),1);
% StateMap(theta) = 1;
% StateMap(sws)   = 2;
% StateMap(rem)   = 3;
% StateMap(AwNtheta) = 4;
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



