function [Data] = isi(Data,Settings)
% Code to fit displacement model to in-scan intervention data stored in the
% struct 'Data'.
% 
% Inputs: 
%   Data: struct holding PET scan data, with the fields: 
%       subID:      string with subject ID
%       t:          column array of times corresponding to blood data [min]
%       inFcn:      arterial input function, same size as 't' 
%       wb:         whole blood activity curve, same size as 't'
%       tPET:       kx1 array of PET frame mid times [min]
%       TACs:       kxR array of time-activity-curves for R regions
%       roiNames:   1xR cell array of region names corresponding to
%                   'TACs'
%       tb:         time of drug intervention [min]
%       scanDur:    time of the end of the PET scan [min]
%
%   Settings: struct holding options to pass to the solver, with the
%   fields: 
%       model:      string specifying the model, currently only '1tcm' is
%                   implemented
%       solver:     string specifying the solver, either 'numerical' or
%                   'singlestep'
%       doPlot:     set to 1 to plot fits, 0 otherwise 
%       weights:    string specifying frame weighting, either 'ones' or
%                   'framedur'. If not specified, the code sets the fit 
%                   parameters. 'ones' is the default. For duration-based 
%                   weights, the fields 'dur' (with frame durations 
%                   corresponding to tPET), or 'startEndTimes' should 
%                   preferrably be added to 'Data'. If not, the frame start 
%                   and end times will be estimated from 'Data.tPET'.  
%       fitParams:  Specifies options for fitting. Sub-fields: stepSize, 
%                   startParamsGlobal, upBoundBlobal, lowBoundGlobal, 
%                   startParamsROI, upBoundROI,lowBoundROI,options, 
%                   options.MaxFunEvals, options.MaxIter, options.TolFun, 
%                   options.TolX. 
%                   If not specified, the code sets the fit parameters.
%       
% Output: the code returns the struct 'Data', with additional field
% 'weights' and 'isi'.
% 'isi' has the sub-fields: 
%   model:          string specifying the model, either '1tcm' or '2tcm'
%   solver:         string specifying the solver, either 'numerical' or
%                   'singlestep'
% 	modelcurves:    kxR array of model fits to the data
%   VND:            estimated V_ND
%   occ:            estimated max occupancy reached during the scan
%   K1:             1xR arrray of estimated K_1s for each region 
%   VS:             1xR array of estimated V_Ss for each region 
%   VT:             1xR array of estimated V_Ts for each region 
%   vB:             1xR array of estimated fractional blood volume for 
%                   each region 
%   te:             estimated t_e (time when the occupancy occ is reached).
%                   Returned only if solver is 'numerical'.
%   ts:             estimated t_s (time of approximated step). Returned
%                   only if the solver is 'singlestep'.
%   exitFlag:       exit flag from the fit 
%
%__________________________________________________________________________
%                                   Gjertrud Louise Laurell & Martin Schain
%                                                   gjertrud.laurell@nru.dk
%                          Neurobiology Research Unit (NRU), Copenhagen, DK
%                                                                      2022

% SET FRAME WEIGHTS
% If nothing else is specified, all frames are equally weighted
if ~isfield(Settings,'weights') || isempty(Settings.weights)
    Settings.weights = 'ones'; 
end
switch Settings.weights
    case 'ones'
        Data.weights = ones(size(Data.tPET)); 
    case 'framedur'
        if isfield(Data,'dur')
            Data.weights = sqrt(Data.dur); 
        else
            if isfield(Data,'startEndTimes')
                startEndTimes = Data.startEndTimes; 
            else
                startEndTimes = getStartEndTimes(Data.tPET,Data.scanDur); 
            end
            frameDur = getFrameDurations(startEndTimes); 
            Data.weights = sqrt(frameDur); 
        end
end

% SET FIT PARAMETERS
if ~isfield(Settings,'fitParams')
    Settings.fitParams.options=optimoptions('lsqnonlin','display','off',...
        'DiffMinChange',1e-2);
    Settings.fitParams.options.MaxFunEvals = 10^4;
    Settings.fitParams.options.MaxIter = 10^4;
    Settings.fitParams.options.TolFun = 1e-6;
    Settings.fitParams.options.TolX = 1e-6;
    % Initiall guesses, upper and lower bounds for global parameters 
    % [occupancy, VND, t_e/t_s]
    Settings.fitParams.startParamsGlobal = [0.5 2 5];
    Settings.fitParams.upBoundGlobal = [1 inf inf];
    Settings.fitParams.lowBoundGlobal = [0 0 0];
    % Initiall guesses, upper and lower bounds for ROI-specific parameters 
    % [K1, VS, vB]
    Settings.fitParams.startParamsROI = [0.4 10 0.05];
    Settings.fitParams.upBoundROI = [inf inf 1];
    Settings.fitParams.lowBoundROI = [0 0 0];    
end

% FIT MODEL TO THE DATA
switch Settings.solver
    case 'numerical'
        Data = isiFit_numerical(Data,Settings); 
    case 'singlestep'
        Data = isiFit_singlestep(Data,Settings); 
end

% ILLUSTRATE RESULTS 
if Settings.doPlot
    R = length(Data.roiNames); 
    cols = lines(R);
    legText = {};
    figure('Position',[50 50 1100 800]), hold on 
    for roi = 1:R
        plot(Data.tPET,Data.TACs(:,roi),'ok',...
            'MarkerFaceColor',cols(roi,:),'MarkerSize',8)
        plot(Data.tPET,Data.isi.modelCurves(:,roi),'--',...
            'Color',cols(roi,:),'LineWidth',2)
        legText{(roi*2)-1} = Data.roiNames{roi}; 
        legText{roi*2} = [Data.roiNames{roi} ' fit']; 
    end
    title({[Settings.model ' displacement model (' Settings.solver ...
        ') fit to ' Data.subID],...
        ['Occupancy = ' num2str(Data.isi.occ*100,3) '%']})
    ylabel('Activity')
    xlabel('Time [min]')
    legend(legText)
    set(gca,'FontSize',18)
end

end

% SUPPORTING FUNCTIONS
function startEndTimes = getStartEndTimes(midTimes,scanDur)
% Approximate the start and end times of PET frames based on the frame
% mid-times. 

startEndTimes = zeros(length(midTimes)+1,1); 
for t = 1:length(midTimes)
    frDur = (midTimes(t) - startEndTimes(t))*2;
    startEndTimes(t+1) = startEndTimes(t) + frDur; 
    if t == length(midTimes)
        startEndTimes(t+1) = scanDur; 
    end
end

end

function frameDurs = getFrameDurations(startEndTimes)
% Calculates frame durations based on array of frame start and end times 
startTimes = startEndTimes(1:end-1);
endTimes = startEndTimes(2:end);
frameDurs = endTimes - startTimes; 
end
