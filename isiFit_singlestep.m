function [Data] = isiFit_singlestep(Data,Settings)
% Fit displacement model to in-scan intervention data analytically, using
% the single-step approximation 
%
%__________________________________________________________________________
%                             Gjertrud Louise Laurell & Martin Schain, 2022

% SET STEP SIZE FOR SINGLE STEP
h = 1/60; 

% UPDATE BLOOD DATA ARRAYS TO SINGLE STEP STEP SIZE (h)
if Data.t(1) > 0   % Ensure that blood data arrays starts at time = 0 min  
    Data.t = [0; Data.t];	
    Data.inFcn = [0; Data.inFcn]; 	
    Data.wb = [0; Data.wb]; 
end
oldTime = Data.t;   oldAIF = Data.inFcn;    oldWB = Data.wb; 
Data.t = (0:h:Data.scanDur)'; 
Data.inFcn = ...
    interp1(oldTime,oldAIF,Data.t,'pchip',mean(oldAIF(end-59:end)));
Data.wb = interp1(oldTime,oldWB,Data.t,'pchip',mean(oldWB(end-59:end)));

switch Settings.model
    case '1tcm'
        [models, K1, VS, vB, occ, VND, ts, exitFlag] = ...
            runSingleStepISI_1TCM(Data,Settings); 
end

% RESTORE BLOOD DATA ARRAYS
Data.t = oldTime;   Data.inFcn = oldAIF;    Data.wb = oldWB; 

% SAVE OUTPUTS IN FIELD 'isi' IN 'Data'
isi.model = Settings.model; 
isi.solver = Settings.solver; 
isi.weights = Settings.weights; 
isi.modelCurves = models; 
isi.occ = occ;          isi.VND = VND;          isi.K1 = K1; 
isi.VS = VS;            isi.VT = VND + VS;      isi.vB = vB; 
isi.ts = ts;            isi.exitFlag = exitFlag; 

Data.isi = isi;
end

% SUPPORTING FUNCTIONS - 1TCM CASE
function [models, K1, VS, vB, occ, VND,ts,exitFlag] = ...
    runSingleStepISI_1TCM(Data,Settings)
% Fit data to 1TC displacement model

% Fit outer layer parameters (occupancy, VND & ts)
[pOccVndTs, ~, ~, exitFlag] = ...
    lsqnonlin(@fitSS_1TC_outer,Settings.fitParams.startParamsGlobal,...
    Settings.fitParams.lowBoundGlobal,Settings.fitParams.upBoundGlobal,...
    Settings.fitParams.options,Data,Settings); 

% Pre-allocate model curves and inner parameter array
models = zeros(size(Data.TACs)); 
innerPars = zeros(size(Data.TACs,2),3); 

% Fit inner layer parameters with outer layer fixed
for roi = 1:size(Data.TACs,2)
    [pK1VsvB, ~, ~, ~] = ...
        lsqnonlin(@fitSS_1TC_inner,Settings.fitParams.startParamsROI,...
        Settings.fitParams.lowBoundROI,Settings.fitParams.upBoundROI,...
        Settings.fitParams.options,Data,roi,pOccVndTs); 
    pTot = [pK1VsvB pOccVndTs]; 
    modelCurve = SS_1TCM_createCurves(pTot,Data); 
    models(:,roi) = interp1(Data.t,modelCurve(:),Data.tPET); 
    innerPars(roi,:) = pK1VsvB; 
end

% Assign out-parameters
K1 = innerPars(:,1);    VS = innerPars(:,2);    vB = innerPars(:,3); 
occ = pOccVndTs(1);     VND = pOccVndTs(2);     ts = pOccVndTs(3); 
end

function [err] = fitSS_1TC_outer(currP_global,Data,Settings)
% Outer layer for 1TCM - fits global parameters  (occupancy, VND & ts)
err = zeros(size(Data.TACs)); 
for roi = 1:size(Data.TACs,2)
    [pK1VsvB, ~, ~, ~] = ...
        lsqnonlin(@fitSS_1TC_inner,Settings.fitParams.startParamsROI,...
        Settings.fitParams.lowBoundROI,Settings.fitParams.upBoundROI,...
        Settings.fitParams.options,Data,roi,currP_global); 
    P = [pK1VsvB currP_global]; 
    modelCurve = SS_1TCM_createCurves(P,Data); 
    modelTAC = interp1(Data.t,modelCurve(:),Data.tPET); 
    err(:,roi) = (modelTAC - Data.TACs(:,roi)).*Data.weights(:);
end
end

function [err] = fitSS_1TC_inner(p,Data,currROI,p_global)
% Inner layer for 1TCM - fits ROI-specific parameters (K1, VS and vB)
P = [p p_global]; 
modelCurve = SS_1TCM_createCurves(P,Data); 
modelTAC = interp1(Data.t,modelCurve(:),Data.tPET); 
err = (modelTAC - Data.TACs(:,currROI)).*Data.weights(:); 
end

function [modelCurve] = SS_1TCM_createCurves(P,Data)
% Draw the model curves for the current choice of parameter values using
% Single step

% Set step size 
h = 1/60; 

% Time data 
t = Data.t;         tb = Data.tb; 

% Current parameter values
k1 = P(1);          vs = P(2);          vb = P(3);
o = P(4);           vnd = P(5);         ts = tb + P(6); 
k2 = k1/vnd;        bp = vs/vnd; 

% Blood data 
wb = Data.wb;       inFcn = Data.inFcn;  

% Find intervention time and split time array in two
[~,T_id] = min(abs(t-ts)); 
t_pre = t(1:T_id);      t_post = t(T_id+1:end); 

% PRE t_s
irfPre = k1*exp(-k2*t_pre/(1+bp)); 
CtPre = h*filter(inFcn(1:T_id),1,irfPre); 

% POST t_s
tau = t_post - ts;
Ct_ts = CtPre(end); 
a = 1 + (1-o)*bp; 
irfPost = k1*exp(-k2*tau/a); 
CtPost = h*filter(inFcn(T_id+1:end),1,irfPost) + Ct_ts*exp(-k2*tau/a); 

% TOTAL TISSUE CURVE
Ct = [CtPre(:); CtPost(:)];

% Correct for fractional blood volume 
modelCurve = (1-vb)*Ct + vb*wb(:); 
end