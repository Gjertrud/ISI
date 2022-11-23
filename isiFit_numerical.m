function [Data] = isiFit_numerical(Data,Settings)
% Fit displacement model to in-scan intervention data numerically, using
% the Euler Forwarth method
%
%__________________________________________________________________________
%                             Gjertrud Louise Laurell & Martin Schain, 2022

% SET STEP SIZE FOR EULER FORWARD
h = 1/120; 

% UPDATE BLOOD DATA ARRAYS TO EULER FORWARD STEP SIZE (h)
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
        [models, K1, VS, vB, occ, VND, te, exitFlag] = ...
            runNumericalISI_1TCM(Data,Settings); 
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
isi.te = te;            isi.exitFlag = exitFlag; 

Data.isi = isi; 
end

% SUPPORTING FUNCTIONS - 1TCM CASE
function [models, K1, VS, vB, occ, VND, te, exitFlag] = ...
    runNumericalISI_1TCM(Data,Settings) 
% Fit data to 1TC displacement model 

% Fit outer layer parameters (occupancy, VND & te)
[pOccVndTe, ~, ~, exitFlag] = ...
    lsqnonlin(@fitEF_1TC_outer,Settings.fitParams.startParamsGlobal,...
    Settings.fitParams.lowBoundGlobal,Settings.fitParams.upBoundGlobal,...
    Settings.fitParams.options,Data,Settings);

% Pre-allocate model curves and inner parameter array
models = zeros(size(Data.TACs));
innerPars = zeros(size(Data.TACs,2),3); 

% Fit inner layer parameters with outer layer fixed
for roi = 1:size(Data.TACs,2)
    [pK1VsvB, ~, ~, ~] = ...
        lsqnonlin(@fitEF_1TC_inner,Settings.fitParams.startParamsROI,...
        Settings.fitParams.lowBoundROI,Settings.fitParams.upBoundROI,...
        Settings.fitParams.options,Data,roi,pOccVndTe);
    pTot = [pK1VsvB pOccVndTe];
    modelCurve = EF_1TCM_createCurves(pTot,Data); 
    models(:,roi) = interp1(Data.t,modelCurve(:),Data.tPET); 
    innerPars(roi,:) = pK1VsvB; 
end

% Assign out-parameters
K1 = innerPars(:,1);    VS = innerPars(:,2);    vB = innerPars(:,3); 
occ = pOccVndTe(1);     VND = pOccVndTe(2);     te = pOccVndTe(3); 
end

function [err] = fitEF_1TC_outer(currP_global,Data,Settings)
% Outer layer for 1TCM - fits global parameters (occupancy, VND & te)
err = zeros(size(Data.TACs)); 
for roi = 1:size(Data.TACs,2)
    [pK1VsvB, ~, ~, ~] = ...
        lsqnonlin(@fitEF_1TC_inner,Settings.fitParams.startParamsROI,...
        Settings.fitParams.lowBoundROI,Settings.fitParams.upBoundROI,...
        Settings.fitParams.options,Data,roi,currP_global);
    P = [pK1VsvB currP_global];
    modelCurve = EF_1TCM_createCurves(P,Data); 
    modelTAC = interp1(Data.t,modelCurve(:),Data.tPET);
    err(:,roi) = (modelTAC - Data.TACs(:,roi)).*Data.weights(:); 
end
end

function [err] = fitEF_1TC_inner(p,Data,currROI,p_global)
% Inner layer for 1TCM - fits ROI-specific parameters (K1, VS and vB)
P = [p p_global]; 
modelCurve = EF_1TCM_createCurves(P,Data); 
modelTAC = interp1(Data.t,modelCurve(:),Data.tPET);
err = (modelTAC - Data.TACs(:,currROI)).*Data.weights(:); 
end

function [modelCurve] = EF_1TCM_createCurves(P,Data)
% Draw the model curves for the current choice of parameter values using
% Euler Forward

% Set step size
h = 1/120; 

% Current parameter values
k1 = P(1);          vs0 = P(2);         vb = P(3); 
o = P(4);           vnd = P(5);         te = Data.tb + P(6); 
tb = Data.tb;       tm = mean([tb te]);  

% Blood data
t = Data.t;         wb = Data.wb;       cp = Data.inFcn; 

% Get occupancy curve
o_curve = getAgriculturalOccupancyCurve(t,tb,tm,te,o); 

% Apply effect of occupancy on VS
vs = vs0*(1-o_curve); 

% Pre-allocate tissue curve
Ct = zeros(size(t)); 

% Draw the tissue curve
for i = 2:length(t)
    Ct(i) = h*k1*cp(i-1)+(1-h*k1/(vnd+vs(i-1)))*Ct(i-1);
end

% Correct for fractional blood volume
modelCurve = (1-vb)*Ct(:) + vb*wb(:); 
end

% SUPPORTING FUNCTION - OCCUPANCY CURVE
function [o_curve] = getAgriculturalOccupancyCurve(t,tb,tm,te,o)
% This function generates an occupancy growth model based on the input
% parameters.
% t: time vector, tb: begin time of growth, tm: time of max derivative,
% te: end time of growth, o: maximal occupancy reached during the scan

agrifun = (1 + (te-t)./(te-tm)).*((t-tb)./(te-tb)).^((te-tb)/(te-tm));
f = stepfun(t,tb).*agrifun;         
y = f + stepfun(t,te).*(1 - f);               
o_curve = o*y;
end