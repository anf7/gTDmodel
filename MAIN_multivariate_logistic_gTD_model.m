function Predictor_Table = MAIN_multivariate_logistic_gTD_model(InputData,Settings,varargin)


close all
drawnow
set(0,'DefaultFigureWindowStyle','normal')
InterceptOn = true;
FigsOn = true;

Q_ValueVector = Settings.QOptMin:Settings.QOptInt:Settings.QOptMax;
Test_Q = Settings.TestQVal;
if Test_Q == 0
    Test_Q = 10^-4;
end
gEUDa = Settings.gEUDa;


%% Create Patient Data Structure

[pData, DVHMatrix, DoseValues] = make_patient_data_structure(InputData);
TumorVol = zeros(length(pData),1);
for QInd = 1:length(pData)
    TumorVol(QInd,1) = pData(QInd).TumorVolume;
end
logTumorVol = log(TumorVol);
if FigsOn
    cmap = jet(256);
    cmap(1,:) = [1,1,1];
    dispDVHMat = ceil(255*DVHMatrix/max(DVHMatrix(:)));
    imagesc(DoseValues,1:size(DVHMatrix,1),dispDVHMat);
    colormap(cmap)
    ylabel('Patient Index')
    xlabel('Dose')
    title('Patient Dose Distribution')
    drawnow
end
%%  Create parallel pool if none exists

if isempty(gcp('nocreate'))
    parpool;
end 

%%  Calculate cell viability per radiation dose for patients and reference

[ReferenceViableVolumeVector,SurvivingFractionMatrix,t_vec] = calc_cell_viability(Settings,pData,DoseValues,DVHMatrix,TumorVol);

%%  Extract user-defined predictor info from settings structure

[Outcome,UserPredictorStr,UserPredictorVar] = get_user_setting_variables(pData,Settings);
Outcome = logical(Outcome);

%%  Calculate Dxx and Vxx predictor values

if ~isempty(Settings.DxVals)
    Dx = make_DX(DVHMatrix,DoseValues,Settings.DxVals);
else
    Dx = [];
end
if ~isempty(Settings.VxVals)
    Vx = make_VX(DVHMatrix,DoseValues,Settings.VxVals);
else
    Vx = [];
end

%%  Calculate gEUD
gEUD = make_gEUD(DVHMatrix,DoseValues,gEUDa);


%% Write predictor variables as a comma seperated string

dsetstr = 'Outcome,';
if Settings.VarOnOff(1,1)
    dsetstr = strcat(dsetstr,'gTD_Q1,');
end
if Settings.VarOnOff(2,1)
    dsetstr = strcat(dsetstr,'gTD_Qtest,');
end
if Settings.VarOnOff(3,1)
    dsetstr = strcat(dsetstr,'gTD_NLLQopt,gTD_AUCQopt,');
end
if Settings.VarOnOff(4,1)
    dsetstr = strcat(dsetstr,'TumorVol,logTumorVol,');
end
for pInd = 1:length(UserPredictorVar);
    eval([UserPredictorStr{pInd},' = ','UserPredictorVar{',num2str(pInd),'};'])
    dsetstr = strcat(dsetstr,[UserPredictorStr{pInd},',']);
end
for dInd = 1:length(Settings.DxVals)
    eval(['D',num2str(Settings.DxVals(dInd)),' = Dx(:,',num2str(dInd),');'])
    dsetstr = strcat(dsetstr,['D',num2str(Settings.DxVals(dInd)),',']);
end
for vInd = 1:length(Settings.VxVals)
    eval(['V',num2str(Settings.VxVals(vInd)),' = Vx(:,',num2str(vInd),');'])
    dsetstr = strcat(dsetstr,['V',num2str(Settings.VxVals(vInd)),',']);
end
if Settings.VarOnOff(end,1)
    dsetstr = strcat(dsetstr,'gEUD,');
end
dsetstr = dsetstr(1:end-1);

%%  Calculate gTD arrays over q-value range

NormalizedLogLikelihood_functA = zeros(1,length(Q_ValueVector));
AUC_functA = zeros(1,length(Q_ValueVector));
AUCci1_functA = zeros(1,length(Q_ValueVector));
AUCci2_functA = zeros(1,length(Q_ValueVector));
maxAUC = -inf;
maxNLL = -inf;
for QInd = 1:length(Q_ValueVector)
    
    clc
    disp(['Optimizing q-value: ',num2str(100*QInd/length(Q_ValueVector),3),'%'])
    
    SetQ = Q_ValueVector(QInd);
    SetQ(SetQ == 0) = 10^-4;
    
    gTD = calculate_gTD_array(pData,DVHMatrix,SurvivingFractionMatrix,SetQ,ReferenceViableVolumeVector,t_vec,Settings.EQdose);
     
    gTD_dataset = table(Outcome,gTD);
    gTD_modelspec = 'Outcome ~ gTD';
    gTD_mdl = fitglm(gTD_dataset,gTD_modelspec,'Distribution','binomial','Link','logit','Intercept',InterceptOn);
    
%     Rsq = mdl1.Rsquared.Ordinary;
    [~,~,~,AUC] = perfcurve(Outcome,gTD_mdl.Fitted.Response,true);
    [~,Aci] = auc([Outcome,gTD_mdl.Fitted.Response],0.05);
    AUC_functA(1,QInd) = AUC;
    AUCci1_functA(1,QInd) = Aci(1);
    AUCci2_functA(1,QInd) = Aci(2);
    if AUC > maxAUC
        maxAUC = AUC;
        gTD_AUCQopt = gTD;
        AUCQopt = SetQ;
    end

    TCPlogit = gTD_mdl.Fitted.Response;
%     NLLlogit = (sum(log(TCPlogit(Outcome)))+sum(log(1-TCPlogit(~Outcome))))/length(Outcome);
    NLLlogit = gTD_mdl.LogLikelihood;
    NormalizedLogLikelihood_functA(1,QInd) = NLLlogit;
    if NLLlogit > maxNLL
        maxNLL = NLLlogit;
        gTD_NLLQopt = gTD;
        NLLQopt = SetQ;
    end
    
end
gTD_Qtest = calculate_gTD_array(pData,DVHMatrix,SurvivingFractionMatrix,...
    Test_Q,ReferenceViableVolumeVector,t_vec,Settings.EQdose);
gTD_Q1 = calculate_gTD_array(pData,DVHMatrix,SurvivingFractionMatrix,...
    1,ReferenceViableVolumeVector,t_vec,Settings.EQdose);


predictor_dataset = [];
eval(['predictor_dataset = table(',dsetstr,');'])

%%  Display shift in gTD values for optimized q-value vs. q=1

% draw_gTD_shift_hist(gTD_AUCQopt,gTD_Q1,AUCQopt,1,Outcome);
% draw_gTD_shift_hist(gTD_Qtest,gTD_Q1,Test_Q,1,Outcome);


set(0,'DefaultFigureWindowStyle','docked')
Followup = nan(length(pData),1);
for i = 1:length(pData)
    Followup(i,1) = pData(i).TimeToEvent;
end


%% Make Predictor Table
FractGrpSize = 1/2;
p_strings_table = cell(size(predictor_dataset,2)-1,1);
logistic_formula_string_table = cell(size(predictor_dataset,2)-1,1);
Rsq_table = cell(size(predictor_dataset,2)-1,1);
AUC_table = cell(size(predictor_dataset,2)-1,1);
NLL_table = cell(size(predictor_dataset,2)-1,1);
HLp_table = cell(size(predictor_dataset,2)-1,1);
KMp_table = cell(size(predictor_dataset,2)-1,1);
MTTFHR_table = cell(size(predictor_dataset,2)-1,1);
MTTFLR_table = cell(size(predictor_dataset,2)-1,1);
Sp_table = cell(size(predictor_dataset,2)-1,1);
Dev_table = cell(size(predictor_dataset,2)-1,1);
p_table = cell(size(predictor_dataset,2)-1,1);
coeff_table = cell(size(predictor_dataset,2)-1,1);

set(0,'DefaultFigureWindowStyle','docked')
bestKMp = inf;
for pInd = 2:size(predictor_dataset,2)   
    
    predictor_subset = predictor_dataset(:,[1,pInd]);
    psub_strings = predictor_subset.Properties.VariableNames;
    p_strings_table{pInd-1} = psub_strings{2};
    
    KMp = [];
    if ~any(isnan(Followup))
        predictor_array = table2array(predictor_dataset(:,pInd));
        predictor_string = p_strings_table{pInd-1};
        if strcmp(predictor_string,'TumorVol') || strcmp(predictor_string,'logTumorVol')
            [KMp,FHR,FLR] = kmplots(predictor_array,~Outcome,Followup,'PredictorName',predictor_string,...
                'AddTitle',true,'TimeUnits','Years','FractionalGroupSize',FractGrpSize,...
                'RiskCorrelation','direct','ShowFigures',FigsOn);
        else
            [KMp,FHR,FLR] = kmplots(predictor_array,~Outcome,Followup,'PredictorName',predictor_string,...
                'AddTitle',true,'TimeUnits','Years','FractionalGroupSize',FractGrpSize,'ShowFigures',FigsOn);
        end
    end
    if ~isempty(KMp) && KMp < bestKMp
        bestKMp = KMp;
    end
    if ~isempty(KMp)
        KMp_table{pInd-1} = num2str(KMp);
        MTTFHR_table{pInd-1} = num2str(FHR);
        MTTFLR_table{pInd-1} = num2str(FLR);
    else
        KMp_table{pInd-1} = 'N/A';
        MTTFHR_table{pInd-1} = 'N/A';
        MTTFLR_table{pInd-1} = 'N/A';
    end
    
    psub_modelspec = [psub_strings{1},' ~ ',psub_strings{2}];
    psub_mdl = fitglm(predictor_subset,psub_modelspec,'Distribution','binomial','Link','logit','Intercept',InterceptOn);
    if length(psub_strings{2}) >= 3 && strcmp(psub_strings{2}(1:3),'gTD')
        [HLp,outstr] = draw_TCP_curve_logitfit(psub_mdl,Outcome,FigsOn,FigsOn,Settings.EQdose,true);
    elseif strcmp(psub_strings{2},'TumorVol')
        [HLp,outstr] = draw_TCP_curve_logitfit(psub_mdl,Outcome,FigsOn,FigsOn,'Volume (cc)',false);
    elseif strcmp(psub_strings{2},'logTumorVol')
        [HLp,outstr] = draw_TCP_curve_logitfit(psub_mdl,Outcome,FigsOn,FigsOn,'log(Volume (cc))',false);
    elseif length(psub_strings{2}) > 1 && strcmpi(psub_strings{2}(1),'D') && isnumeric(str2double(psub_strings{2}(2:end)))
        [HLp,outstr] = draw_TCP_curve_logitfit(psub_mdl,Outcome,FigsOn,FigsOn,'Dose (Gy)',false);
    elseif length(psub_strings{2}) > 1 && strcmpi(psub_strings{2}(1),'V') && isnumeric(str2double(psub_strings{2}(2:end)))
        [HLp,outstr] = draw_TCP_curve_logitfit(psub_mdl,Outcome,FigsOn,FigsOn,'Volume (cc)',false);
    else
        [HLp,outstr] = draw_TCP_curve_logitfit(psub_mdl,Outcome,FigsOn,FigsOn,'Predictor Units',false);
    end

    logistic_formula_string_table{pInd-1} = outstr;
    HLp_table{pInd-1} = num2str(HLp);
    Rsq = psub_mdl.Rsquared.Ordinary;

    Rsq_table{pInd-1} = num2str(Rsq);
    [~,~,~,AUC] = perfcurve(Outcome,psub_mdl.Fitted.Response,true);

    AUC_table{pInd-1} = num2str(AUC);
    TCPlogit = psub_mdl.Fitted.Response;
%     NLL = (sum(log(TCPlogit(Outcome)))+sum(log(1-TCPlogit(~Outcome))))/length(Outcome);
    NLL = psub_mdl.LogLikelihood;
    
    NLL_table{pInd-1} = num2str(NLL);
    
    Dev_table{pInd-1} = num2str(psub_mdl.Deviance);
%     p_table{pInd-1} = num2str(psub_mdl.Coefficients{2,4});
    pval = coefTest(psub_mdl);
    p_table{pInd-1} = num2str(pval);
    coeff_table{pInd-1} = num2str(psub_mdl.Coefficients{2,1});
    Sp_table{pInd-1} = num2str(corr(Outcome,psub_mdl.Fitted.Response,'type','Spearman'));
end


if nargin > 2
    for i = 1:length(varargin)
        if ischar(varargin{i})
            p_modelspec = [psub_strings{1},' ~ ',varargin{i}];
            try
                [~] = fitglm(predictor_dataset,p_modelspec,'Distribution','binomial','Link','logit','Intercept',InterceptOn);
            catch
                error('User specified predictor(s) do not exist')
            end
            
            if ~isempty(strfind(p_modelspec,'gTD_NLLQopt')) || ~isempty(strfind(p_modelspec,'gTD_AUCQopt'))
                maxAUC = -inf;
                maxNLL = -inf;
                for QInd = 1:length(Q_ValueVector)
                    clc
                    disp(['Optimizing q-value (multivariate model): ',num2str(100*QInd/length(Q_ValueVector),3),'%'])

                    SetQ = Q_ValueVector(QInd);
                    SetQ(SetQ == 0) = 10^-4;

                    gTD = calculate_gTD_array(pData,DVHMatrix,SurvivingFractionMatrix,SetQ,ReferenceViableVolumeVector,t_vec,Settings.EQdose);
                    gTD_plus_predictor_mdl = fitglm(predictor_dataset,p_modelspec,'Distribution','binomial','Link','logit','Intercept',InterceptOn);

                    if strfind(p_modelspec,'gTD_AUCQopt')
                        [~,~,~,AUC] = perfcurve(Outcome,gTD_plus_predictor_mdl.Fitted.Response,true);
                        if AUC > maxAUC
                            maxAUC = AUC;
                            user_predictor_mdl = gTD_plus_predictor_mdl;
                        end
                    end

                    if strfind(p_modelspec,'gTD_NLLQopt')
                        TCPlogit = gTD_plus_predictor_mdl.Fitted.Response;
                        NLLlogit = (sum(log(TCPlogit(Outcome)))+sum(log(1-TCPlogit(~Outcome))))/length(Outcome);
                        if NLLlogit > maxNLL
                            maxNLL = NLLlogit;
                            user_predictor_mdl = gTD_plus_predictor_mdl;
                        end
                    end

                end
            else
                user_predictor_mdl = fitglm(predictor_dataset,p_modelspec,'Distribution','binomial','Link','logit','Intercept',InterceptOn);
            end

            newpstr = varargin{i};
            if strfind(p_modelspec,'gTD_NLLQopt')
                strind = strfind(varargin{i},'gTD_NLLQopt');
                if ~isempty(strind)
                    p_string_gTD_2 = varargin{i}(strind+11:end);
                    newpstr = [varargin{i}(1:strind+2),' (best NLL: q = ',num2str(maxNLL,2),')',p_string_gTD_2];
                end
            elseif strfind(p_modelspec,'gTD_AUCQopt')
                strind = strfind(varargin{i},'gTD_AUCQopt');
                if ~isempty(strind)
                    p_string_gTD_2 = varargin{i}(strind+11:end);
                    newpstr = [varargin{i}(1:strind+2),' (best AUC: q = ',num2str(maxAUC,2),')',p_string_gTD_2];
                end
            end
                              
            p_strings_table{end+1} = newpstr;
            [HLp,outstr] = draw_TCP_curve_logitfit(user_predictor_mdl,Outcome,FigsOn,FigsOn,'Predictor Function Value',false);
            
            logistic_formula_string_table{end+1} = outstr;
            HLp_table{end+1} = num2str(HLp);
            Rsq = user_predictor_mdl.Rsquared.Ordinary;
            Rsq_table{end+1} = num2str(Rsq);
            [~,~,~,AUC] = perfcurve(Outcome,user_predictor_mdl.Fitted.Response,true);
            AUC_table{end+1} = num2str(AUC);
            TCPlogit = user_predictor_mdl.Fitted.Response;
%             NLL = (sum(log(TCPlogit(Outcome)))+sum(log(1-TCPlogit(~Outcome))))/length(Outcome);
            NLL = user_predictor_mdl.LogLikelihood;
            NLL_table{end+1} = num2str(NLL);
            KMp_table{end+1} = 'N/A';
            MTTFHR_table{end+1} = 'N/A';
            MTTFLR_table{end+1} = 'N/A';
            Sp_table{end+1} = num2str(corr(Outcome,user_predictor_mdl.Fitted.Response,'type','Spearman'));
            Dev_table{end+1} = num2str(user_predictor_mdl.Deviance);
            pval = coefTest(user_predictor_mdl);
%             p_table{end+1} = num2str(user_predictor_mdl.Coefficients{2,4});
            p_table{end+1} = num2str(pval);
            coeff_table{end+1} = num2str(user_predictor_mdl.Coefficients{2,1});
                       
            if ~any(isnan(Followup))
                [KMp_m,FHR_m,FLR_m] = kmplots(user_predictor_mdl.Fitted.Response,~Outcome,Followup,'PredictorName',varargin{i},...
                    'AddTitle',true,'TimeUnits','Years','FractionalGroupSize',FractGrpSize);
            end
            if exist('KMp_m','var')
                KMp_table{end} = num2str(KMp_m);
                MTTFHR_table{end} = num2str(FHR_m);
                MTTFLR_table{end} = num2str(FLR_m);
            end
        end
    end
end
            
%% Full model
% fullmdl_string = [];
% for n = 1:size(predictor_dataset,2)
%     vars = predictor_dataset(:,n);
%     var_str2 = vars.Properties.VariableNames;
%     if ~isempty(strfind(var_str2{1},'gTD')) && ~strcmp(var_str2{1},'gTD_NLLQopt')
%         continue
%     end
%     if n == 1
%         fullmdl_string = strcat(fullmdl_string,var_str2{1},'~');
%     else
%         fullmdl_string = strcat(fullmdl_string,var_str2{1},'+');
%     end
% end
% fullmdl_string = fullmdl_string(1:end-1);
% 
% fullmdl = fitglm(predictor_dataset,fullmdl_string,'Distribution','binomial','Link','logit','Intercept',InterceptOn);
% 
% newpstr = 'All_predictors';
% 
% p_strings_table{end+1} = newpstr;
% [HLp,outstr] = draw_TCP_curve_logitfit(fullmdl,Outcome,TCPFigsOn,TCPFigsOn,'Predictor Function Value',false);
% 
% logistic_formula_string_table{end+1} = outstr;
% HLp_table{end+1} = num2str(HLp);
% Rsq = fullmdl.Rsquared.Ordinary;
% Rsq_table{end+1} = num2str(Rsq);
% [~,~,~,AUC] = perfcurve(Outcome,fullmdl.Fitted.Response,true);
% AUC_table{end+1} = num2str(AUC);
% TCPlogit = fullmdl.Fitted.Response;
% %             NLL = (sum(log(TCPlogit(Outcome)))+sum(log(1-TCPlogit(~Outcome))))/length(Outcome);
% NLL = fullmdl.LogLikelihood;
% NLL_table{end+1} = num2str(NLL);
% KMp_table{end+1} = 'N/A';
% Sp_table{end+1} = num2str(corr(Outcome,fullmdl.Fitted.Response,'type','Spearman'));
% Dev_table{end+1} = num2str(fullmdl.Deviance);
% pval = coefTest(fullmdl);
% %             p_table{end+1} = num2str(user_predictor_mdl.Coefficients{2,4});
% p_table{end+1} = num2str(pval);
% 
% if ~any(isnan(Followup))
%     KMp_m = kmplots(fullmdl.Fitted.Response,~Outcome,Followup/365,'PredictorName',varargin{i},...
%         'AddTitle',true,'TimeUnits','Years','FractionalGroupSize',FractGrpSize);
% end
% if exist('KMp_m','var')
%     KMp_table{end} = num2str(KMp_m);
% end




%% Make Tables


set(0,'DefaultFigureWindowStyle','normal')
cnames = {'Predictors','Logistic Formula','AUC','normalized LL','R squared','HL p-value','KM log rank p-value','Median Time To Failure High Risk','Median Time To Failure High Risk','Deviance','Coeff. p-value','Coeff. Value'};
rnames = [];
table_data = [p_strings_table,logistic_formula_string_table,AUC_table,NLL_table,Rsq_table,HLp_table,KMp_table,MTTFHR_table,MTTFLR_table,Dev_table,p_table,coeff_table];
for i = 1:size(table_data,1)
    k = strfind(table_data{i,1},'gTD_Q1');
    if ~isempty(k) 
        table_data{i,1} = [table_data{i,1}(1:k(1)-1),'gTD (q = 1)',table_data{i,1}(k(1)+6:end)];
    end
    k = strfind(table_data{i,1},'gTD_Qtest');
    if ~isempty(k) 
        table_data{i,1} = [table_data{i,1}(1:k(1)-1),'gTD (q = ',num2str(Test_Q),')',table_data{i,1}(k(1)+9:end)];
    end
    k = strfind(table_data{i,1},'gTD_NLLQopt');
    if ~isempty(k)
        table_data{i,1} = [table_data{i,1}(1:k(1)-1),'gTD (best NLL: q = ',num2str(NLLQopt),')',table_data{i,1}(k(1)+11:end)];
    end
    k = strfind(table_data{i,1},'gTD_AUCQopt');
    if ~isempty(k)
        table_data{i,1} = [table_data{i,1}(1:k(1)-1),'gTD (best AUC: q = ',num2str(AUCQopt),')',table_data{i,1}(k(1)+11:end)];
    end
    k = strfind(table_data{i,1},'TumorVol');
    if ~isempty(k)
        table_data{i,1} = [table_data{i,1}(1:k(1)-1),'Tumor Volume',table_data{i,1}(k(1)+8:end)];
    end
    k = strfind(table_data{i,1},'logTumorVol');
    if ~isempty(k)
        table_data{i,1} = [table_data{i,1}(1:k(1)-1),'log(Tumor Volume)',table_data{i,1}(k(1)+11:end)];
    end
end
table_data_nohtml = table_data;
for i = 1:size(table_data,1)
    for j = 1:size(table_data,2)
        table_data{i,j} = ['<html><font color="black">',table_data{i,j},'</body></html>'];
    end
end
% table_data{bestAUCInd,1} = ['<html><font color="black"><b>',table_data_nohtml{bestAUCInd,1},'</b></body></html>'];
% table_data{bestAUCInd,2} = ['<html><font color="black"><b>',table_data_nohtml{bestAUCInd,2},'</b></body></html>'];
% table_data{bestAUCInd,3} = ['<html><font color="black"><b>',table_data_nohtml{bestAUCInd,3},'</b></body></html>'];
% 
% table_data{bestNLLInd,1} = ['<html><font color="black"><b>',table_data_nohtml{bestNLLInd,1},'</b></body></html>'];
% table_data{bestNLLInd,2} = ['<html><font color="black"><b>',table_data_nohtml{bestNLLInd,2},'</b></body></html>'];
% table_data{bestNLLInd,4} = ['<html><font color="black"><b>',table_data_nohtml{bestNLLInd,4},'</b></body></html>'];
% 
% table_data{bestRsqInd,1} = ['<html><font color="black"><b>',table_data_nohtml{bestRsqInd,1},'</b></body></html>'];
% table_data{bestRsqInd,2} = ['<html><font color="black"><b>',table_data_nohtml{bestRsqInd,2},'</b></body></html>'];
% table_data{bestRsqInd,5} = ['<html><font color="black"><b>',table_data_nohtml{bestRsqInd,5},'</b></body></html>'];
% 
% table_data{bestHLpInd,1} = ['<html><font color="black"><b>',table_data_nohtml{bestHLpInd,1},'</b></body></html>'];
% table_data{bestHLpInd,2} = ['<html><font color="black"><b>',table_data_nohtml{bestHLpInd,2},'</b></body></html>'];
% table_data{bestHLpInd,6} = ['<html><font color="black"><b>',table_data_nohtml{bestHLpInd,6},'</b></body></html>'];
% 
% if ~isempty(bestKMpInd)
%     table_data{bestKMpInd,1} = ['<html><font color="black"><b>',table_data_nohtml{bestKMpInd,1},'</b></body></html>'];
%     table_data{bestKMpInd,2} = ['<html><font color="black"><b>',table_data_nohtml{bestKMpInd,2},'</b></body></html>'];
%     table_data{bestKMpInd,7} = ['<html><font color="black"><b>',table_data_nohtml{bestKMpInd,7},'</b></body></html>'];
% end

[AUC_table_out,NLL_table_out,Rsq_table_out,HLp_table_out,KMp_table_out,MTTFHR_table_out,MTTFLR_table_out,Dev_table_out,p_table_out,coeff_table_out]...
    = deal(nan(size(table_data,1),1));
for i = 1:size(table_data,1)
    AUC_table_out(i,1) = str2double(AUC_table{i});
    NLL_table_out(i,1) = str2double(NLL_table{i});
    Rsq_table_out(i,1) = str2double(Rsq_table{i});
    HLp_table_out(i,1) = str2double(HLp_table{i});
    Sp_table_out(i,1) = str2double(Sp_table{i});
    Dev_table_out(i,1) = str2double(Dev_table{i});
    p_table_out(i,1) = str2double(p_table{i});
    coeff_table_out(i,1) = str2double(coeff_table{i});

    KMnum = str2double(KMp_table{i});
    FHRnum = str2double(MTTFHR_table{i});
    FLRnum = str2double(MTTFLR_table{i});
    if isempty(KMnum)
        KMp_table_out(i,1) = NaN;
        MTTFHR_table_out(i,1) = NaN;
        MTTFLR_table_out(i,1) = NaN;
    else
        KMp_table_out(i,1) = KMnum;
        MTTFHR_table_out(i,1) = FHRnum;
        MTTFLR_table_out(i,1) = FLRnum;
    end
end
Predictor_Table = table(table_data_nohtml(:,1),table_data_nohtml(:,2),AUC_table_out,NLL_table_out,Rsq_table_out,HLp_table_out,KMp_table_out,MTTFHR_table_out,MTTFLR_table_out,Dev_table_out,p_table_out,coeff_table_out,...
    'VariableNames',{'Predictor','Logistic_Formula','AUC','norm_LL','R_squared','Hosmer_Lemeshow_p_Value','KM_LogRank_p_Value','Median_Time_To_Failure_High_Risk','Median_Time_To_Failure_Low_Risk','Deviance','Coeff_p_value','Coeff'});

if FigsOn
    
    f = figure;

    t = uitable(f,'Data',table_data,'ColumnName',cnames,'RowName',rnames);

    set(t,'Position',[20, 20, 1164, (size(table_data,1) + 1)*19])
    set(f,'Position',[27, 50, 1205, (size(table_data,1) + 1)*19 + 40]);
    jscrollpane = findjobj_u(t);
    jtable = jscrollpane.getViewport.getView;
    jtable.setSortable(true);
    jtable.setAutoResort(true);
    jtable.setMultiColumnSortable(true);
    jtable.setPreserveSelectionsAfterSorting(true);
    set(t,'ColumnWidth',{300,360,90,90,90,90,140})
    set(t,'RowStriping','off')
    set(f,'Name','Predictor Performance Table (sortable)')


    
    %% Plot AUC vs. q-value

    corder = [0,0.447,0.741;...
        0.850,0.325,0.098;
        0.929,0.694,0.125;
        0.494,0.184,0.556;
        0.466,0.674,0.188;
        0.301,0.745,0.933;
        0.635,0.078,0.184;
        1.000,0.000,0.000;
        0.000,1.000,0.000;
        0.000,0.000,1.000;
        1.000,1.000,0.000;
        0.000,1.000,1.000;
        1.000,0.000,1.000;
        ];

    h = figure;
    lnhgtd = plot(Q_ValueVector,AUC_functA,'-k','linewidth',2);
    hold on
    plot(Q_ValueVector,AUCci1_functA,'k:','linewidth',2)
    plot(Q_ValueVector,AUCci2_functA,'k:','linewidth',2)

    [~,mi] = max(AUC_functA);
    lnh = [];
    legendstr = cell(0);
    Predictor_Val = [];
    for i = 1:size(Predictor_Table,1)
        if length(table_data_nohtml{i,1}) < 3 || ~strcmp(table_data_nohtml{i,1}(1:3),'gTD')
            lnh(end+1) = plot([min(Q_ValueVector),max(Q_ValueVector)],[Predictor_Table{i,3},Predictor_Table{i,3}],'-','LineWidth',2,'Color',corder(length(legendstr)+1,:));
            legendstr{end+1} = table_data_nohtml{i,1};
            Predictor_Val(end+1,1) = Predictor_Table{i,3};
        end
    end

    [predictorsort,psortind] = sort(Predictor_Val,1,'descend');
    lhnsort = [lnhgtd, lnh(psortind)];
    legenstrsort = [{'gTD'},legendstr(psortind)];


    zz = get(gca,'ylim');

    plot([Q_ValueVector(mi),Q_ValueVector(mi)],[max(AUC_functA),-10000],'k--','linewidth',1)
    patch([Q_ValueVector,fliplr(Q_ValueVector)],[AUCci1_functA,fliplr(AUCci2_functA)],[1,0,0],'FaceAlpha',0.15,'LineWidth',0.2)

    xlabel('q')
    ylabel('AUC')
    lh = legend(lhnsort,legenstrsort);

    for i = 1:length(predictorsort)
        ii = length(predictorsort) - i + 1;
        if sum(predictorsort == predictorsort(ii)) > 1
            set(lhnsort(ii + 1),'LineStyle','--')
            predictorsort(ii) = nan;
        end
    end

    customfigset(h,18)
    ylim(zz)
    set(lh,'fontsize',12)
    hold off
    set(h,'Position',[197 262 669 536])


    %% Plot nLL vs. q-value


    [~,MaxAIndex] = max(NormalizedLogLikelihood_functA);

    CI95 = false(size(Q_ValueVector));
    % CI95_nLL = (length(Outcome)*NormalizedLogLikelihood_functA(1,MaxAIndex) - 1.92)/length(Outcome);
    CI95_nLL = (NormalizedLogLikelihood_functA(1,MaxAIndex) - 1.92);
    CI99 = false(size(Q_ValueVector));
    % CI99_nLL = (length(Outcome)*NormalizedLogLikelihood_functA(1,MaxAIndex) - 3)/length(Outcome);
    CI99_nLL = (NormalizedLogLikelihood_functA(1,MaxAIndex) - 3);

    for n = 1:length(Q_ValueVector)
    %     if (length(Outcome)*NormalizedLogLikelihood_functA(1,n)) >=...
    %             (length(Outcome)*NormalizedLogLikelihood_functA(1,MaxAIndex) - 1.92)
        if (NormalizedLogLikelihood_functA(1,n)) >=...
                    (NormalizedLogLikelihood_functA(1,MaxAIndex) - 1.92)
            CI95(1,n) = true;
        else
            CI95(1,n) = false;
        end
    %     if (length(Outcome)*NormalizedLogLikelihood_functA(1,n)) >=...
    %             (length(Outcome)*NormalizedLogLikelihood_functA(1,MaxAIndex) - 3)
        if (NormalizedLogLikelihood_functA(1,n)) >=...
                    (NormalizedLogLikelihood_functA(1,MaxAIndex) - 3)
            CI99(1,n) = true;
        else
            CI99(1,n) = false;
        end
    end


    fHandle = figure;


    CI95_x = Q_ValueVector(CI95);
    CI95_y = CI95_nLL*ones(size(CI95_x));
    nLL95_y = NormalizedLogLikelihood_functA(CI95);

    CI99_x = Q_ValueVector(CI99);
    CI99_y = CI99_nLL*ones(size(CI99_x));
    nLL99_y = NormalizedLogLikelihood_functA(CI99);

    lnhgtd = plot(Q_ValueVector,NormalizedLogLikelihood_functA,'k-','linewidth',2);
    hold on
    plot([Q_ValueVector(1),Q_ValueVector(end)],[CI95_nLL,CI95_nLL],'k:','linewidth',1);
    plot([Q_ValueVector(1),Q_ValueVector(end)],[CI99_nLL,CI99_nLL],'k-.','linewidth',0.5);


    lnh = [];
    legendstr = cell(0);
    Predictor_Val = [];
    for i = 1:size(Predictor_Table,1)
        if length(table_data_nohtml{i,1}) < 3 || ~strcmp(table_data_nohtml{i,1}(1:3),'gTD')
            lnh(end+1) = plot([min(Q_ValueVector),max(Q_ValueVector)],[Predictor_Table{i,4},Predictor_Table{i,4}],...
                '-','LineWidth',2,'Color',corder(length(legendstr)+1,:));
            legendstr{end+1} = table_data_nohtml{i,1};
            Predictor_Val(end+1,1) = Predictor_Table{i,4};
        end
    end

    [predictorsort,psortind] = sort(Predictor_Val,1,'descend');
    lhnsort = [lnhgtd, lnh(psortind)];
    legenstrsort = [{'gTD'},legendstr(psortind)];


    zz = get(gca,'ylim');

    patch([CI95_x,fliplr(CI95_x)],[CI95_y,fliplr(nLL95_y)],[1,0,0],'FaceAlpha',0.25,'LineStyle','none')
    patch([CI99_x,fliplr(CI99_x)],[CI99_y,fliplr(nLL99_y)],[1,0,0],'FaceAlpha',0.15,'LineStyle','none')

    if ~CI95(1)
        plot([CI95_x(1),CI95_x(1)],[nLL95_y(1),-10000],'k:','linewidth',1)
    end
    if ~CI95(end)
        plot([CI95_x(end),CI95_x(end)],[nLL95_y(end),-10000],'k:','linewidth',1)
    end

    if ~CI99(1)
        plot([CI99_x(1),CI99_x(1)],[nLL99_y(1),-10000],'k-.','linewidth',0.5)
    end
    if ~CI99(end)
        plot([CI99_x(end),CI99_x(end)],[nLL99_y(end),-10000],'k-.','linewidth',0.5)
    end


    maxx = Q_ValueVector(MaxAIndex);
    maxy = NormalizedLogLikelihood_functA(MaxAIndex);
    plot([maxx,maxx],[maxy,-10000],'k--','linewidth',1)

    xlabel('q')
    ylabel('Normalized Log Likelihood')
    ylim(zz)
    grid on
    xlim([min(Q_ValueVector),max(Q_ValueVector)])
    ax = gca;
    set(ax,'xtick',min(Q_ValueVector):1:max(Q_ValueVector))
    customfigset(fHandle,18)

    lh = legend(lhnsort,legenstrsort);
    set(lh,'fontsize',12)

    for i = 1:length(predictorsort)
        ii = length(predictorsort) - i + 1;
        if sum(predictorsort == predictorsort(ii)) > 1
            set(lhnsort(ii + 1),'LineStyle','--')
            predictorsort(ii) = nan;
        end
    end

    xscl = max(Q_ValueVector) - min(Q_ValueVector);
    yscl = zz(2) - zz(1);

    text((min(Q_ValueVector) + 0.15*(xscl)),...
        CI95_nLL + 0.01*yscl,'Log-ratio 95% Conf. Int.','FontSize',12,'VerticalAlignment','baseline')
    text((min(Q_ValueVector) + 0.25*(xscl)),...
        CI99_nLL - 0.01*yscl,'Log-ratio 99% Conf. Int.','FontSize',12,'VerticalAlignment','top')
    set(fHandle,'Position',[197 262 669 536])
    hold off

end



