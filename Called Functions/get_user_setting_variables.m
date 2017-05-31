function [Outcome,UserPredictorStr,UserPredictorVar] = get_user_setting_variables(pData,Settings)

OutcomeArray = Settings.VarOnOff(:,2);
OutcomeInd = find(OutcomeArray);
if length(OutcomeInd) ~= 1 || OutcomeInd < 5 
    error('Select one binary response array variable')
else
    OutcomeStr = Settings.VarNames{OutcomeInd};
    if ~islogical(eval(strcat('pData(1).',OutcomeStr)))
        error('Select only one binary response array variable')
    else
        Outcome = false(length(pData),1);
        for n = 1:length(Outcome)
            Outcome(n,1) = eval(strcat('pData(',num2str(n),').',OutcomeStr));
        end
    end
end


PredictorArray = Settings.VarOnOff(:,1);
UserPredictorArray = PredictorArray(5:end-1,1);
UserPredictorInd = 4 + find(UserPredictorArray);
UserPredictorVar = cell(length(UserPredictorInd),1);
UserPredictorStr = cell(length(UserPredictorInd),1);
for m = 1:length(UserPredictorInd)
    UserPredictorStr{m,1} = Settings.VarNames{UserPredictorInd(m)};
    UserPredictor = zeros(length(Outcome),1);
    for n = 1:length(UserPredictor)
        UserPredictor(n,1) = eval(strcat('pData(',num2str(n),').',UserPredictorStr{m,1}));
    end
    UserPredictorVar{m,1} = UserPredictor;
end