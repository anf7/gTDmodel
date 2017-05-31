function [ReferenceViableVolumeVector,SurvivingFractionMatrix,t_vec] = calc_cell_viability(Settings,pData,DoseValues,DVHMatrix,TumorVol)



SameTreatmentSchedule = zeros(length(pData),1);
TreatmentScheduleCell = cell(length(pData),1);
NumberFractionsCell = cell(length(pData),1);
for PatientIndex = 1:length(pData) 
    if PatientIndex == 1
        SameTreatmentSchedule(PatientIndex) = 1;
    elseif isequal(pData(PatientIndex).TreatmentSchedule,pData(PatientIndex - 1).TreatmentSchedule)
        SameTreatmentSchedule(PatientIndex) = SameTreatmentSchedule(PatientIndex - 1);
    else
        SameTreatmentSchedule(PatientIndex) = PatientIndex;
    end
    TreatmentScheduleCell{PatientIndex,1} = pData(PatientIndex).TreatmentSchedule;
    NumberFractionsCell{PatientIndex,1} = pData(PatientIndex).NumberFractions;
end
    
UniqueIndices = unique(SameTreatmentSchedule);

UniqueEmptyDVH = false(length(UniqueIndices),length(DoseValues));
for CaseIndex = 1:length(UniqueIndices) 
    RedundancyIndices = SameTreatmentSchedule == UniqueIndices(CaseIndex);
    UniqueEmptyDVH(CaseIndex,:) = ~logical(sum(DVHMatrix(RedundancyIndices,:),1));
end

UniqueSurvivingFractionMatrix = zeros(length(UniqueIndices),length(DoseValues));

UniqeTreatmentScheduleCell = TreatmentScheduleCell(UniqueIndices,1);
UniqeNumberFractionsCell = NumberFractionsCell(UniqueIndices,1);
AlphaVal = Settings.Alpha;
AlphaBetaVal = Settings.AlphaBeta;
parfor CaseIndex = 1:length(UniqueIndices)  
    EmptyDVHbin = UniqueEmptyDVH(CaseIndex,:);
    TreatmentSchedule = find(UniqeTreatmentScheduleCell{CaseIndex,1});
    DosePerFraction = DoseValues/UniqeNumberFractionsCell{CaseIndex,1};
    [SurvivingFraction,~,~,~] = surviving_fraction_model(DosePerFraction,...
        TreatmentSchedule,AlphaVal,AlphaBetaVal,EmptyDVHbin);    
    UniqueSurvivingFractionMatrix(CaseIndex,:) = SurvivingFraction;
end

SurvivingFractionMatrix = zeros(length(pData),length(DoseValues));
for PatientIndex = 1:length(pData) 
    UniqueIndex = SameTreatmentSchedule(PatientIndex);
    SurvivingFractionMatrix(PatientIndex,:) = UniqueSurvivingFractionMatrix(UniqueIndices == UniqueIndex,:);
end

NanSurvivingFractionMatrix = SurvivingFractionMatrix;
NanSurvivingFractionMatrix(NanSurvivingFractionMatrix == 0) = NaN;
MinSurvivingFraction = 0.1*(min(TumorVol)/Settings.ReferenceVolume)*nanmin(NanSurvivingFractionMatrix(:));



[~,SurvivingFraction_t,t_vec,trt_pts] = surviving_fraction_model(Settings.EQdose,1,AlphaVal,AlphaBetaVal,0,MinSurvivingFraction);
Y = SurvivingFraction_t(logical(trt_pts));
X = t_vec(logical(trt_pts));
X = X(logical(Y));
Y = Y(logical(Y));
ii = find(trt_pts,length(X),'first');
t_vec = t_vec(1:ii(end));
Y_interp = exp(interp1(X,log(Y),t_vec));
ReferenceViableVolumeVector = Settings.ReferenceVolume*Y_interp;

% SurvivingFractionMatrix = log(SurvivingFractionMatrix);