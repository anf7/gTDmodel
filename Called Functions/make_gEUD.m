function gEUD = make_gEUD(DVHMatrix,DoseValues,gEUDa)

if gEUDa == 0
    gEUDa = 10^-4;
end
gEUD = zeros(size(DVHMatrix,1),1);
for PatientIndex = 1:size(DVHMatrix,1)
    PatientDVH = DVHMatrix(PatientIndex,:);
    TotalDVH = sum(DVHMatrix(PatientIndex,:));
    gEUD(PatientIndex,1) = nansum((PatientDVH/TotalDVH).*DoseValues.^gEUDa)^(1/gEUDa);
end