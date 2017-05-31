function Vx = make_VX(DVHMatrix,DoseValues,X)

Vx = zeros(size(DVHMatrix,1),length(X));
for Xind = 1:length(X)
    x = X(Xind);
    [~,VXind] = min(abs(DoseValues - x));
    for PatientIndex = 1:size(DVHMatrix,1)
        PatientDVH = DVHMatrix(PatientIndex,:);
        TotalDVH = sum(DVHMatrix(PatientIndex,:));
        Vx(PatientIndex,Xind) = sum(PatientDVH(VXind:end))/TotalDVH;
    end
end