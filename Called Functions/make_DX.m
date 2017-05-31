function Dx = make_DX(DVHMatrix,DoseValues,X)


Dx = zeros(size(DVHMatrix,1),length(X));
for Xind = 1:length(X)
    x = X(Xind);
    for PatientIndex = 1:size(DVHMatrix,1)
        PatientDVH = DVHMatrix(PatientIndex,:);
        TotalDVH = sum(DVHMatrix(PatientIndex,:));
        for Dind = 1:length(DoseValues);
            if sum(PatientDVH(1:Dind)) <= TotalDVH*((100-x)/100)
                Dx(PatientIndex,Xind) = DoseValues(Dind);
            end
        end
    end
end