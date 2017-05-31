function [pData, DVHmatrix, DoseValues] = make_patient_data_structure(TableData)

DVH_DoseInterval = 1;


StructData = table2struct(TableData);

StructDataFields = fieldnames(StructData);
CellData = (struct2cell(StructData))';
SortColIndex = find(strcmp('TreatmentDuration',StructDataFields));
SortedCellData = (sortrows(CellData,SortColIndex))';
pData = cell2struct(SortedCellData, StructDataFields, 1);

NumberPatients = length(StructData);

DVH_MaxDose = 0;
for PatientIndex = 1:NumberPatients
    DVH_MaxDose = max([DVH_MaxDose,max(pData(PatientIndex).DVH(1,:))]);
end

DoseValues = 0:DVH_DoseInterval:ceil(DVH_MaxDose);

DVHmatrix = zeros(NumberPatients,length(DoseValues));
for PatientIndex = 1:NumberPatients
    for DVHDoseIndex = 1:length(pData(PatientIndex).DVH(1,:))
        Dose = pData(PatientIndex).DVH(1,DVHDoseIndex);
        FractionalVolume = pData(PatientIndex).DVH(2,DVHDoseIndex);
        DoseIndex = round(Dose*(1/DVH_DoseInterval))+1;
        DVHmatrix(PatientIndex,DoseIndex) = DVHmatrix(PatientIndex,DoseIndex) + FractionalVolume;
    end
    pData(PatientIndex).TumorVolume = sum(DVHmatrix(PatientIndex,:));
end

