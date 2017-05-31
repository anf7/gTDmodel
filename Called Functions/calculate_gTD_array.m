function gTD = calculate_gTD_array(pData,DVHMatrix,PIH_SurvivingFraction,SetA,ReferenceViableVolumeVector,t_vec_cont,RefEQDx)
 

gTD = zeros(length(pData),1);
for PatientIndex = 1:length(pData)
   
    PatientDVH = DVHMatrix(PatientIndex,:);
    DoseBinsFractionalVolume = PatientDVH/sum(PatientDVH(:));
    DoseBinsSurvivingFraction = PIH_SurvivingFraction(PatientIndex,:);
    
    DosePrt = DoseBinsFractionalVolume.*DoseBinsSurvivingFraction.^SetA;
%     DoseBinsFractionalVolume(DoseBinsFractionalVolume == 0) = NaN;
%     DosePrt = exp(SetA*log(DoseBinsFractionalVolume.*DoseBinsSurvivingFraction));
%     DosePrt(isinf(DosePrt)) = NaN;
    
    WeightedViableVolume = ...
        pData(PatientIndex).TumorVolume*(nansum(DosePrt)).^(1/SetA);

    testvec = (ReferenceViableVolumeVector > WeightedViableVolume);
    LastTimepointIndex = find(testvec,1,'last');
    if LastTimepointIndex == length(testvec)
        error('Minimum reference viable volume exceeds viable tumor subvolume')
    end
    
    if isempty(LastTimepointIndex)
        TrtTime = 0;
    else
        TrtTime = t_vec_cont(LastTimepointIndex);
    end
    gTD(PatientIndex,1) = RefEQDx*TrtTime/24;
end
    
