WU_1_vol = zeros(1,0);
WU_0_vol = zeros(1,0);
NKI_1_vol = zeros(1,0);
NKI_0_vol = zeros(1,0);

for n = 1:size(NSCLC_Data,1)
    if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
        WU_1_vol(1,end + 1) = sum(NSCLC_Data{n,1}{1}(2,:));
    end
end

for n = 1:size(NSCLC_Data,1)
    if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
        WU_0_vol(1,end + 1) = sum(NSCLC_Data{n,1}{1}(2,:));
    end
end

for n = 1:size(NSCLC_Data,1)
    if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
        NKI_1_vol(1,end + 1) = sum(NSCLC_Data{n,1}{1}(2,:));
    end
end

for n = 1:size(NSCLC_Data,1)
    if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
        NKI_0_vol(1,end + 1) = sum(NSCLC_Data{n,1}{1}(2,:));
    end
end

% med_0 = median(NKI_0_vol)
% med_1 = median(NKI_1_vol)
% min_0 = min(NKI_0_vol)
% max_0 = max(NKI_0_vol)
% min_1 = min(NKI_1_vol)
% max_1 = max(NKI_1_vol)
% [~,p_vol] = ttest2(NKI_1_vol,NKI_0_vol)


% WU_1_md = zeros(1,0);
% WU_0_md = zeros(1,0);
% NKI_1_md = zeros(1,0);
% NKI_0_md = zeros(1,0);
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
%         dvh = NSCLC_Data{n,1}{1};
%         WU_1_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:));  
%     end
% end
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
%         dvh = NSCLC_Data{n,1}{1};
%         WU_0_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:));  
%     end
% end
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
%         dvh = NSCLC_Data{n,1}{1};
%         NKI_1_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:));
%     end
% end
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
%         dvh = NSCLC_Data{n,1}{1};
%         NKI_0_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:));
%     end
% end
% 
% med_1 = median(WU_1_md)
% min_1 = min(WU_1_md)
% max_1 = max(WU_1_md)
% 
% med_0 = median(WU_0_md)
% min_0 = min(WU_0_md)
% max_0 = max(WU_0_md)
% 
% 
% [~,p_md] = ttest2(WU_1_md,WU_0_md)




% WU_1_md = zeros(1,0);
% WU_0_md = zeros(1,0);
% NKI_1_md = zeros(1,0);
% NKI_0_md = zeros(1,0);
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
%         dvh = NSCLC_Data{n,1}{1};
%         WU_1_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:))/NSCLC_Data{n,6};  
%     end
% end
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
%         dvh = NSCLC_Data{n,1}{1};
%         WU_0_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:))/NSCLC_Data{n,6};  
%     end
% end
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
%         dvh = NSCLC_Data{n,1}{1};
%         NKI_1_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:))/NSCLC_Data{n,6};
%     end
% end
% 
% for n = 1:size(NSCLC_Data,1)
%     if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
%         dvh = NSCLC_Data{n,1}{1};
%         NKI_0_md(1,end + 1) = sum(dvh(1,:).*dvh(2,:))/sum(dvh(2,:))/NSCLC_Data{n,6};
%     end
% end
% 
% med_1 = median([WU_1_md,NKI_1_md])
% min_1 = min([WU_1_md,NKI_1_md])
% max_1 = max([WU_1_md,NKI_1_md])
% 
% med_0 = median([WU_0_md,NKI_0_md])
% min_0 = min([WU_0_md,NKI_0_md])
% max_0 = max([WU_0_md,NKI_0_md])
% 
% 
% [~,p_md] = ttest2([WU_1_md,NKI_1_md],[WU_0_md,NKI_0_md])


WU_1_f = zeros(1,0);
WU_0_f = zeros(1,0);
NKI_1_f = zeros(1,0);
NKI_0_f = zeros(1,0);

for n = 1:size(NSCLC_Data,1)
    if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
        WU_1_f(1,end + 1) = NSCLC_Data{n,5};  
    end
end

for n = 1:size(NSCLC_Data,1)
    if strcmp('WUSTL',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
        WU_0_f(1,end + 1) = NSCLC_Data{n,5};  
    end
end

for n = 1:size(NSCLC_Data,1)
    if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 1
        NKI_1_f(1,end + 1) = NSCLC_Data{n,5};
    end
end

for n = 1:size(NSCLC_Data,1)
    if strcmp('NKI',NSCLC_Data{n,2}) && NSCLC_Data{n,3} == 0
        NKI_0_f(1,end + 1) = NSCLC_Data{n,5};
    end
end

med_1 = median(WU_1_f)
min_1 = min(WU_1_f)
max_1 = max(WU_1_f)

med_0 = median(WU_0_f)
min_0 = min(WU_0_f)
max_0 = max(WU_0_f)


[~,p_f] = ttest2(WU_1_f,WU_0_f)