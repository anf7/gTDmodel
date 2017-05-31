function [HLp,outstr] = draw_TCP_curve_logitfit(logisticmodel,Outcome,showHL,showTCP,RefUnits,showJJ)


if ~iscell(logisticmodel)
    mdl = {logisticmodel};
else
    mdl = logisticmodel;
end


xpts = mdl{1}.Fitted.LinearPredictor;
ypts = mdl{1}.Fitted.Response;

showtextdata = true;

bins = 5;
binoverlap = 1;


C = table2cell(mdl{1}.Coefficients);

if size(C,1) > 1
    c0 = C{1,1};
    c1 = C{2,1};
else
    c0 = 0;
    c1 = C{1,1};
end

xpts0 = (xpts - c0)/c1;


[sX,sortind] = sort(xpts);
% [sX0,~] = sort(xpts0);
sX0 = xpts0(sortind);
sO = Outcome(sortind);
sY = ypts(sortind);

% OLthresh = 2*std(sX);
OLthresh = inf;

sX_cut_ind = find((sX < mean(sX) + OLthresh) & (sX > mean(sX) - OLthresh));    
sX_cut_ind0 = find((sX0 < mean(sX0) + OLthresh) & (sX0 > mean(sX0) - OLthresh));
sO_cut = sO(sX_cut_ind);   
sX_cut = sX(sX_cut_ind);
sX0_cut = sX0(sX_cut_ind0);

rStart = round(linspace(1,length(Outcome)-binoverlap*length(Outcome)/bins,bins));
rEnd = round(linspace(1+binoverlap*length(Outcome)/bins,length(Outcome),bins));


for n = 1:bins
    if binoverlap == 1
        r1 = round(length(sO_cut)*((n-1)/bins))+1;
        r2 = round(length(sO_cut)*((n)/bins));
%         r10 = round(length(sX0_cut)*((n-1)/bins))+1;
%         r20 = round(length(sX0_cut)*((n)/bins));
    else
        r1 = rStart(n);
        r2 = rEnd(n);
%         r10 = rStart(n);
%         r20 = rEnd(n);
    end
    Obin{n} = sO_cut(r1:r2);
    Xbin{n} = sX_cut(r1:r2);
    X0bin{n} = sX0_cut(r1:r2);
end
joint(1) = min(X0bin{1});
for n = 1:bins-1
    joint(n+1) = 0.5*(max(X0bin{n})+min(X0bin{n+1}));
end
joint(bins+1) = max(X0bin{bins});

for n = 1:bins
    meanO(n) = mean(Obin{n});
    medianX(n) = median(Xbin{n});
    predictedO(n) = 1/(1 + exp(-medianX(n)));
    binsize(n) = numel(Obin{n});
end


sX00 = sX0(~sO);
sX01 = sX0(sO);
sY0 = sY(~sO);
sY1 = sY(sO);

for n = 1:length(sX00)
    ysO0(n) = sY0(n) - 0.005 - 0.03*rand(1);
end

ysO0(ysO0 < 0.005) = 0.005;
for n = 1:length(sX01)
    ysO1(n) = sY1(n) + 0.005 + 0.03*rand(1);
end
ysO1(ysO1 > 0.995) = 0.995;

rangeX = max(sX0_cut) - min(sX0_cut);

X0fit = linspace(min(sX0_cut)-0.15*rangeX,max(sX0_cut)+0.15*rangeX,1000);
Yfit = 1./(1 + exp(-(c1*X0fit + c0)));

addYfit = cell(length(mdl)-1);
for n = 2:length(mdl)
    a1C = table2cell(mdl{n}.Coefficients);
    a1c0 = a1C{1,1};
    a1c1 = a1C{2,1};
    addYfit{n-1} = 1./(1 + exp(-(a1c1*X0fit + a1c0)));
end

%% Hosmer-Lemeshaw Test

    
Z = nan(1,bins);
PIg = nan(1,bins);
Orate = nan(1,bins);
for n = 1:bins
    Og = sum(Obin{n});
    Ng = numel(Obin{n});
    PIg(1,n) = mean(1./(1 + exp(-Xbin{n})));
    Orate(1,n) = Og/Ng;
    
    Z(1,n) = ((Og - Ng*PIg(1,n))^2)/(Ng*PIg(1,n)*(1-PIg(1,n)));    
end

ChiSqHL = sum(Z);
HLp = chi2cdf(ChiSqHL,bins-2,'upper');



eq = [mdl{1}.Formula.LinearPredictor,' '];
eqb = eq + 0;
eqb(eqb == 42) = 58;
eq = char(eqb);
T{1} = num2str(mdl{1}.Coefficients{1,1});
count = 2;
for n = 1:length(mdl{1}.CoefficientNames)
    ckstr = [' ',mdl{1}.CoefficientNames{1,n},' '];
    if strfind(eq,ckstr);
        T{count} = [num2str(mdl{1}.Coefficients{n,1}),'*',mdl{1}.CoefficientNames{1,n}];
        count = count + 1;
    end
end
TT{1} = T{1};
for n = 2:length(T)
    Tq = T{n} + 0;
    Tq(Tq == 58) = 42;
    if Tq(1) == 45
        Tq(4:end+2) = Tq(2:end);
        Tq(1:3) = [32,45,32];
    else
        Tq(4:end+3) = Tq(1:end);
        Tq(1:3) = [32,43,32];
    end
    TT{n} = char(Tq);
end
outstr = ['Y = '];
for n = 1:length(TT)
    outstr = [outstr,TT{n}];
end



if showHL
    %% Draw Calibration Plot
    
    for n = 1:bins
        [d, e] = binofit(round(meanO(n)*binsize(n)),binsize(n));
        mv(n) = d;
        ci1(n) = mv(n) - e(1);
        ci2(n) = e(2) - mv(n);
    end

    hHL = figure;
%     plot(PIg,Orate,'ok','markersize',15);
    errorbar(PIg,Orate,ci1,ci2,'ok','markersize',15);
    
    hold on
    plot([-0.05,1.05],[-0.05,1.05],'--k','linewidth',2)
%     xlim([max([0,floor(10*(min([min(meanO),min(predictedO)])-0.025))/10]),...   
%         min([1,ceil(10*(max([max(meanO),max(predictedO)])+0.025))/10])])
%     ylim([max([0,floor(10*(min([min(meanO),min(predictedO)])-0.025))/10]),...
%         min([1,ceil(10*(max([max(meanO),max(predictedO)])+0.025))/10])])
    
    xlim([min([min(PIg),min(Orate - ci1)])-0.05,max([max(PIg),max(Orate + ci2)])+0.05])
    ylim([min([min(PIg),min(Orate - ci1)])-0.05,max([max(PIg),max(Orate + ci2)])+0.05])
    
    ylabel('Observed Control Rate')
    xlabel('Predicted Control Rate')
    title(['Hosmer-Lemeshow Test p = ',num2str(HLp)])
    customfigset(hHL,15)
    set(hHL,'Name',outstr)
end
    
if showTCP 
    %% Draw TCP Figure

    h = figure;
    drawnow
    for n = 1:bins
        ph = patch([joint(n),joint(n),joint(n+1),joint(n+1)],...
            [0,1,1,0],[1,1,1],'linewidth',2);
        set(ph,'EdgeColor',[0.7,0.7,0.7]);
        t(n) = text((joint(n)+joint(n+1))/2,0.05,['Q_',num2str(n)]);
        set(t(n),'VerticalAlignment','baseline','HorizontalAlignment','center','Color',[0.7,0.7,0.7])
    end
    hold on
    for n = 1:bins
        [d, e] = binofit(round(meanO(n)*binsize(n)),binsize(n));
        mv(n) = d;
        ci1(n) = mv(n) - e(1);
        ci2(n) = e(2) - mv(n);
    end
    medianX0 = (medianX - c0)/c1;
    errorbar(medianX0,mv,ci1,ci2,'ok','markersize',15)



    %%
    nhbins = 100;
    spacing = (max(X0fit) - min(X0fit))/nhbins;
    nx = linspace(min(X0fit)-spacing,max(X0fit),101);
    for n = 2:length(nx)
        if nx(n) == ceil(max(sX00))
            dbin0(n) = sum((sX00 >= nx(n-1)) & (sX00 <= nx(n)));
        else
            dbin0(n) = sum((sX00 >= nx(n-1)) & (sX00 < nx(n)));
        end
    end
    maxdbin0 = max(dbin0);
    for n = 2:length(nx)
        if nx(n) == ceil(max(sX01))
            dbin1(n) = sum((sX01 >= nx(n-1)) & (sX01 <= nx(n)));
        else
            dbin1(n) = sum((sX01 >= nx(n-1)) & (sX01 < nx(n)));
        end
    end
    maxdbin1 = max(dbin1);
    maxdbin = max([maxdbin0,maxdbin1]);

    for n = 1:length(nx)
        if dbin0(n) ~= 0
            if n > 1
                patch([nx(n-1),nx(n-1),nx(n),nx(n)],[0,-0.1*(dbin0(n)/maxdbin),-0.1*(dbin0(n)/maxdbin),0],...
                    [0.95,0,0],'linestyle','none');
            else
                nx0 = nx(1) - (nx(end) - nx(1))/(length(nx) - 1);
                patch([nx0,nx0,nx(n),nx(n)],[0,-0.1*(dbin0(n)/maxdbin),-0.1*(dbin0(n)/maxdbin),0],...
                    [0.95,0,0],'linestyle','none');
            end
        end
    end
    for n = 1:length(nx)
        if dbin1(n) ~= 0
            if n > 1
                patch([nx(n-1),nx(n-1),nx(n),nx(n)],[1,1+0.1*(dbin1(n)/maxdbin),1+0.1*(dbin1(n)/maxdbin),1],...
                    [0.2,0.2,1],'linestyle','none');
            else
                nx0 = nx(1) - (nx(end) - nx(1))/(length(nx) - 1);
                patch([nx0,nx0,nx(n),nx(n)],[0,-0.1*(dbin0(n)/maxdbin),-0.1*(dbin0(n)/maxdbin),0],...
                    [0.95,0,0],'linestyle','none');
            end
        end
    end

    plot(X0fit,Yfit,'-k','linewidth',2);

    for n = 1:length(addYfit)
        plot(X0fit,addYfit{n},':k','linewidth',2);
    end
    
    lh0 = plot(sX00,ysO0,'.','markersize',7,'linewidth',1,'color',[0.95,0,0]);
    
    plot(X0fit,zeros(size(X0fit)),'-k','linewidth',2)
    plot(X0fit,ones(size(X0fit)),'-k','linewidth',2)
    lh1 = plot(sX01,ysO1,'.','markersize',7,'linewidth',1,'color',[0.2,0.2,1]);
    
    
    %% Show JJ fit
    if showJJ
        JJYfit = 0.95./(1+(61.2./X0fit).^(4*1.5));
        plot(X0fit,JJYfit,':k','linewidth',2)
    end
    %%
    

    ylim([-0.11,1.11])
    xlim([X0fit(1),X0fit(end)])

    % xlablestr = ['gTD (EQD_{',num2str(RefEQDx),'Gy})'];
    % xlabel(xlablestr)
    ylabel('TCP')
    if isnumeric(RefUnits)
        xlabel(['gTD (EQD',num2str(RefUnits),' Gy)'])
    else
        xlabel(RefUnits)
    end
    customfigset(h,25)
    lh = legend([lh1,lh0],'Local Control','Local Failure');

    set(lh,'fontsize',18,'position',[0.143,0.735,0.164,0.076])

    set(h,'Name',outstr)
    hold off

end







% function [HLp,outstr] = draw_TCP_curve_logitfit(logisticmodel,Outcome,showHL,showTCP,RefUnits)
% 
% 
% if ~iscell(logisticmodel)
%     mdl = {logisticmodel};
% else
%     mdl = logisticmodel;
% end
% 
% 
% xpts = mdl{1}.Fitted.LinearPredictor;
% ypts = mdl{1}.Fitted.Response;
% 
% showtextdata = true;
% 
% bins = 5;
% binoverlap = 1;
% 
% 
% xpts0 = xpts;
% C = table2cell(mdl{1}.Coefficients);
% c0 = C{1,1};
% c1 = C{2,1};
% 
% xpts = (xpts - c0)/c1;
% 
% 
% [sX,sortind] = sort(xpts);
% [sX0,~] = sort(xpts0);
% sO = Outcome(sortind);
% sY = ypts(sortind);
% 
% % OLthresh = 2*std(sX);
% OLthresh = inf;
% 
% sX_cut_ind = find((sX < mean(sX) + OLthresh) & (sX > mean(sX) - OLthresh));    
% sX_cut_ind0 = find((sX0 < mean(sX0) + OLthresh) & (sX0 > mean(sX0) - OLthresh));
% sO_cut = sO(sX_cut_ind);   
% sX_cut = sX(sX_cut_ind);
% sX0_cut = sX0(sX_cut_ind0);
% 
% rStart = round(linspace(1,length(Outcome)-binoverlap*length(Outcome)/bins,bins));
% rEnd = round(linspace(1+binoverlap*length(Outcome)/bins,length(Outcome),bins));
% 
% 
% for n = 1:bins
%     if binoverlap == 1
%         r1 = round(length(sO_cut)*((n-1)/bins))+1;
%         r2 = round(length(sO_cut)*((n)/bins));
%         r10 = round(length(sX0_cut)*((n-1)/bins))+1;
%         r20 = round(length(sX0_cut)*((n)/bins));
%     else
%         r1 = rStart(n);
%         r2 = rEnd(n);
%         r10 = rStart(n);
%         r20 = rEnd(n);
%     end
%     Obin{n} = sO_cut(r1:r2);
%     Xbin{n} = sX_cut(r1:r2);
%     X0bin{n} = sX0_cut(r10:r20);
% end
% joint(1) = min(Xbin{1});
% for n = 1:bins-1
%     joint(n+1) = 0.5*(max(Xbin{n})+min(Xbin{n+1}));
% end
% joint(bins+1) = max(Xbin{bins});
% 
% for n = 1:bins
%     meanO(n) = mean(Obin{n});
%     medianX(n) = median(Xbin{n});
%     predictedO(n) = 1/(1 + exp(-medianX(n)));
%     binsize(n) = numel(Obin{n});
% end
% 
% 
% sX0 = sX(~sO);
% sX1 = sX(sO);
% sY0 = sY(~sO);
% sY1 = sY(sO);
% 
% for n = 1:length(sX0)
%     ysO0(n) = sY0(n) - 0.005 - 0.03*rand(1);
% end
% 
% ysO0(ysO0 < 0.005) = 0.005;
% for n = 1:length(sX1)
%     ysO1(n) = sY1(n) + 0.005 + 0.03*rand(1);
% end
% ysO1(ysO1 > 0.995) = 0.995;
% 
% rangeX = max(sX_cut) - min(sX_cut);
% 
% Xfit = linspace(min(sX_cut)-0.15*rangeX,max(sX_cut)+0.15*rangeX,1000);
% Yfit = 1./(1 + exp(-(c1*Xfit + c0)));
% 
% addYfit = cell(length(mdl)-1);
% for n = 2:length(mdl)
%     a1C = table2cell(mdl{n}.Coefficients);
%     a1c0 = a1C{1,1};
%     a1c1 = a1C{2,1};
%     addYfit{n-1} = 1./(1 + exp(-(a1c1*Xfit + a1c0)));
% end
% 
% %% Hosmer-Lemeshaw Test
% 
%     
% Z = nan(1,bins);
% PIg = nan(1,bins);
% Orate = nan(1,bins);
% for n = 1:bins
%     Og = sum(Obin{n});
%     Ng = numel(Obin{n});
%     PIg(1,n) = mean(1./(1 + exp(-X0bin{n})));
%     Orate(1,n) = Og/Ng;
%     
%     Z(1,n) = ((Og - Ng*PIg(1,n))^2)/(Ng*PIg(1,n)*(1-PIg(1,n)));    
% end
% 
% ChiSqHL = sum(Z);
% HLp = chi2cdf(ChiSqHL,bins-2,'upper');
% 
% 
% 
% eq = [mdl{1}.Formula.LinearPredictor,' '];
% eqb = eq + 0;
% eqb(eqb == 42) = 58;
% eq = char(eqb);
% T{1} = num2str(mdl{1}.Coefficients{1,1});
% count = 2;
% for n = 1:length(mdl{1}.CoefficientNames)
%     ckstr = [' ',mdl{1}.CoefficientNames{1,n},' '];
%     if strfind(eq,ckstr);
%         T{count} = [num2str(mdl{1}.Coefficients{n,1}),'*',mdl{1}.CoefficientNames{1,n}];
%         count = count + 1;
%     end
% end
% TT{1} = T{1};
% for n = 2:length(T)
%     Tq = T{n} + 0;
%     Tq(Tq == 58) = 42;
%     if Tq(1) == 45
%         Tq(4:end+2) = Tq(2:end);
%         Tq(1:3) = [32,45,32];
%     else
%         Tq(4:end+3) = Tq(1:end);
%         Tq(1:3) = [32,43,32];
%     end
%     TT{n} = char(Tq);
% end
% outstr = ['Y = '];
% for n = 1:length(TT)
%     outstr = [outstr,TT{n}];
% end
% 
% 
% 
% if showHL
%     %% Draw Calibration Plot
%     
%     for n = 1:bins
%         [d, e] = binofit(round(meanO(n)*binsize(n)),binsize(n));
%         mv(n) = d;
%         ci1(n) = mv(n) - e(1);
%         ci2(n) = e(2) - mv(n);
%     end
% 
%     hHL = figure;
% %     plot(PIg,Orate,'ok','markersize',15);
%     errorbar(PIg,Orate,ci1,ci2,'ok','markersize',15);
%     
%     hold on
%     plot([0,1],[0,1],'--k','linewidth',2)
%     xlim([max([0,floor(10*(min([min(meanO),min(predictedO)])-0.025))/10]),...   
%         min([1,ceil(10*(max([max(meanO),max(predictedO)])+0.025))/10])])
%     ylim([max([0,floor(10*(min([min(meanO),min(predictedO)])-0.025))/10]),...
%         min([1,ceil(10*(max([max(meanO),max(predictedO)])+0.025))/10])])
%     
%     xlim([0.5,1]),ylim([0.5,1])
%     
%     ylabel('Observed Control Rate')
%     xlabel('Predicted Control Rate')
%     title(['Hosmer-Lemeshow Test p = ',num2str(HLp)])
%     customfigset(hHL,15)
%     set(hHL,'Name',outstr)
% end
%     
% if showTCP 
%     %% Draw TCP Figure
% 
%     h = figure;
%     drawnow
%     for n = 1:bins
%         ph = patch([joint(n),joint(n),joint(n+1),joint(n+1)],...
%             [0,1,1,0],[1,1,1],'linewidth',2);
%         set(ph,'EdgeColor',[0.7,0.7,0.7]);
%         t(n) = text((joint(n)+joint(n+1))/2,0.05,['Q_',num2str(n)]);
%         set(t(n),'VerticalAlignment','baseline','HorizontalAlignment','center','Color',[0.7,0.7,0.7])
%     end
%     hold on
%     for n = 1:bins
%         [d, e] = binofit(round(meanO(n)*binsize(n)),binsize(n));
%         mv(n) = d;
%         ci1(n) = mv(n) - e(1);
%         ci2(n) = e(2) - mv(n);
%     end
%     errorbar(medianX,mv,ci1,ci2,'ok','markersize',15)
% 
% 
% 
%     %%
%     nhbins = 100;
%     spacing = (max(Xfit) - min(Xfit))/nhbins;
%     nx = linspace(min(Xfit)-spacing,max(Xfit),101);
%     for n = 2:length(nx)
%         if nx(n) == ceil(max(sX0))
%             dbin0(n) = sum((sX0 >= nx(n-1)) & (sX0 <= nx(n)));
%         else
%             dbin0(n) = sum((sX0 >= nx(n-1)) & (sX0 < nx(n)));
%         end
%     end
%     maxdbin0 = max(dbin0);
%     for n = 2:length(nx)
%         if nx(n) == ceil(max(sX1))
%             dbin1(n) = sum((sX1 >= nx(n-1)) & (sX1 <= nx(n)));
%         else
%             dbin1(n) = sum((sX1 >= nx(n-1)) & (sX1 < nx(n)));
%         end
%     end
%     maxdbin1 = max(dbin1);
%     maxdbin = max([maxdbin0,maxdbin1]);
% 
%     for n = 1:length(nx)
%         if dbin0(n) ~= 0
%             if n > 1
%                 patch([nx(n-1),nx(n-1),nx(n),nx(n)],[0,-0.1*(dbin0(n)/maxdbin),-0.1*(dbin0(n)/maxdbin),0],...
%                     [0.95,0,0],'linestyle','none');
%             else
%                 nx0 = nx(1) - (nx(end) - nx(1))/(length(nx) - 1);
%                 patch([nx0,nx0,nx(n),nx(n)],[0,-0.1*(dbin0(n)/maxdbin),-0.1*(dbin0(n)/maxdbin),0],...
%                     [0.95,0,0],'linestyle','none');
%             end
%         end
%     end
%     for n = 1:length(nx)
%         if dbin1(n) ~= 0
%             if n > 1
%                 patch([nx(n-1),nx(n-1),nx(n),nx(n)],[1,1+0.1*(dbin1(n)/maxdbin),1+0.1*(dbin1(n)/maxdbin),1],...
%                     [0.2,0.2,1],'linestyle','none');
%             else
%                 nx0 = nx(1) - (nx(end) - nx(1))/(length(nx) - 1);
%                 patch([nx0,nx0,nx(n),nx(n)],[0,-0.1*(dbin0(n)/maxdbin),-0.1*(dbin0(n)/maxdbin),0],...
%                     [0.95,0,0],'linestyle','none');
%             end
%         end
%     end
% 
%     plot(Xfit,Yfit,'-k','linewidth',2);
% 
%     for n = 1:length(addYfit)
%         plot(Xfit,addYfit{n},':k','linewidth',2);
%     end
%     
%     lh0 = plot(sX0,ysO0,'.','markersize',7,'linewidth',1,'color',[0.95,0,0]);
%     
%     plot(Xfit,zeros(size(Xfit)),'-k','linewidth',2)
%     plot(Xfit,ones(size(Xfit)),'-k','linewidth',2)
%     lh1 = plot(sX1,ysO1,'.','markersize',7,'linewidth',1,'color',[0.2,0.2,1]);
%     
%     
%     %% Show JJ fit
%     JJYfit = 0.95./(1+(61.2./Xfit).^(4*1.5));
%     plot(Xfit,JJYfit,':k','linewidth',2)
%     %%
%     
% 
%     ylim([-0.11,1.11])
%     xlim([Xfit(1),Xfit(end)])
% 
%     % xlablestr = ['gTD (EQD_{',num2str(RefEQDx),'Gy})'];
%     % xlabel(xlablestr)
%     ylabel('TCP')
%     if isnumeric(RefUnits)
%         xlabel(['gTD (EQD',num2str(RefUnits),' Gy)'])
%     else
%         xlabel(RefUnits)
%     end
%     customfigset(h,25)
%     lh = legend([lh1,lh0],'Local Control','Local Failure');
% 
%     set(lh,'fontsize',18,'position',[0.143,0.735,0.164,0.076])
% 
%     set(h,'Name',outstr)
%     hold off
% 
% end
% 
% 






