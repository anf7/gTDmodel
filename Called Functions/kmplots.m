function [p_value,MedFailureHighRisk,MedFailureLowRisk] = kmplots(Predictor,FailureEvent,Followup,varargin)
% kmplots creates Kaplan-Meier plots for the given data and reports statistical signifigance.
% 
% Code written by Andrew Fontanella for Joe Deasy's lab at MSKCC
% 
% Function form:
%     kmplots(Predictor,FailureEvent,Followup,'Parameter1','Value1','Parameter2','Value2',...)
% 
% Required inputs (these must be equal in size):
%     Predictor - a one dimensional array of numerical predictor values for each patient
%     FailureEvent - a logical array indicating whether a failure event occured for that patient
%     Followup - a numerical array indicating the time to observed event (FailureEvent = 1) or right censoring (FailureEvent = 0)
% 
% Optional Parameter/Value Pairs:
%     'PredictorName'         A string containing the name of the predictor variable
%     'ConfidenceInterval'    A numerical scalar value indicating the percent confidence interval to display 
%                               [Allowed values are between 0 and 100 (not inclusive)]
%                               [Default = 95]
%     'TimeUnits'             Units of time for the Followup variable 
%                               [Allowed values are |'days'|, |'weeks'|, |'months'|, or |'years'|] 
%                               [Default = 'years']
%     'AddTitle'              A logical scalar indicating whether to include a title in the plot
%                               [Allowed values are |true| or |false|] 
%                               [Default = false]
%     'CutoffTime'            A numerical scalar value representing the maximum followup time cutoff (data with Followup 
%                             values greater than this value will be excluded) 
%                               [Default = inf] (i.e. no data cutoff)
%     'WeibulFit'             A logical scalar indicating whether a Weibul Fit is superimposed over the survival curves
%                               [Allowed values are |true| or |false|]
%                               [Default = false]
%     'FractionalGroupSize'   A numerical scalar value indicating the fraction predictor-ranked cases to include in 
%                             each risk group
%                               [Allowed values are between 0 and 1 (not inclusive)]
%                               [Default = 1/3]
%     'RiskCorrelation'       A string indicating the proposed correlation between the predictor and the risk of failure
%                               [Allowed values are |'direct'| or |'inverse'|]
%                               [Default = 'inverse']

p = inputParser;

defaultPredictorName = 'Predictor Variable';
defaultCI = 95;
defaultTimeUnits = 'years';
defaultTitleFlag = false;
defaultCutoffTime = inf;
defaultWeibul = false;
defaultGroupFraction = 1/3;
defaultCorrelation = 'inverse';
defalutShowFigures = true;

expectedTimeUnits = {'days','weeks','months','years'};
expectedCorrelation = {'direct','inverse'};

addRequired(p,'Predictor',@isnumeric);
addRequired(p,'FailureEvent',@islogical);
addRequired(p,'Followup',@isnumeric);

addParameter(p,'PredictorName',defaultPredictorName,@ischar);
addParameter(p,'ConfidenceInterval',defaultCI,@isnumeric);
addParameter(p,'TimeUnits',defaultTimeUnits,...
    @(x) any(validatestring(x,expectedTimeUnits)));
addParameter(p,'AddTitle',defaultTitleFlag,@islogical)
addParameter(p,'CutoffTime',defaultCutoffTime,@isnumeric);
addParameter(p,'AddWeibulFit',defaultWeibul,@islogical)   
addParameter(p,'FractionalGroupSize',defaultGroupFraction,@isnumeric);
addParameter(p,'RiskCorrelation',defaultCorrelation,...
    @(x) any(validatestring(x,expectedCorrelation))); 
addParameter(p,'ShowFigures',defalutShowFigures,@islogical)

parse(p,Predictor,FailureEvent,Followup,varargin{:});

PredictorNameString =  p.Results.PredictorName;
KMConfInterval = p.Results.ConfidenceInterval/100;
if KMConfInterval >= 1 || KMConfInterval <= 0;
    error('ConfidenceInterval parameter value must be >0 and <100')
end
TimeUnits = p.Results.TimeUnits;
TitleFlag = p.Results.AddTitle;
ShowFigures = p.Results.ShowFigures;
CutoffTime = p.Results.CutoffTime;
if CutoffTime <= 0
    error('CutoffTime parameter value must be positive')
end
WeibulFit = p.Results.AddWeibulFit;
FractionalGroupSize = p.Results.FractionalGroupSize;
if FractionalGroupSize >= 1 || FractionalGroupSize <= 0;
    error('FractionalGroupSize parameter value must be >0 and <1')
end
PedictorRiskCorrelation = p.Results.RiskCorrelation;

[~,PredictorRankIndex] = sort(Predictor);
PredictorRanked_Censored = ~FailureEvent(PredictorRankIndex);
PredictorRanked_Followup = Followup(PredictorRankIndex);

PredictorRanked_Censored(PredictorRanked_Followup > CutoffTime) = [];
PredictorRanked_Followup(PredictorRanked_Followup > CutoffTime) = [];

lowI = [1,round(FractionalGroupSize*length(PredictorRanked_Censored))];
highI = [round((1-FractionalGroupSize)*length(PredictorRanked_Censored)),length(PredictorRanked_Censored)];

if strcmp(PedictorRiskCorrelation,'inverse')
    Followup_highrisk = PredictorRanked_Followup(lowI(1):lowI(2));
    Censored_highrisk = PredictorRanked_Censored(lowI(1):lowI(2));
    Followup_lowrisk = PredictorRanked_Followup(highI(1):highI(2));
    Censored_lowrisk = PredictorRanked_Censored(highI(1):highI(2));
else
    Followup_highrisk = PredictorRanked_Followup(highI(1):highI(2));
    Censored_highrisk = PredictorRanked_Censored(highI(1):highI(2));
    Followup_lowrisk = PredictorRanked_Followup(lowI(1):lowI(2));
    Censored_lowrisk = PredictorRanked_Censored(lowI(1):lowI(2));
end

[~,TimeSortIndex_high] = sort(Followup_highrisk);
TimeRanked_Followup_highrisk = Followup_highrisk(TimeSortIndex_high);
TimeRanked_Censored_highrisk = Censored_highrisk(TimeSortIndex_high);

[~,TimeSortIndex_low] = sort(Followup_lowrisk);
TimeRanked_Followup_lowrisk = Followup_lowrisk(TimeSortIndex_low);
TimeRanked_Censored_lowrisk = Censored_lowrisk(TimeSortIndex_low);


MedFailureHighRisk = median(TimeRanked_Followup_highrisk(~TimeRanked_Censored_highrisk));
MedFailureLowRisk = median(TimeRanked_Followup_lowrisk(~TimeRanked_Censored_lowrisk));

% [~,~,~,stats] = coxphfit(TimeRanked_Group,TimeRanked_Followup,'censoring',TimeRanked_Censored);


p_value = logrank([TimeRanked_Followup_lowrisk, TimeRanked_Censored_lowrisk],...
    [TimeRanked_Followup_highrisk, TimeRanked_Censored_highrisk],1-KMConfInterval);


if ShowFigures

    fHandle = figure;
    hold on

    [fl,xl,fl_lower,fl_upper] = ecdf([TimeRanked_Followup_lowrisk;...
        max([max(TimeRanked_Followup_highrisk),max(TimeRanked_Followup_lowrisk)])+1],...
        'censoring',[TimeRanked_Censored_lowrisk;0],'function','survivor','alpha',1-KMConfInterval);

    y_upper = fl_upper;
    y_upper(1) = 1;
    y_upper(end) = y_upper(end-1);
    y_lower = fl_lower;
    y_lower(1) = 1;
    y_lower(end) = y_lower(end-1);
    counter = 1;
    x1 = nan(1,2*length(fl_upper)-1);
    y1 = nan(1,2*length(fl_upper)-1);
    for n = 1:length(fl_upper);
        if n == 1
            x1(counter) = xl(n);
            y1(counter) = y_upper(n);
            counter = counter + 1;
        else
            x1(counter) = xl(n);
            y1(counter) = y_upper(n-1);
            counter = counter + 1;
            x1(counter) = xl(n);
            y1(counter) = y_upper(n);
            counter = counter + 1;
        end
    end
    counter = 1;
    x2 = nan(1,2*length(fl_upper)-1);
    y2 = nan(1,2*length(fl_upper)-1);
    for n = 1:length(fl_lower);
        if n == 1
            x2(counter) = xl(n);
            y2(counter) = y_lower(n);
            counter = counter + 1;
        else
            x2(counter) = xl(n);
            y2(counter) = y_lower(n-1);
            counter = counter + 1;
            x2(counter) = xl(n);
            y2(counter) = y_lower(n);
            counter = counter + 1;
        end
    end
    xb = [x1,fliplr(x2)];
    yb = [y1,fliplr(y2)];
    % patch(xb,yb,[0,0,1],'FaceAlpha',0.15,'LineStyle','none')
    patch(xb,yb,[0,0,1],'FaceAlpha',0.15,'LineWidth',0.5,'EdgeColor',[0.5,0.5,1])

    counter = 0;
    Clx = nan(1,length(TimeRanked_Censored_lowrisk));
    Cly = nan(1,length(TimeRanked_Censored_lowrisk));
    for m = 1:length(TimeRanked_Censored_lowrisk)
        if TimeRanked_Censored_lowrisk(m)
            counter = counter + 1;
            for n = 1:length(xl)
                if n == 1 && TimeRanked_Followup_lowrisk(m) <= xl(n)
                    Cly(1,counter) = 1;
                    Clx(1,counter) = TimeRanked_Followup_lowrisk(m);
                elseif n > 1 && TimeRanked_Followup_lowrisk(m) > xl(n-1) && TimeRanked_Followup_lowrisk(m) < xl(n)
                    Cly(1,counter) = fl(n-1);
                    Clx(1,counter) = TimeRanked_Followup_lowrisk(m);
                end
            end
        end
    end
    Clx(isnan(Clx)) = [];
    Cly(isnan(Cly)) = [];

    [fh,xh,fh_lower,fh_upper] = ecdf([TimeRanked_Followup_highrisk;max([max(TimeRanked_Followup_highrisk),max(TimeRanked_Followup_lowrisk)])+1],'censoring',[TimeRanked_Censored_highrisk;0],'function','survivor','alpha',1-KMConfInterval);
    y_upper = fh_upper;
    y_upper(1) = 1;
    y_upper(end) = y_upper(end-1);
    y_lower = fh_lower;
    y_lower(1) = 1;
    y_lower(end) = y_lower(end-1);
    counter = 1;
    x1 = nan(1,2*length(fh_upper)-1);
    y1 = nan(1,2*length(fh_upper)-1);
    for n = 1:length(fh_upper);
        if n == 1
            x1(counter) = xh(n);
            y1(counter) = y_upper(n);
            counter = counter + 1;
        else
            x1(counter) = xh(n);
            y1(counter) = y_upper(n-1);
            counter = counter + 1;
            x1(counter) = xh(n);
            y1(counter) = y_upper(n);
            counter = counter + 1;
        end
    end
    counter = 1;
    x2 = nan(1,2*length(fh_upper)-1);
    y2 = nan(1,2*length(fh_upper)-1);
    for n = 1:length(fh_lower);
        if n == 1
            x2(counter) = xh(n);
            y2(counter) = y_lower(n);
            counter = counter + 1;
        else
            x2(counter) = xh(n);
            y2(counter) = y_lower(n-1);
            counter = counter + 1;
            x2(counter) = xh(n);
            y2(counter) = y_lower(n);
            counter = counter + 1;
        end
    end
    xb = [x1,fliplr(x2)];
    yb = [y1,fliplr(y2)];
    % patch(xb,yb,[1,0,0],'FaceAlpha',0.15,'LineStyle','none')
    patch(xb,yb,[1,0,0],'FaceAlpha',0.15,'LineWidth',0.5,'EdgeColor',[1,0.5,0.5])


    counter = 0;
    Chx = nan(1,length(TimeRanked_Censored_lowrisk));
    Chy = nan(1,length(TimeRanked_Censored_lowrisk));
    for m = 1:length(TimeRanked_Censored_highrisk)
        if TimeRanked_Censored_highrisk(m)
            counter = counter + 1;
            for n = 1:length(xh)
                if n == 1 && TimeRanked_Followup_highrisk(m) <= xh(n)
                    Chy(1,counter) = 1;
                    Chx(1,counter) = TimeRanked_Followup_highrisk(m);
                elseif n > 1 && TimeRanked_Followup_highrisk(m) > xh(n-1) && TimeRanked_Followup_highrisk(m) < xh(n)
                    Chy(1,counter) = fh(n-1);
                    Chx(1,counter) = TimeRanked_Followup_highrisk(m);
                end
            end
        end
    end
    Chx(isnan(Chx)) = [];
    Chy(isnan(Chy)) = [];

    ax1 = stairs(xl,fl,'-b','linewidth',2);
    plot(Clx,Cly,'+b','MarkerSize',7)
    ax2 = stairs(xh,fh,'-r','linewidth',2);
    plot(Chx,Chy,'+r','MarkerSize',7)


    if WeibulFit
        x = 0:max([max(TimeRanked_Followup_highrisk),max(TimeRanked_Followup_lowrisk)])+1;

        parmhat1 = wblfit(TimeRanked_Followup_lowrisk,[],TimeRanked_Censored_lowrisk);
        wblsurv = 1-cdf('weibull',x,parmhat1(1),parmhat1(2));
        ax3 = plot(x,wblsurv,'linewidth',2,'color',[0.65,0.65,1]);

        parmhat2 = wblfit(TimeRanked_Followup_highrisk,[],TimeRanked_Censored_highrisk);
        wblsurv = 1-cdf('weibull',x,parmhat2(1),parmhat2(2));
        ax4 = plot(x,wblsurv,'linewidth',2,'color',[1,0.65,0.65]);

        lh = legend([ax1,ax2,ax3,ax4],{'Low Risk Group','High Risk Group','Low Risk Weibull Fit','High Risk Weibull Fit'},'location','southwest');
    else
        lh = legend([ax1,ax2],{'Low Risk Group','High Risk Group'},'location','southwest');
    end

    ylim([0,1.0])
    xlim([0,max([max(TimeRanked_Followup_highrisk),max(TimeRanked_Followup_lowrisk)])])

    titlestr = ['p = ',num2str(p_value,3)];
    set(fHandle,'Name',PredictorNameString)
    xlabel(['Time (',TimeUnits,')'])
    ylabel('Survival')
    customfigset(fHandle,18)
    set(lh,'fontsize',13)
    grid on
    hold off
    if TitleFlag
        title(titlestr)
    else
        msgbox(['p = ',num2str(p_value,3)],'Statistical Significance')
    end


    if WeibulFit

        fHandle = figure;
        hold on
        hazard1 = pdf('wbl',x,parmhat1(1),parmhat1(2))./(1-cdf('wbl',x,parmhat1(1),parmhat1(2)));
        hazard2 = pdf('wbl',x,parmhat2(1),parmhat2(2))./(1-cdf('wbl',x,parmhat2(1),parmhat2(2)));
        ax1 = plot(x,hazard1,'-b','linewidth',2);
        ax2 = plot(x,hazard2,'-r','linewidth',2);
        lh = legend([ax1,ax2],{'Low Risk Group','High Risk Group'},'location','northwest');
        xlim([0,max([max(TimeRanked_Followup_highrisk),max(TimeRanked_Followup_lowrisk)])])
        titlestr = ['Hazard Functions for Weibull Fitted ',PredictorNameString,' Groupings'];
        if TitleFlag
            title(titlestr)
        end
        xlabel(['Time (',TimeUnits,')'])
        ylabel('Hazzard')

        customfigset(fHandle,18)
        set(lh,'fontsize',13)
        hp = get(fHandle,'Position');
        nhp = [hp(1),hp(2),670,540];
        set(fHandle,'Position',nhp)
        hold off
    end
end