function set_options(data)

close all
f = figure('Position',[250 500 870 350]);

fpos = get(f,'Position');


tablevarnames_p = data.Properties.VariableNames;

count = 1;
for n = 1:length(tablevarnames_p)
   if isnumeric(data{1,n}) || islogical(data{1,n})
        tablevarnames{count} = tablevarnames_p{n};
        count = count + 1;
   end
end

% 
% cformat = {'numeric','numeric','numeric',...
%     'numeric','numeric','numeric',...
%     'numeric'};
% %     'bank','logical',{'Fixed' 'Adjustable'}};


t1 = uitable('Data', 75,...
             'ColumnName', 'Reference Volume');

t2 = uitable('Data', 2,...
             'ColumnName', 'Equivalent Dose');
         
t3 = uitable('Data', 0.39,...
             'ColumnName', 'Radiosensitivity Alpha');

t4 = uitable('Data', 6.8,...
             'ColumnName', 'Radiosensitivity Alpha/Beta');         
         
         

t5 = uitable('Data', 1,...
             'ColumnName', 'Test q-Value'); 
         
t6 = uitable('Data', -3,...
             'ColumnName', 'q Optimization Min');

t7 = uitable('Data', 2,...
             'ColumnName', 'q Optimization Max');
         
t8 = uitable('Data', 0.1,...
             'ColumnName', 'q Optimization Interval');
         
         
tarray = [t1,t2,t3,t4,t5,t6,t7,t8];
set(tarray,'RowName',[])
set(tarray,'ColumnWidth',{150})
set(tarray,'ColumnFormat',{'numeric'})
set(tarray,'ColumnEditable',true)


for n = 1:4
    set(tarray(n),'Position',[20, fpos(4) - 20 - 50*(n - 1) - tarray(n).Extent(4),...
        tarray(n).Extent(3), tarray(n).Extent(4)]);
end
for n = 5:8
    set(tarray(n),'Position',[180, fpos(4) - 20 - 50*(n - 5) - tarray(n).Extent(4),...
        tarray(n).Extent(3), tarray(n).Extent(4)]);
end

t10 = uitable('Data', {'Provide comma seperated list (e.g. "D90,V30,D5")'},...
             'ColumnName', 'Dx/Vx Predictors',...
             'RowName',[],...
             'ColumnWidth',{156 + tarray(1).Extent(3)},...
             'ColumnFormat',{'char'},...
             'ColumnEditable',true,...
             'Position',[20, fpos(4) - 220 - tarray(1).Extent(4),...
                160 + tarray(1).Extent(3), tarray(1).Extent(4)]);
            
t11 = uitable('Data', 9,...
             'ColumnName', 'gEUD a-Value',...
             'RowName',[],...
             'ColumnWidth',{150},...
             'ColumnFormat',{'numeric'},...
             'ColumnEditable',true);
            
set(t11,'Position',[20, fpos(4) - 20 - 50*(5) - tarray(n).Extent(4),...
        tarray(n).Extent(3), tarray(n).Extent(4)]);

rnames = [{'gTD (q=1)','gTD (q=test value)','gTD (optimized q)','Tumor Volume'},tablevarnames,{'gEUD'}];
cformat = cell(1,length(rnames));
ix = cellfun('isempty',cformat);
cformat(ix) = {'logical'};
t9 = uitable('Data', false(length(rnames),2),...
             'ColumnName',{'Predictors','Response Var.'},...
             'RowName',rnames,...
             'ColumnWidth',{150},...
             'ColumnFormat',cformat,...
             'ColumnEditable',true(1,length(rnames)),...
             'RowStriping','off');

set(t9,'Position',[350, fpos(4) - 20 - t9.Extent(4), t9.Extent(3),...
                t9.Extent(4)])

ts = {t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11};

tarray2 = [t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11];

set(tarray2,'CellSelectionCallback',{@deslectfn,ts})
set(tarray2,'CellEditCallback',{@deslectfn,ts})
            
uicontrol('Style','Pushbutton','Position',[180, fpos(4) - 20 - 50*5 - tarray(n).Extent(4),...
        tarray(n).Extent(3), tarray(n).Extent(4)],'String','Set Options','FontSize',11,...
        'FontWeight','Bold','FontAngle','Italic','Callback',{@getsettings,ts},'BackgroundColor',[1,0.5,0.5])
end
    
function settings = getsettings(~,~,ts)
    settings.ReferenceVolume = get(ts{1},'Data');
    settings.EQdose = get(ts{2},'Data');
    settings.Alpha = get(ts{3},'Data');
    settings.AlphaBeta = get(ts{4},'Data');
    settings.TestQVal = get(ts{5},'Data');
    settings.QOptMin = get(ts{6},'Data');
    settings.QOptMax = get(ts{7},'Data');
    settings.QOptInt = get(ts{8},'Data');
    settings.VarOnOff = get(ts{9},'Data');
    settings.VarNames = get(ts{9},'RowName');
    settings.gEUDa = get(ts{11},'Data');
    DVcell = get(ts{10},'Data');
    DVstr = DVcell{:}; 
    V = [];
    D = [];
    if strcmp(DVstr,'Provide comma seperated list (e.g. "D90,V30,D5")')
        settings.DVdata = [];
    else
        DVstr = strrep(DVstr, ' ', '');
        C = strsplit(DVstr,',');
        vcount = 1;
        dcount = 1;
        for n = 1:length(C)
            if strcmpi(C{n}(1),'V')
                V(vcount) = str2double(C{n}(2:end));
                vcount = vcount + 1;
            end
            if strcmpi(C{n}(1),'D')
                D(dcount) = str2double(C{n}(2:end));
                dcount = dcount + 1;
            end
        end
    end
    if ~isempty(D)
        D(D > 100) = [];
        D(D < 0) = [];
    end
    if ~isempty(V)
        V(V > 100) = [];
        V(V < 0) = [];
    end
    settings.DxVals = D;
    settings.VxVals = V;
    [filename, pathname, ~] = uiputfile('*.mat','Save Settings File');
    save(strcat(pathname,filename),'settings')
    close all
end
     
function deslectfn(hO,~,ts)
    hname = get(hO,'ColumnName');
    for n = 1:length(ts)
        tsname = get(ts{n},'ColumnName');
        if ~strcmp(hname,tsname)
            d = get(ts{n},'Data');
            set(ts{n}, 'Data', []);
            set(ts{n}, 'Data', d);
        end
    end
    if strcmp(get(hO,'ColumnName'),'Dx/Vx Predictors'),
        DVstr = get(hO,'Data');
        if isempty(DVstr{:})
            set(hO,'Data',{'Provide comma seperated list (e.g. "D90, V30, D5")'})
        end
    end
end
    