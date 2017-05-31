function [table1, table12] = km_sf(x)

x=sortrows(x,1);
%table of patients observed for each survival time
%the TABULATE function sets up this matrix:
%table1=[time count percent(on total)]
table1=[0 size(x,1) 1; tabulate(x(:,1))];
%if all observed time are integers remove not observed time added by
%TABULATE function
table1(table1(:,3)==0,:)=[];

%Table of censored data
table12=tabulate(x(x(:,2)==1));
if ~isempty(table12)
    % remove not observed time added by TABULATE function
    table12(table12(:,3)==0,:)=[];
    % setup the vector of the censored data
    [cens,loc]=ismember(table1(:,1),table12(:,1)); %find censored data
end    

%the percents stored in the the third column are unuseful;
%so, place in the third column how many subjects are still alive at the
%beginning of the i-th interval.
a1=[table1(1,2); -1.*table1(2:end,2)];
table1(:,3)=cumsum(a1); table1(2:end,3)=table1(1:end-1,3);
%number of deaths in the intervals (don't take in account the censored
%data)
if ~isempty(table12)
    table1(cens,2)=table1(cens,2)-table12(loc(cens),2);
end
%finally, delete the first row that is now useless
table1(1,:)=[];
