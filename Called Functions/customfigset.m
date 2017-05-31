function customfigset(h,fsize)

if nargin < 2
    fsize = 12;
end

% ss = get(0,'monitorpositions');
% ss2 = ss(1,:);
% ss2(3) = ss2(3) - ss2(1) + 1;
% ss2(4) = ss2(4) - ss2(2) + 1;
% set(h,'position',ss2)
a = findall(h);
b = isprop(a,'fontsize');
set(a(b),'fontsize',fsize)
% hh = gca;
% set(hh,'position',[0.065, 0.055, 0.8875, 0.9075])
% 
% 
% set(hh,'position',[0.08, 0.08, 0.87, 0.87])