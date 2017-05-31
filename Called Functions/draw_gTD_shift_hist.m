function draw_gTD_shift_hist(gTD1,gTD2,q1,q2,Outcome)

cc = gTD1 - gTD2;
if ~isequal(cc,zeros(size(cc)))
    alp = 0.9;
    mincc = floor(min(cc));
    maxcc = ceil(max(cc));
    cc_1 = cc(Outcome);
    cc_0 = cc(~Outcome);
    f = figure;
    hold on
    inc = (maxcc - mincc)/100;
    hist(cc_1,mincc:inc:maxcc);
    hist(cc_0,mincc:inc:maxcc);

    ph = get(gca,'Children');
    if strcmp('cc_0',ph(1).DisplayName)
        set(ph(1),'FaceColor',[1,0.3,0.3],'FaceAlpha',alp)
        set(ph(2),'FaceColor',[0.3,0.3,1],'FaceAlpha',alp)
        lh = legend('Local Control','Local Failure');
    else
        set(ph(2),'FaceColor',[1,0.3,0.3],'FaceAlpha',alp)
        set(ph(1),'FaceColor',[0.3,0.3,1],'FaceAlpha',alp)
        lh = legend('Local Control','Local Failure');
    end
    xlim([mincc-0.025*(maxcc-mincc),maxcc+0.025*(maxcc-mincc)])
    xlabel(['gTD_{(q=',num2str(q1),')} - gTD_{(q=',num2str(q2),')}'])
    customfigset(f,18)
    set(lh,'FontSize',14)
end