gEUDv = -100:10:100;
gEUDstats = zeros(size(Predictor_Table,2),length(gEUDv));
gEUDmodel = cell(0);
for n = 1:length(gEUDv)
    settings.gEUDa = gEUDv(n);
    Predictor_Table = MAIN_multivariate_logistic_gTD_model(BrainMet_Data_Table,settings);
    for m = 1:size(Predictor_Table,2)
        if isnumeric(Predictor_Table{1,m})
            gEUDstats(m,n) = Predictor_Table{11,m};
        else
            gEUDmodel{n} = Predictor_Table{11,m};
        end
    end
    disp(gEUDv(n))
end


plot(gEUDv,gEUDstats(12,:),'r',gEUDv,gEUDstats(11,:),'g',[gEUDv(1),gEUDv(end)],[0,0],'k'),legend('Predictor Coefficent','p-value')

figure,plot(gEUDv,gEUDstats(4,:)),title('LogLikelihood')
figure,plot(gEUDv,gEUDstats(7,:)),title('KM p-value')
figure,plot(gEUDv,gEUDstats(8,:)),title('Deviance')

keyboard