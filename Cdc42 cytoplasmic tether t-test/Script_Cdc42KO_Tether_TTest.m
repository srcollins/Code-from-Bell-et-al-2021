%% plot cytoplasmic bridge data

%Ctrl for 2/6/21=1/32
% KO for 2/6/21=17/37

 %Ctrl for 2/9/21=2/56
 % KO for 2/9/21=14/26
 
 %Ctrl for 2/11/21 = 4/89;
 %KO for 2/11/21 = 20/60;
 
 %Ctrl for 2/12/21 = 5/154
 %KO for 2/21/21= 18/71;

Ctrl=[1/32,2/56, 4/89, 5/154];
KO=[17/37, 14/26, 20/60, 18/71];

meanCtrl=mean(Ctrl);
meanKO=mean(KO);
semCtrl=std(Ctrl)/sqrt(numel(Ctrl));
semKO=std(KO)/sqrt(numel(KO));


xAx=categorical({'Ctrl', 'Cdc42-KO'});
xAx=reordercats(xAx,{'Ctrl', 'Cdc42-KO'});
yVals=[meanCtrl, meanKO];
erVals=[semCtrl, semKO];
bar(xAx,yVals,'FaceColor','none'); hold on;
er=errorbar(xAx,yVals,erVals, 'LineWidth',1.5);
er.Color='k';
er.LineStyle='none';
ylabel('Relavtive Frequency');
 
%scatter ctrl
xC=repmat(xAx(1),1,4);
scatter(xC,Ctrl,50,'k','filled','jitter','on','jitterAmount',0.25);
xKo=repmat(xAx(2),1,4);
scatter(xKo,KO,50,'k','filled','jitter','on','jitterAmount',0.25);
%% t-test output

[h,p,ci,stats]=ttest2(Ctrl,KO,'Tail','both','Vartype','unequal')