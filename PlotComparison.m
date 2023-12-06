clear

load('Experiment1_Data.mat')
APCA_C1 = APCA_C;
NP21 = NP2;
th1 = 6000;

load('Experiment2_Data.mat')
APCA_C2 = APCA_C;
NP22 = NP2;
th2 = 5000;

load('Experiment3_Data_th12000.mat')
APCA_C3 = APCA_C;
NP23 = NP2;
th3 = 12000;

% WT, Omicron, Omicron-A484K, Omicron-L452R
% Each row is from an independent experiment (biological replicate)
APCA_mean = [mean(APCA_C1(14,1:NP21(14))), mean(APCA_C1(20,1:NP21(20))), mean(APCA_C1(26,1:NP21(26))), mean(APCA_C1(32,1:NP21(32))); 
    mean(APCA_C2(12,1:NP22(12))), mean(APCA_C2(24,1:NP22(24))), mean(APCA_C2(36,1:NP22(36))), mean(APCA_C2(48,1:NP22(48)));
    mean(APCA_C3(10,1:NP23(10))), mean(APCA_C3(22,1:NP23(22))), mean(APCA_C3(34,1:NP23(34))), mean(APCA_C3(46,1:NP23(46)))];

APCAth_mean = [mean(APCA_C1(14,APCA_C1(14,1:NP21(14))>th1)), mean(APCA_C1(20,APCA_C1(20,1:NP21(20))>th1)), mean(APCA_C1(26,APCA_C1(26,1:NP21(26))>th1)), mean(APCA_C1(32,APCA_C1(32,1:NP21(32))>th1)); 
    mean(APCA_C2(12,APCA_C2(12,1:NP22(12))>th2)), mean(APCA_C2(24,APCA_C2(24,1:NP22(24))>th2)), mean(APCA_C2(36,APCA_C2(36,1:NP22(36))>th2)), mean(APCA_C2(48,APCA_C2(48,1:NP22(48))>th2));
    mean(APCA_C3(10,APCA_C3(10,1:NP23(10))>th3)), mean(APCA_C3(22,APCA_C3(22,1:NP23(22))>th3)), mean(APCA_C3(34,APCA_C3(34,1:NP23(34))>th3)), mean(APCA_C3(46,APCA_C3(46,1:NP23(46))>th3))];

APCA_median = [median(APCA_C1(14,1:NP21(14))), median(APCA_C1(20,1:NP21(20))), median(APCA_C1(26,1:NP21(26))), median(APCA_C1(32,1:NP21(32))); 
    median(APCA_C2(12,1:NP22(12))), median(APCA_C2(24,1:NP22(24))), median(APCA_C2(36,1:NP22(36))), median(APCA_C2(48,1:NP22(48)));
    median(APCA_C3(10,1:NP23(10))), median(APCA_C3(22,1:NP23(22))), median(APCA_C3(34,1:NP23(34))), median(APCA_C3(46,1:NP23(46)))];

APCAth_median = [median(APCA_C1(14,APCA_C1(14,1:NP21(14))>th1)), median(APCA_C1(20,APCA_C1(20,1:NP21(20))>th1)), median(APCA_C1(26,APCA_C1(26,1:NP21(26))>th1)), median(APCA_C1(32,APCA_C1(32,1:NP21(32))>th1)); 
    median(APCA_C2(12,APCA_C2(12,1:NP22(12))>th2)), median(APCA_C2(24,APCA_C2(24,1:NP22(24))>th2)), median(APCA_C2(36,APCA_C2(36,1:NP22(36))>th2)), median(APCA_C2(48,APCA_C2(48,1:NP22(48))>th2));
    median(APCA_C3(10,APCA_C3(10,1:NP23(10))>th3)), median(APCA_C3(22,APCA_C3(22,1:NP23(22))>th3)), median(APCA_C3(34,APCA_C3(34,1:NP23(34))>th3)), median(APCA_C3(46,APCA_C3(46,1:NP23(46))>th3))];
 
% figure
% for ii = 1:4
% bar(ii,mean(APCA_mean(:,ii)))
% hold on
% errorbar(ii,mean(APCA_mean(:,ii)),std(APCA_mean(:,ii)),'')
% end
% set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})

% figure
% bar(1:4,mean(APCA_mean))
% hold on
% plot(1:4,APCA_mean,'.','MarkerSize',20)
% errorbar(1:4,mean(APCA_mean),std(APCA_mean),'.','MarkerSize',1)
% set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
% title('mean')
% 
% figure
% bar(1:4,mean(APCA_median))
% hold on
% plot(1:4,APCA_median,'.','MarkerSize',20)
% errorbar(1:4,mean(APCA_median),std(APCA_median),'.','MarkerSize',1)
% set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
% title('median')
% 
figure
bar(1:4,mean(APCAth_mean))
hold on
plot(1:4,APCAth_mean,'.','MarkerSize',20)
errorbar(1:4,mean(APCAth_mean),std(APCAth_mean),'.','MarkerSize',1)
set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
ylabel('Fluorescence')
title('mean w/ th')

figure
bar(1:4,[595 820 1050 825])
set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
ylabel('Binding energy (kcal/mol.)')
ylim([0 1200])

figure
bar(1:4,mean(APCAth_median))
hold on
plot(1:4,APCAth_median,'.','MarkerSize',20)
errorbar(1:4,mean(APCAth_median),std(APCAth_median),'.','MarkerSize',1)
set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
title('median w/ th')

figure
bar(1:4,mean(APCAth_mean./APCAth_mean(:,1)))
% hold on
% plot(1:4,APCAth_mean./APCAth_mean(:,1),'.','MarkerSize',20)
% errorbar(1:4,mean(APCAth_mean./APCAth_mean(:,1)),std(APCAth_mean./APCAth_mean(:,1)),'.','MarkerSize',1)
% set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
title('mean w/ th')

figure
bar(1:4,(APCAth_mean./APCAth_mean(:,1)))
% hold on
% plot(1:4,APCAth_mean./APCAth_mean(:,1),'.','MarkerSize',20)
% errorbar(1:4,mean(APCAth_mean./APCAth_mean(:,1)),std(APCAth_mean./APCAth_mean(:,1)),'.','MarkerSize',1)
% set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
title('mean w/ th')

figure
bar(1:3,(APCAth_mean./APCAth_mean(:,1))')
% hold on
% plot(1:4,APCAth_mean./APCAth_mean(:,1),'.','MarkerSize',20)
% errorbar(1:4,mean(APCAth_mean./APCAth_mean(:,1)),std(APCAth_mean./APCAth_mean(:,1)),'.','MarkerSize',1)
% set(gca,'XTick',1:4,'XTickLabel',{'Wuhan', 'Omicron', 'Omicron-A484K', 'Omicron-L452R'})
title('mean w/ th')

