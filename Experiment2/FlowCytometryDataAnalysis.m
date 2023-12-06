clear

filename = '071522BF_1_001.fcs';
[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(filename);

FSCA = fcsdat(:,1);
FSCH = fcsdat(:,2);
SSCA = fcsdat(:,3);
FITCA = fcsdat(:,4);
APCA = fcsdat(:,7);
figure
plot(FSCA,SSCA,'.');
xlabel('FSC-A')
ylabel('SSC-A')
xlim([0 max(fcsdat(:,1))])
ylim([0 max(fcsdat(:,3))])

%% First gate
% x = 250000 * rand(1,5000);
% y = 250000 * rand(1,5000);
x = FSCA;
y = SSCA;

p1 = 1000*[60,28];
p2 = 1000*[120,10];
p3 = 1000*[230,83];
p4 = 1000*[240,181];
p5 = 1000*[171,188];
p6 = 1000*[73,110];


l1 = -y + p1(2) + (p2(2) - p1(2))/(p2(1)-p1(1))*(x - p1(1));
l2 = -y + p2(2) + (p3(2) - p2(2))/(p3(1)-p2(1))*(x - p2(1));
l3 = -y + p3(2) + (p4(2) - p3(2))/(p4(1)-p3(1))*(x - p3(1));
l4 = -y + p4(2) + (p5(2) - p4(2))/(p5(1)-p4(1))*(x - p4(1));
l5 = -y + p5(2) + (p6(2) - p5(2))/(p6(1)-p5(1))*(x - p5(1));
l6 = -y + p6(2) + (p1(2) - p6(2))/(p1(1)-p6(1))*(x - p6(1));

FSCA_P1 = x((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0));
FSCH_P1 = FSCH((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0));
SSCA_P1 = y((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0));
FITCA_P1 = FITCA((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0));
APCA_P1 = APCA((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0));
% hold on
% plot(x((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0)),y((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6>0)),'r.')
% plot(x((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6<0)),y((l1<0)&(l2<0)&(l3<0)&(l4>0)&(l5>0)&(l6<0)),'b.')
hold on
plot(FSCA_P1,SSCA_P1,'r.');

%% Second gate
x = FSCA_P1;
y = FSCH_P1;

q1 = 1000*[61,14];
q2 = 1000*[258,72];
l1 = -y + q1(2) + (q2(2) - q1(2))/(q2(1)-q1(1))*(x - q1(1));

FSCA_P2 = FSCA_P1((l1<0));
FSCH_P2 = FSCH_P1((l1<0));
APCA_P2 = APCA_P1((l1<0));
FITCA_P2 = FITCA_P1((l1<0));

figure
plot(FSCA_P1,FSCH_P1,'.');
xlabel('FSC-A')
ylabel('FSC-H')
xlim([0 250000])
ylim([0 250000])
hold on
plot(FSCA_P2,FSCH_P2,'b.');

figure
plot(FITCA_P2,APCA_P2,'.');
xlabel('FITC-A')
ylabel('APC-A')
xlim([0 1000])
ylim([0 1000])
hold on
plot(FSCA_P2,FSCH_P2,'r.');

CountNeg = size(APCA_P2(APCA_P2<5000),1);
CountPos = size(APCA_P2(APCA_P2>5000),1);

MeanAPCA = mean(APCA_P2(APCA_P2>5000));
MedianAPCA = median(APCA_P2(APCA_P2>5000));

disp([CountNeg,CountPos])
disp([MeanAPCA, MedianAPCA])

