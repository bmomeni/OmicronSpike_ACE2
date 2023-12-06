clear

listing = dir('*.fcs');
N = length(listing);
APCA_C = zeros(45,20000);

for n = 1:N
    
    filename = listing(n).name;
    [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(filename);
    
    FSCA = fcsdat(:,1);
    FSCH = fcsdat(:,2);
    SSCA = fcsdat(:,3);
    FITCA = fcsdat(:,4);
    APCA = fcsdat(:,7);
    
    %% First gate
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
    
    CountNeg(n) = size(APCA_P2(APCA_P2<5000),1);
    CountPos(n) = size(APCA_P2(APCA_P2>5000),1);
    
    MeanAPCA(n) = mean(APCA_P2(APCA_P2>5000));
    MedianAPCA(n) = median(APCA_P2(APCA_P2>5000));
    
    fl = length(filename);
    sample(n) = str2num(filename(fl-6:fl-4));
    NP2(sample(n)) = length(FSCA_P2);
    APCA_C(sample(n),1:NP2(sample(n))) = APCA_P2';
    
end

figure
bar(sample,MeanAPCA)
title('Mean')

figure
bar(sample,MedianAPCA)
title('Median')

disp('WT vs Omicron')
p1 = ranksum(APCA_C(12,1:NP2(12)),APCA_C(24,1:NP2(24)))

disp('WT vs Omicron-A484K')
p2 = ranksum(APCA_C(12,1:NP2(12)),APCA_C(36,1:NP2(36)))

disp('WT vs Omicron-L452R')
p3 = ranksum(APCA_C(12,1:NP2(12)),APCA_C(48,1:NP2(48)))

disp('Omicron vs Omicron-L452R')
p4 = ranksum(APCA_C(24,1:NP2(24)),APCA_C(48,1:NP2(48)))

CombinedData = [APCA_C(12,1:NP2(12)), APCA_C(24,1:NP2(24)), APCA_C(36,1:NP2(36)), APCA_C(48,1:NP2(48))];
GroupData = [ones(1,NP2(12)), 2*ones(1,NP2(24)), 3*ones(1,NP2(36)), 4*ones(1,NP2(48))];
[pkw,tbl,stats] = kruskalwallis(CombinedData,GroupData);

save('Experiment2_Data.mat')

