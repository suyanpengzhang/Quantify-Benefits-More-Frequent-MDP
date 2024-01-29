data = [
39826.6	59202.8	9696.67	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000;
50729.8	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000;
54973.7	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000;
46783.9	78097.9	79937.2	21920.1	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000;
];
%% basecase
x = 1:2:28;
y11 = data(1,1:2);
y12 = data(1,2:3);
y13 = data(1,3:end);
y21 = data(2,1);
y23 = data(2,1:end);
y32 = data(3,1);
y33 = data(3,1:end);

figure
p = plot(x(1:2),y11,x(2:3),y12,'--',x(3:end),y13,':', ...
    x(1),y21,x(1:end),y23,':', ...
    x(1),y32,'--',x(1:end),y33,':');
set(gca,"FontSize",18)
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(1).Marker = '>';
p(1).MarkerIndices = [1,2,3];
p(1).MarkerSize = 10;
p(2).Marker = '+';
p(2).MarkerIndices = [2,3];
p(2).MarkerSize = 10;
p(3).Marker = '*';
p(3).MarkerIndices = 2:1:length(y13);
p(3).MarkerSize = 10;

p(4).LineWidth = 3;
p(5).LineWidth = 3;
p(4).Marker = '>';
p(4).MarkerIndices = [1];
p(4).MarkerSize = 10;
p(5).Marker = '*';
p(5).MarkerIndices = 2:1:length(y23);
p(5).MarkerSize = 10;

p(6).LineWidth = 3;
p(7).LineWidth = 3;
p(6).Marker = '>';
p(6).MarkerIndices = [1];
p(6).MarkerSize = 10;
p(7).Marker = '*';
p(7).MarkerIndices = 2:1:length(y33);
p(7).MarkerSize = 10;

p(1).MarkerEdgeColor = [0,0,1];
p(2).MarkerEdgeColor = [0,0,1];
p(3).MarkerEdgeColor = [0,0,1];
p(4).MarkerEdgeColor = [1,0,0];
p(5).MarkerEdgeColor = [1,0,0];
p(6).MarkerEdgeColor = [0,0,0];
p(7).MarkerEdgeColor = [0,0,0];
p(1).Color = [0,0,1];
p(2).Color = [0,0,1];
p(3).Color = [0,0,1];
p(4).Color = [1,0,0];
p(5).Color = [1,0,0];
p(6).Color = [0,0,0];
p(7).Color = [0,0,0];

legend({'ACLF2, L: wait, M: wait ','ACLF2, L: accept, M: wait ','ACLF2, L: accept, M: accept ', ...
    'ACLF=3OF, L: wait, M: wait ','ACLF=3OF, L: accept, M: accept ', ...
    'ACLF>3OF, L: wait, M: wait ','ACLF>3OF, L: accept, M: accept '},FontSize=18)
xlabel('Day',FontSize=16) 
ylabel('D: Difference in Lifetime Total Expected Value ($)',FontSize=16) 

title('rr = 0.8, o = 70%',FontSize=18 )
%% fig2
x = 1:2:28;
y11 = data(4,1:2);
y12 = data(4,2:4);
y13 = data(4,4:end);
figure
p = plot(x(1:2),y11,x(2:4),y12,'--',x(4:end),y13,':');
set(gca,"FontSize",18)
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(1).Marker = '>';
p(1).MarkerIndices = [1,2];
p(1).MarkerSize = 10;
p(2).Marker = '+';
p(2).MarkerIndices = [2,3];
p(2).MarkerSize = 10;
p(3).Marker = '*';
p(3).MarkerIndices = 2:1:length(y13);
p(3).MarkerSize = 10;


p(1).MarkerEdgeColor = [0,0,0];
p(2).MarkerEdgeColor = [0,0,0];
p(3).MarkerEdgeColor = [0,0,0];

p(1).Color = [0,0,0];
p(2).Color = [0,0,0];
p(3).Color = [0,0,0];


legend({'ACLF2, L: wait, M: wait ','ACLF2, L: accept, M: wait ','ACLF2, L: accept, M: accept '},FontSize=18)
xlabel('Day',FontSize=16) 
ylabel('D: Difference in Lifetime Total Expected Value ($)',FontSize=16) 
ylim([-5000 90000])
title('rr = 0.7, o = 70%',FontSize=18 )


%% over state
data =[
51244.8	90414.8	121248	136759	128928	95759.3	54187.8	23251.4	-2000	-2000	-2000	-2000	-2000	-2000;
90240.2	127508	98490.7	19461.5	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000;
104851	134884	27856.4	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000	-2000;

];
x = 1:1:3;
%2,4,6
y11 = data(1:3,1);
y21 = data(1:3,2);
y31 = data(1,3);
y32 = data(1:3,3);


figure
p = plot(x(1:3),y11,x(1:3),y21,x(1),y31,x(1:3),y32,'--');

set(gca,"FontSize",18)
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(1).Marker = '>';
p(1).MarkerIndices = 1:1:length(y11);
p(1).MarkerSize = 10;
p(2).Marker = '>';
p(2).MarkerIndices = 1:1:length(y21);
p(2).MarkerSize = 10;
p(3).Marker = '>';
p(3).MarkerIndices = [1];
p(3).MarkerSize = 10;

p(4).LineWidth = 3;

p(4).Marker = '+';
p(4).MarkerIndices = 2:1:length(y32);
p(4).MarkerSize = 10;


p(1).MarkerEdgeColor = [0,0,1];
p(2).MarkerEdgeColor = [1,0,0];
p(3).MarkerEdgeColor = [0,0,0];
p(4).MarkerEdgeColor = [0,0,0];

p(1).Color = [0,0,1];
p(2).Color = [1,0,0];
p(3).Color = [0,0,0];
p(4).Color = [0,0,0];

legend({'Day 2, L: wait, M: wait ','Day 4, L: wait, M: wait ', ...
    'Day 6, L: wait, M: wait ','Day 6, L: accept, M: wait'},FontSize=12,location = 'southwest')

xlabel('State',FontSize=16) 
ylabel('D: Difference in Lifetime Total Expected Value ($)',FontSize=14) 
set(gca,'xtick',1:3);
set(gca,'xticklabel',{'ACLF2', 'ACLF=3OF', 'ACLF>3OF'},'fontsize',14)
title('rr = 0.5, o = 70%',FontSize=16)
%% CKD
data = [
-199.371	-165.23	-130.487	-94.8354	-57.9444	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20;
-166.185	-131.305	-95.5395	-58.5587	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20	-20;
];
x = 45:1:64;
y11 = data(1,1:5);
y12 = data(1,5:end);
y21 = data(2,1:4);
y22 = data(2,4:end);

figure

p = plot(x(1:5),y11,x(5:end),y12,'--', ...
    x(1:4),y21,x(4:end),y22,'--');
set(gca,"FontSize",18)
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
p(4).LineWidth = 3;
p(1).Marker = '>';
p(1).MarkerIndices = [1,2,3,4,5,6];
p(1).MarkerSize = 10;
p(2).Marker = '+';
p(2).MarkerIndices = 2:1:length(y12);
p(2).MarkerSize = 10;
p(3).Marker = '>';
p(3).MarkerIndices = [1,2,3,4,5];
p(3).MarkerSize = 10;
p(4).Marker = '+';
p(4).MarkerIndices = 2:1:length(y22);
p(4).MarkerSize = 10;


p(1).MarkerEdgeColor = [0,0,1];
p(2).MarkerEdgeColor = [0,0,1];
p(3).MarkerEdgeColor = [1,0,0];
p(4).MarkerEdgeColor = [1,0,0];

p(1).Color = [0,0,1];
p(2).Color = [0,0,1];
p(3).Color = [1,0,0];
p(4).Color = [1,0,0];


legend({'CKD Stage 1, L: wait, M: wait ','CKD Stage 1, L: accept, M: accept ', ...
    'CKD Stage 2, L: wait, M: wait ','CKD Stage 2, L: accept, M: accept '},FontSize=18,location = 'southeast')
xlabel('Age',FontSize=16) 
ylabel('D: Difference in Lifetime Total Expected Value ($)',FontSize=16) 

%title('rr = 0.8, o = 70%',FontSize=18 )
