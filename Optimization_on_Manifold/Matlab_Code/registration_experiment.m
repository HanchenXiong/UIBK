shifts=[80,50,50];
angles=[0.2*pi,0.3*pi, 0.8*pi];
datagenerator('RedCup.txt',angles(1),angles(2),angles(3),shifts(1), shifts(2),shifts(3),0);

source=load('source.txt');
target=load('target.txt');

num_point=size(source,1);

h=figure;
set(gcf,'position',[400,200,1000,500]);               
subplot(1,2,1);
scatter3(target(:,1),target(:,2),target(:,3),10,'filled','black'); hold on;
scatter3(source(:,1),source(:,2),source(:,3),10,'red');
set(gca,'xticklabel',[],'xtick',[]);
set(gca,'yticklabel',[],'ytick',[]);
set(gca,'zticklabel',[],'ztick',[]);
axis equal; axis vis3d;

grid off;
view(60,40);
title('Before registration','fontsize',25);

      
target=target-repmat(shifts,num_point,1);

%[optimal_rotation,trace_1]=rotation_search(source,target,'N-SA');
[optimal_rotation,trace_2]=rotation_search(source,target,'SMC');
%[optimal_rotation,trace_3]=rotation_search(source,target,'ASMC');




subplot(1,2,2);
scatter3(target(:,1),target(:,2),target(:,3),10,'filled','black'); hold on;
scatter3(source(:,1),source(:,2),source(:,3),10,'red');
set(gca,'xticklabel',[],'xtick',[]);
set(gca,'yticklabel',[],'ytick',[]);
set(gca,'zticklabel',[],'ztick',[]);
axis equal; axis vis3d;
grid off;
view(60,40);
title('After registration','fontsize',25);

figure;
%plot(trace_1,'b'); hold on;
plot(trace_2,'r'); hold on;
%plot(trace_3,'g');
%legend('N-SA','SMC','ASMC');

