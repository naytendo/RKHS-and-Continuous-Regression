load('HumanTreadmillData.mat')

hip = markerData(:,1:3);
knee = markerData(:,4:6);
ankle = markerData(:,7:9);
theta1 = zeros(1,length(hip));
theta2 = zeros(1,length(hip));

for ii = 1:length(theta1)
    [theta1(ii),theta2(ii)] = getJointAngles(hip(ii,:),knee(ii,:),ankle(ii,:));
end
% figure()
% plot(xSamp(:,1),xSamp(:,2));
m0 = 9000;
mEnd = length(theta1);
xSamp = [theta1(m0:mEnd)',theta2(m0:mEnd)']; %input data
tData = 0:1/150:(length(xSamp)-1)/150;
ySamp = Output_Q(m0:end,2); %output data


maxGrid = [65,65];
minGrid = [-35,-5];


centers = zeros(2,length(xSamp));
N = 1;
centers(:,N) = xSamp(1,:)';
separation = 10;
indexes = zeros(1,length(xSamp));
indexes(1) = 1;

for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:N
        if norm((centers(:,jj)-xSamp(mm,:)')) > separation
            check = check +1;
        end
    end
    if check == N
        centers(:,N+1) = xSamp(mm,:)';
        N = N+1;
        indexes(N) = mm;
 
        
    end
end
% centers = xSamp(:,1:3:end)';
centers = centers(:,1:N-1)';

indexes = indexes(1:N);

tNi = tData(indexes);
beta = 10;
K = zeros(N-1,N-1);
figure()

hold on
plot(xSamp(:,1),xSamp(:,2))
plot(centers(:,1),centers(:,2),'ko','linewidth',3)
grid on
title('Location of centers')
legend('$\phi(t)$','$\Xi_N$','interpreter','latex')


type = 'matern32';

for pp = 1:N-1
    for jj = 1:N-1
        K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
    end
end
cond(K)




%%
w0 = tNi(2:end)-tNi(1:end-1);
w = [w0(1),2*w0(2:end-1),w0(end)]/2;
W = diag(w);


gamma = 1e-3;
figstr = sprintf('%g',tend);
y = ySamp(indexes(1:end-1));

alpha = (K*W*K + gamma*K)\(K*W*y);

regressionCoefs = alpha';
res = 0.5;
x1Range = minGrid(1):res:maxGrid(1);
x2Range = minGrid(2):res:maxGrid(2);

[X1,X2] = meshgrid(x1Range,x2Range);





kernEstimate = zeros(size(X1));

for ii = 1:size(kernEstimate,1)
    for jj = 1:size(kernEstimate,2)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),[X1(ii,jj),X2(ii,jj)],beta);
        end
        kernEstimate(ii,jj) = regressionCoefs*kernVector;
    end
end
figure()

grid on
hold on

surf(X1,X2,kernEstimate,'FaceAlpha',0.5,'EdgeColor','none')
spacing = 3;  % play around so it fits the size of your data set
for i = 1 : spacing : length(X1(:,1))
    plot3(X1(i,:), X2(i,:), kernEstimate(i,:),'-','color',[0.5 0.5 0.5]);
end

for i = 1 : spacing : length(X1(1,:))
    plot3(X1(:,i), X2(:,i), kernEstimate(:,i),'-','color',[0.5 0.5 0.5]);  
end

plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp(:,1)),1),'color',[0 0 0], 'LineWidth', 0.25)
plot3(xSamp(:,1),xSamp(:,2),ySamp,'r', 'LineWidth', 0.25)



view(-8,45)
ax = gca;
font_sz_ticks = 28;
font_sz_labels = 48;
zlim([0 300])
ax.FontSize = font_sz_ticks;
xlabel('$\theta^{(1)}$','interpreter','latex','fontsize',font_sz_labels)
ylabel('$\theta^{(2)}$','fontsize',font_sz_labels,'interpreter','latex')
zlabel('$\hat{g}_N$','fontsize',36,'interpreter','latex')
text(40,10,0,'$\phi(t)$','fontsize',36,'interpreter','latex','color',[0 0 0])
% text(-30,30,450,strcat('$\gamma = $',figstr),'fontsize',font_sz_labels,'interpreter','latex')

set(gcf,'Position',[100 100 800 600])
saveas(gcf,'kernelEstimateMocapApprox1','epsc')