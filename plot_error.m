close all
centers = generateCenters('along orbit',3);

figure()
plot(centers(:,1),centers(:,2),'o')
title('centers in domain')

beta = 4;
M = length(centers);
K = zeros(M,M);
type = 'matern32';
phi0 = [0;-5;1];
    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end
 
%%
tspanNew = [0 100];
gamma = 0.1;
figstr = sprintf('%g',tspanNew(end));
alpha0 = zeros(length(centers),1);
[t2,states2] = ode45(@(t,states)updateRegressionCoefficients(t,states,K,centers,type,beta,gamma),tspanNew,[phi0;alpha0]);
phiNew = states2(:,1:3);
alphaNew = states2(:,4:end);



figure()
title('evolution of the coefficients')
hold on
for kk = 1:length(centers)
    plot(t2,alphaNew(:,kk))
end
    
%%
maxGrid = [30,30];
minGrid = [-35,-35];
regressionCoefs = alphaNew(end,:);
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
% surf(X1,X2,(-10*sin(1/10*X1)+1/1000*(X2+X1).^3)+300,'EdgeColor','none','FaceAlpha',0.4,'FaceColor','green');
grid on
hold on
error = abs((-10*sin(1/10*X1)+1/2000*(X2+X1).^3)+300-kernEstimate);
surf(X1,X2,error,'FaceAlpha',0.5,'EdgeColor','none')
spacing = 4;  % play around so it fits the size of your data set
% for i = 1 : spacing : length(X1(:,1))
%     plot3(X1(:,i), X2(:,i), error(:,i),'-','color',[0.3 0.3 0.3]);
%     plot3(X1(i,:), X2(i,:), error(i,:),'-','color',[0.3 0.3 0.3]);
% end
% plot3(phiNew(:,1),phiNew(:,2),zeros(length(phiNew(:,1)),1),'k', 'LineWidth', 0.25)
G = @(x)(-10*sin(1/10*x(1))+1/1000*(x(2)+x(1))^3+300);

xSamp = phiNew(:,1:2);

ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G((xSamp(ii,:)));
end
plot3(xSamp(:,1), xSamp(:,2),zeros(length(xSamp(:,1)),1),'.', 'color','r', 'LineWidth', 2, 'MarkerSize', 3);
view(-20,28)
ax = gca;
font_sz_ticks = 28;
font_sz_labels = 48;
zticks(200:200:400);
ax.FontSize = font_sz_ticks;
xlabel('$x$','interpreter','latex','fontsize',font_sz_labels)
ylabel('$y$','fontsize',font_sz_labels,'interpreter','latex')
zlabel('$| G-\hat{g}_N |$','fontsize',36,'interpreter','latex')
% text(0,20,0,'$\phi(t)$','fontsize',font_sz_labels,'interpreter','latex','color','black')
% text(-20,20,400,'$G$','fontsize',font_sz,'interpreter','latex','color','green')
% text(10,-20,50,strcat('$\gamma = $',figstr),'fontsize',font_sz_labels,'interpreter','latex')
text(0,-20,20,strcat('After  $',figstr, '$ seconds'),'fontsize',28,'interpreter','latex')
set(gcf,'Position',[100 100 800 600])

cmap = colormap;

% colormap(1-colormap)
%