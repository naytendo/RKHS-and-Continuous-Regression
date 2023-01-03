
[centers,~,~] = generateCenters('along orbit',3,[1;1;1],[0 500]);

% figure()
% plot(centers(:,1),centers(:,2),'o')
% title('centers in domain')
G = @(x,y)(-10*sin(1/10*x)+1/2000*(x+y).^3+200);
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
tend = 200;
tspanNew = [0 tend];

for gamma = [1e-4]
timestr = sprintf('%g',tspanNew(end));
gammaStr = sprintf('%g',gamma);
alpha0 = zeros(length(centers),1);
[t2,states2] = ode45(@(t,states)updateRegressionCoefficients(t,states,K,centers,type,beta,gamma,G),tspanNew,[phi0;alpha0]);
phiNew = states2(:,1:3);
alphaNew = states2(:,4:end);



% figure()
% title('evolution of the coefficients')
% hold on
% for kk = 1:length(centers)
%     plot(t2,alphaNew(:,kk))
% end
    
%%
maxGrid = [30,30];
minGrid = [-30,-30];
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
surf(X1,X2,G(X1,X2),'EdgeColor','none','FaceAlpha',0.01,'FaceColor','green');
grid on
hold on

surf(X1,X2,kernEstimate,'EdgeColor','none','FaceAlpha',0.6)
spacing2 = 3;  % play around so it fits the size of your data set

gridColor2 = [0.3 0.3 0.3];
for i = 1 : spacing2 : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), kernEstimate(:,i),'-','color',gridColor2);
    plot3(X1(i,:), X2(i,:), kernEstimate(i,:),'-','color',gridColor2);
end

gridColor = [0.1 0.1 0.1];
gridSpacing = 40;  % play around so it fits the size of your data set
for i = 1 : gridSpacing : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), G(X1(:,i),X2(:,i)),'-','color',gridColor);
    plot3(X1(i,:), X2(i,:), G(X1(i,:),X2(i,:)),'-','color',gridColor);
end

plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp(:,1)),1),'k', 'LineWidth', 0.25)

ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G(xSamp(ii,1),xSamp(ii,2));
end
plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r', 'LineWidth', 0.5, 'MarkerSize', 1);

view(11,30)
ax = gca;
font_sz_ticks = 28;
font_sz_labels = 48;
zticks(200:200:400);

ax.FontSize = font_sz_ticks;
xlabel('$x$','interpreter','latex','fontsize',font_sz_labels)
ylabel('$y$','fontsize',font_sz_labels,'interpreter','latex')
zlabel('$\hat{g}_N$','fontsize',36,'interpreter','latex')
text(10,-10,-15,'$\phi(t)$','fontsize',36,'interpreter','latex','color',[0 0 0])
% text(-30,30,450,strcat('$\gamma = $',figstr),'fontsize',font_sz_labels,'interpreter','latex')
text(0,20,370,strcat('$\gamma  =$',gammaStr),'fontsize',32,'interpreter','latex')
text(-30,25,260,'$G$','fontsize',font_sz_labels,'interpreter','latex','color',gridColor)
set(gcf,'Position',[100 100 800 600])

xlim([-30,30])
ylim([-30,30])
p1 = [30,10,20];                         % First Point
p2 = [17,7,0];                         % Second Point
dp = p2-p1;                         % Difference
% quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),'linewidth',2,'color',[0.8500 0.3250 0.0980])

if gamma < 1
    saveas(gcf,strcat('estimateGammaP',gammaStr(3:end),'-Method2'),'epsc')
else
    saveas(gcf,strcat('estimateGamma',gammaStr,'-Method2'),'epsc')
end

end
%