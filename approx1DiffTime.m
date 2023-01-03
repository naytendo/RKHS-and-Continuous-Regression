for tend = [10 20 50 100 200 500]
phi0 = [1;1;1];
tspan = [0 tend];
maxGrid = [30,30];
minGrid = [-30,-30];
[centers,xSamp,tNi] = generateCenters('along orbit',1,phi0,tspan);

%%% making quadrature points same as centers 
%%% (future results could make these different)
quadPts = centers;
tQuad = tNi;

G = @(x,y)(-5*sin(1/10*x)+1/2000*(x+y).^3+200);
% figure()
% plot(centers(:,1),centers(:,2),'o')
% title('centers in domain')
centers = centers(1:end-1,:);
quadPts = quadPts(1:end-1,:);
beta = 4;
M = length(centers);
Nw = length(quadPts);
K = zeros(M,Nw);
type = 'matern32';
    for jj = 1:M
        for pp = 1:Nw
            K(jj,pp) = kernel(type,centers(jj,:),quadPts(pp,:),beta);
        end
    end
    
 cond(K)
 
%%
w0 = tQuad(2:end)-tQuad(1:end-1);
w = [w0(1),2*w0(2:end-1)',w0(end)]/2;

W = diag(w);

gamma = 0.1;
figstr = sprintf('%g',tend);

y = zeros(length(quadPts),1);
for ii = 1:length(quadPts)
    y(ii) = G(quadPts(ii,1),quadPts(ii,2));
end

alpha = (K*W*K + gamma*(K))\(K*W*y);


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

surf(X1,X2,G(X1,X2),'EdgeColor','none','FaceAlpha',0.01,'FaceColor','green');
view(11,30)



surf(X1,X2,kernEstimate,'FaceAlpha',0.5,'EdgeColor','none')
spacing2 = 3;  % play around so it fits the size of your data set
% 
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


plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp(:,1)),1),'k', 'LineWidth',1)
ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G(xSamp(ii,1),xSamp(ii,2));
end
plot3(xSamp(:,1),xSamp(:,2),ySamp,'r', 'LineWidth', 0.1)


xlim([-30,30])
ylim([-30,30])
ax = gca;
font_sz_ticks = 28;
font_sz_labels = 48;
zticks(200:200:400);
ax.FontSize = font_sz_ticks;
xlabel('$x$','interpreter','latex','fontsize',font_sz_labels)
ylabel('$y$','fontsize',font_sz_labels,'interpreter','latex')
zlabel('$\hat{g}_N$','fontsize',36,'interpreter','latex')
text(10,-10,-12,'$\phi(t)$','fontsize',36,'interpreter','latex','color',[0 0 0])
% text(-30,30,450,strcat('$\gamma = $',figstr),'fontsize',font_sz_labels,'interpreter','latex')
text(-20,20,370,strcat('After  $',figstr, '$ seconds'),'fontsize',32,'interpreter','latex')
text(-20,25,270,'$G$','fontsize',font_sz_labels,'interpreter','latex','color',gridColor)
set(gcf,'Position',[100 100 800 600])


saveas(gcf,strcat('estimate',figstr,'sec-Method1'),'epsc')
end