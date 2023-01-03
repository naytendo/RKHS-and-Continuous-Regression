tspan = 0:0.01:5;
phi0 = [0;0];
xSamp = zeros(length(tspan),2);
for ii = 1:length(tspan)
    xSamp(ii,1) = 15*sin(3*tspan(ii))+7*sin(7*tspan(ii));
    xSamp(ii,2) = 10*tspan(ii) ;
end

G = @(x)(6*sin(1/10*x(2))-1/2000*(x(2)+x(1))^2-1/2)+20;



ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G((xSamp(ii,:)));
end

depth = 4;
numChild = 2;
maxGrid = [40,60];
minGrid = [-30,-10];

approx = zeros(1,length(xSamp));

centers = zeros(length(xSamp),2);
M = 1;
centers(M,:) = xSamp(1,:);
separation = 25;
indexes = zeros(1,length(xSamp));
indexes(1) = 1;
for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:M
        if norm((centers(jj,:)-xSamp(mm,:))') > separation
            check = check +1;
        end
    end
    if check == M
        centers(M+1,:) = xSamp(mm,:);
        M = M+1;
        indexes(M) = mm;
    end
end
% centers = xSamp(:,1:3:end)';
centers = centers(1:M,:);
indexes = indexes(1:M);

beta = 5;
K = zeros(M,M);
Kf = zeros(M,M);
type = 'matern32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

regressionCoefs = ySamp(indexes)*pinv(K);
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

surf(X1,X2,(6*sin(1/10*X2)-1/2000*(X2+X1).^2-1/2)+20,'EdgeColor','none','FaceAlpha',0.4,'FaceColor','green');
grid on
hold on
% surf(X1,X2,kernEstimate,'FaceAlpha',0.1,'FaceColor','blue','EdgeColor','none')
% 
% spacing = 8;  % play around so it fits the size of your data set
% for i = 1 : spacing : length(X1(:,1))
%     plot3(X1(:,i), X2(:,i), kernEstimate(:,i),'-','color',[0.5 0.5 0.5]);
%     plot3(X1(i,:), X2(i,:), kernEstimate(i,:),'-','color',[0.5 0.5 0.5]);
% end
plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp),1),'k', 'LineWidth', 2)
plot3(xSamp(:,1),xSamp(:,2),ySamp,'LineWidth', 2,'color','red')

ylabel('$X$','fontsize',20,'interpreter','latex')
text(17,32,4,'$\Gamma_{\phi_0}$','fontsize',20,'interpreter','latex')
% text(-30,-10,4,'$\Pi_{\nu} G$','fontsize',20,'interpreter','latex','color','blue')
 text(-27,-9,19,'$G$','fontsize',20,'interpreter','latex','color','green')
view(75,36)

% plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r.', 'LineWidth', 1, 'MarkerSize', 5);
c_index = 3;
plot3(centers(c_index,1),centers(c_index,2),0,'ko','markerSize',10,'linewidth',2)
text(centers(c_index,1)-5,centers(c_index,2),28,'$y(\tau) = G(\phi(\tau))$','fontsize',20,'interpreter','latex','color','red')
text(centers(c_index,1)+15,centers(c_index,2)-10,2,'$\phi(\tau)$','fontsize',20,'interpreter','latex')
% text(centers(c_index,1)+10,centers(c_index,2)-10,2,'$\xi_{M,i}$','fontsize',20,'interpreter','latex')
plot3([centers(c_index,1) centers(c_index,1)],[centers(c_index,2) centers(c_index,2)],[0 ySamp(indexes(c_index))],'k-.','linewidth',2)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'ZTickLabel', [])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])

