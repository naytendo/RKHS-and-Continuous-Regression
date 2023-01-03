function [alphadot] = updateMoCap(t,alpha,K,centers,type,beta,gamma,tData,thetas,y)
    tol = 1/300;
    if t <= 2/120
        k1 = 1;
        k2 = 10;
    else
        k1 = find(abs(tData- (t-5/120))<tol);
        k2 = find(abs(tData- (t+5/120))<tol);
    end
    theta1t = spline(tData(k1:k2),thetas(k1:k2,1),t);
    theta2t = spline(tData(k1:k2),thetas(k1:k2,2),t);
    yt = spline(tData(k1:k2),y(k1:k2),t);
%     theta1t = spline(tData,thetas(:,1),t);
%     theta2t = spline(tData,thetas(:,2),t);
%     yt = spline(tData,y,t);
    kVec = zeros(length(centers),1);    
    for cc = 1:length(centers)
        kVec(cc) = kernel(type,centers(cc,:),[theta1t,theta2t],beta);
    end
    gK = gamma*K;
    alphadot = gK\(-kVec*kVec'*alpha+kVec*yt);
    t
end