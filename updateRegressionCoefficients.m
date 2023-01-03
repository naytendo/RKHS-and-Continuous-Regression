function [dstatedt] = updateRegressionCoefficients(t,states,K,centers,type,beta,gamma,G)
    phi = states(1:3);
    alpha = states(4:end);
    sigma = 10;
    r = 28;
    b = 8/3;
%     G = @(x,y)(-10*sin(1/10*x)+1/2000*(y+x)^3+200);
    % Define states
    x = phi(1);
    y = phi(2);
    z = phi(3);

    dxdt = sigma*(y - x);
    dydt = x*(r- z)-y;
    dzdt = x*y -  b*z;
    phidot = [dxdt;dydt;dzdt];
    kVec = zeros(length(centers),1);    
    for cc = 1:length(centers)
        kVec(cc) = kernel(type,centers(cc,:),[x,y],beta);
    end
    gK = gamma*K;
    alphadot = gK\(-kVec*kVec'*alpha+kVec*G(x,y));
    dstatedt = [phidot;alphadot];
    t
end