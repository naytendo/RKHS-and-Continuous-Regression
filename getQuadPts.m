function [quadPts,xSamp,tQuad] = getQuadPts(phi0,tspan,dt)
        %% Collect quadPts from orbit
            [t,phi1] = ode45(@lorenzSys,tspan,phi0);
            xSamp = phi1(:,1:2);

            quadPts = zeros(length(xSamp),2);
            tQuad = zeros(length(xSamp),1);
            M = 1;
            quadPts(M,:) = xSamp(1,:);
            tA = t(1);
            for mm = 1:length(xSamp)
                if abs(t(mm)-tA) > dt
                    quadPts(M+1,:) = xSamp(mm,:);
                    tQuad(M+1) = t(mm);
                    M = M+1;
                    tA = t(mm);
                end
            end
            quadPts = quadPts(1:M,:);
            tQuad = tQuad(1:M);
    
end