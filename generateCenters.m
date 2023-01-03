function [centers,xSamp,tNi] = generateCenters(option,sep,phi0,tspan)
    switch option
        case 'along orbit'
        %% Collect centers from orbit
            

            [t,phi1] = ode45(@lorenzSys,tspan,phi0);
            xSamp = phi1(:,1:2);

            centers = zeros(length(xSamp),2);
            tNi = zeros(length(xSamp),1);
            M = 1;
            centers(M,:) = xSamp(1,:);
            
            for mm = 1:length(xSamp)
                check = 0;
                for jj = 1:M
                    if norm((centers(jj,:)-xSamp(mm,:))') > sep
                        check = check +1;
                    end
                end
                if check == M
                    centers(M+1,:) = xSamp(mm,:);
                    tNi(M+1) = t(mm);
                    M = M+1;
                end
            end
            centers = centers(1:M,:);
            tNi = tNi(1:M);
    %% Uniformly spaced centers
        case 'uniform grid'
            maxGrid = [30,30];
            minGrid = [-30,-30];
            spacing = sep;
            x1Range = minGrid(1):spacing:maxGrid(1);
            x2Range = minGrid(2):spacing:maxGrid(2);

            [X1,X2] = meshgrid(x1Range,x2Range);
            sz = size(X1);

            centers = zeros(sz(1)*sz(2),2);
            idx = 1;
            for ii = 1:sz(1)
                for jj = 1:sz(2)
                    centers(idx,:) = [X1(ii,jj),X2(ii,jj)];
                    idx = idx+1;
                end
            end

        otherwise 
   
        error('Generated centers option is not supported!');
    end
end