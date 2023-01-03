tspan = [0 100];
phi01 = [1;0;0];
phi02 = [1;2;0];
phi03 = [1;0;3];
[t1,phi1] = ode45(@lorenzSys,tspan,phi01);
[t2,phi2] = ode45(@lorenzSys,tspan,phi02);
[t3,phi3] = ode45(@lorenzSys,tspan,phi03);
plot3(phi1(:,1),phi1(:,2),phi1(:,3))
grid on
hold on
plot3(phi2(:,1),phi2(:,2),phi2(:,3))
plot3(phi3(:,1),phi3(:,2),phi3(:,3))

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex','fontsize',36)
ylabel('$y$','fontsize',36,'interpreter','latex')
zlabel('$z$','fontsize',36,'interpreter','latex')

view(28,22)

figure()
%%
grid on 
hold on 
plot(phi3(:,1),phi3(:,2))
plot(phi2(:,1),phi2(:,2))
plot(phi1(:,1),phi1(:,2))

set(gca,'FontSize',20)
xlabel('$x$','interpreter','latex','fontsize',36)
ylabel('$y$','fontsize',36,'interpreter','latex')
