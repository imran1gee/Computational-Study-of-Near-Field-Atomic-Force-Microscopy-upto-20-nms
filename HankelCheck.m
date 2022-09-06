clear all
tic
x = -2:4/(69):2;
zl = -2:2/(69):0;
zr = 0:2/(69):2;
t = -39:78/(999):39;
x0 = 1;
z0 = 2;
beta = 1;
Dmu = 1;
Dep = 1;
b = 0.05;
dt = (t(2)-t(1));
for i = 1:1:length(x)
    for j = 1:1:length(zl)
        ai = 0;
        for k = 1:1:length(t)
            ai= ai+(beta/(4*pi))*1*exp(1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*(1+1i*b)*cosh(t(k)));
        end
        AI1(i,j) = dt*ai;
        AI2(i,j) = (1i)*(beta/4)*besselh(0,1,sqrt((x(i)+x0)^2+(zl(j)+z0)^2));
    end
end


figure
tiledlayout(2,2)
nexttile
contourf(x,zl,imag(AI1)')
xlabel('x')
ylabel('z')
nexttile
contourf(x,zl,real(AI1)')
xlabel('x')
ylabel('z')

nexttile
contourf(x,zl,imag(AI2)')
xlabel('x')
ylabel('z')
nexttile
contourf(x,zl,real(AI2)')
xlabel('x')
ylabel('z')
set(findall(gcf,'type','axes'),'fontsize',24)


% tiledlayout(2,2)
% nexttile
% mesh(x,zl,imag(AI1)')
% xlabel('x')
% ylabel('z')
% zlabel('Im[E^{in}_{y}]')
% nexttile
% mesh(x,zl,real(AI1)')
% xlabel('x')
% ylabel('z')
% zlabel('Re[E^{in}_{y}]')
% nexttile
% mesh(x,zl,imag(AI2)')
% xlabel('x')
% ylabel('z')
% zlabel('Im[(i/4)H^{1}_{0}]')
% nexttile
% mesh(x,zl,real(AI2)')
% xlabel('x')
% ylabel('z')
% zlabel('Re[(i/4)H^{1}_{0}]')
% set(findall(gcf,'type','axes'),'fontsize',24)
toc