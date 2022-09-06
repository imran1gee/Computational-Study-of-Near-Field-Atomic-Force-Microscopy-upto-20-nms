clear all
tic
x = -2:4/(69):2;% This is the direction of Fourier transformation 
zl = -2:2/(69):0;% The left side of the material
zr = 0:2/(69):2;% The right side of the materials
t = -38:78/(999):38; % The integral over kparallel
dt = (t(2)-t(1));
x0 = 1; % The location of the source
z0 = 2;
beta = 1;
Dmu = 0.9; % The permeability
Dep = 0.5; % The permitivity
b = 0.05; % The convergence factor
% The incident electric field AI
for i = 1:1:length(x)
    for j = 1:1:length(zl)
        ai = 0;
        % The integration is replaced by sum
        for k = 1:1:length(t)
            ai= ai+(beta/(4*pi))*1*exp(1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*(1+1i*b)*cosh(t(k)));
        end
        AI(i,j) = dt*ai;
    end
end
% The reflected electric field AR
for i = 1:1:length(x)
    for j = 1:1:length(zl)
        ar = 0;
        % The integration is replaced by sum
        for k = 1:1:length(t)
            ar= ar+(beta/(4*pi))*((cosh(t(k))-Dmu*sqrt((1/(Dmu*Dep))+sinh(t(k))^2))/(cosh(t(k))+Dmu*sqrt((1/(Dmu*Dep))+sinh(t(k))^2)))*1*exp(1i*sqrt((x(i)+x0)^2+(-zl(j)+z0)^2)*(1+1i*b)*cosh(t(k)));
        end
        AR(i,j) = dt*ar;
    end
end
% The transmitted electric field AT
for i = 1:1:length(x)
    for j = 1:1:length(zr)
        at = 0;
        % The integration is replaced by sum
        for k = 1:1:length(t)
            at= at+(beta/(4*pi))*2*cosh(t(k))/(cosh(t(k)+Dmu*sqrt((1/(Dmu*Dep))+sinh(t(k))^2)))*1*exp(1i*sqrt((x(i)+x0)^2+(zr(j)+z0)^2)*(1+1i*b)*sqrt(1+cosh(t(k))^2)*(1+...
                sqrt((1/(Dmu*Dep*sinh(t(k))^2))+1)));
        end
        AT(i,j) = dt*at;
    end
end
% Plotting imaginary and real part of incident, reflected and transmitted
% electric fields
tiledlayout(3,2)
nexttile
contourf(x,zl,imag(AI)')
nexttile
contourf(x,zl,real(AI)')
nexttile
contourf(x,zl,imag(AR)')
nexttile
contourf(x,zl,real(AR)')
nexttile
contourf(x,zr,imag(AT)')
nexttile
contourf(x,zr,real(AT)')
set(findall(gcf,'type','axes'),'fontsize',24)
toc