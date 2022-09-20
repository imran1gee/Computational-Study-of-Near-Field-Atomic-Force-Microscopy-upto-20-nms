clear all
tic
% x and z are real space variables
x = -20:40/(69):20;
zl = -2:22/(69):20;
%t1, t2, and t3 are integration variables
t1 = -99:99/(999):-0.001;
t2 = 0:pi/(999):pi;
t3 = 0:99/(999):99;
%x0, z0 are the position of the line source
x0 = 1;
z0 = 2;
%beta, Dmu, Dme are current strength, relative permeability, and relative
%permittivity
beta = 1;
Dmu = 0.9;
Dme = 0.01;
b = -0.0;
dt1 = (t1(2)-t1(1));
dt2 = (t2(2)-t2(1));
dt3 = (t3(2)-t3(1));
%The incident electric field Calculation
%%
count = 0;
for i = 1:1:length(x)
    for j = 1:1:length(zl)
        %Non oscillatory integration
        ai1 = 0;
        for k1 = 1:1:length(t1)
            ai1= ai1+1*exp(sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sinh(t1(k1)));
        end
        AI1(i,j) = dt1*ai1;
        ai3 = 0;
        for k3 = 1:1:length(t3)
            ai3= ai3+1*exp(-sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sinh(t3(k3)));
        end
        AI3(i,j) = dt3*ai3;
        %Oscillatory Integration
        ai2 = 0;
        for k2 = 1:1:length(t2)
            ai2= ai2+1*exp(-1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sin(t2(k2)));
        end
        AI2(i,j) = (-1i)*dt2*ai2;
        AI(i,j) = (beta/(1i*4*pi))*(1*AI1(i,j)+1*AI2(i,j)+1*AI3(i,j));
        count = count + 1;
        R(count) = sqrt((z0+zl(j))^2+(x0+x(i))^2);
        ANI(count) = AI(i,j);
    end
end
%The Reflected Electric Field Calculations
%%
clear AI1 AI2 AI3 AI
count = 0;
for i = 1:1:length(x)
    for j = 1:1:length(zl)
        %Non oscillatory integration
        ai1 = 0;
        for k1 = 1:1:length(t1)
            ai1= ai1+1*exp(sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sinh(t1(k1)))...
                *(sinh(t1(k1))-Dmu*sqrt(cosh(t1(k1))^2-1/(Dmu*Dme)))/(sinh(t1(k1))+Dmu*sqrt(cosh(t1(k1))^2-1/(Dmu*Dme)));
        end
        AI1(i,j) = dt1*ai1;
        ai3 = 0;
        for k3 = 1:1:length(t3)
            ai3= ai3+1*exp(-sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sinh(t3(k3)))...
                *((sinh(t3(k3))-Dmu*sqrt(cosh(t3(k3))^2-1/(Dmu*Dme)))/(sinh(t3(k3))+Dmu*sqrt(cosh(t3(k3))^2-1/(Dmu*Dme))));
        end
        AI3(i,j) = dt3*ai3;
        %Oscillatory integration
        ai2 = 0;
        for k2 = 1:1:length(t2)
            ai2= ai2+1*exp(-1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sin(t2(k2)))...
                *((1i*sin(t2(k2))-Dmu*sqrt(cos(t2(k2))^2-1/(Dmu*Dme)))/(1i*sin(t2(k2))+Dmu*sqrt(cos(t2(k2))^2-1/(Dmu*Dme))));
        end
        AI2(i,j) = (-1i)*dt2*ai2;
        AI(i,j) = (beta/(1i*4*pi))*(1*AI1(i,j)+1*AI2(i,j)+1*AI3(i,j));
        count = count + 1;
        ANR(count) = AI(i,j);
    end
end
%The Transmission Section
%%
clear AI1 AI2 AI3 AI
count = 0;
for i = 1:1:length(x)
    for j = 1:1:length(zl)
        %Non oscillatory integration
        ai1 = 0;
        for k1 = 1:1:length(t1)
            ai1= ai1+1*exp(sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sinh(t1(k1)))...
                *(2*Dmu*sqrt(cosh(t1(k1))^2-1/(Dmu*Dme)))/(sinh(t1(k1)+0.5)+Dmu*sqrt(cosh(t1(k1))^2-1/(Dmu*Dme)));
        end
        AI(i,j) = dt1*ai1;
        ai3 = 0;
        for k3 = 1:1:length(t3)
            ai3= ai3+1*exp(-sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sinh(t3(k3)))...
                *(2*Dmu*sqrt(cosh(t3(k3))^2-1/(Dmu*Dme)))/(sinh(t3(k3)+0.5)+Dmu*sqrt(cosh(t3(k3))^2-1/(Dmu*Dme)));
        end
        AI3(i,j) = dt3*ai3;
        %Oscillatory integration
        ai2 = 0;
        for k2 = 1:1:length(t2)
            ai2= ai2+1*exp(-1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sin(t2(k2)))...
                *(2*1i*Dmu*sqrt(-cos(t2(k2))^2+1/(Dmu*Dme))/(+0.00000001*(1-1i)+1i*sin(t2(k2))+1i*Dmu*sqrt(-cos(t2(k2))^2+1/(Dmu*Dme))));
        end
        AI2(i,j) = (-1i)*dt2*ai2;


        count = count + 1;
        R(count) = sqrt((z0+zl(j))^2+(x0+x(i))^2);
        ANT(count) = (beta/(1i*4*pi))*(1*AI(i,j)+1*AI3(i,j)+1*AI2(i,j));
    end
end

figure('DefaultAxesFontSize',24)

plot(R,real(ANI),'or',R,real(ANR),'ok',R,real(ANT),'ob',R,real(ANR+ANT),'.g')
xlabel('rk')
ylabel('Re[E[N/C]]')
% nexttile
% plot(R,imag(AN1),'or',R,imag(AN2),'ok')
% xlabel('rk')
% ylabel('Im[H_{0}^{1}]')
 legend('E^{in}','E^{re}','E^{tr}','E^{re}+E^{tr}')

toc