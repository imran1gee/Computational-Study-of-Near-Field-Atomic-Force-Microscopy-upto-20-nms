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
        %Oscillatory integration Hankel First Kind
        ai2 = 0;
        for k2 = 1:1:length(t2)
            ai2= ai2+1*exp(1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sin(t2(k2)));
        end
        %Hankel Second Kind
        ai2b = 0;
        for k2 = 1:1:length(t2)
            ai2b= ai2b+1*exp(-1i*sqrt((x(i)+x0)^2+(zl(j)+z0)^2)*sin(t2(k2)));
        end
        AI2(i,j) = (1i)*dt2*ai2;
        AI2b(i,j) = (-1i)*dt2*ai2b;
        AI(i,j) = (1/(1i*pi))*(1*AI1(i,j)+1*AI2(i,j)+1*AI3(i,j));
        AIb(i,j) = -(1/(1i*pi))*(1*AI1(i,j)+1*AI2b(i,j)+1*AI3(i,j));
        AAI(i,j) = besselh(0,1,sqrt((x(i)+x0)^2+(zl(j)+z0)^2));
        ABI(i,j) = besselh(0,2,sqrt((x(i)+x0)^2+(zl(j)+z0)^2));
        count = count + 1;
        R(count) = sqrt((z0+zl(j))^2+(x0+x(i))^2);
        AN1(count) = AI(i,j);
        AN1b(count) = AIb(i,j);
        AN2(count) = AAI(i,j);
        AN2b(count) = ABI(i,j);

    end
end

plot(R,real(AN1),'or',R,real(AN2),'ok')


% figure
% tiledlayout(1,2)
% nexttile
plot(R,real(AN1),'or',R,real(AN2),'ok')
xlabel('rk')
ylabel('Re[H_{0}^{1}]')
% nexttile
% plot(R,imag(AN1),'or',R,imag(AN2),'ok')
% xlabel('rk')
% ylabel('Im[H_{0}^{1}]')
legend('Integrated','H_{0}^{1}')
set(findall(gcf,'type','axes'),'fontsize',24)
toc