clear all

tic
% x and z are real space variables
x = 0;
zl = -0.01;

%t1, t2, and t3 are integration variables
t1 = -99:99/(999):-0.000;
t2 = 0:pi/(999):pi;
t3 = 0.001:99/(999):99;
%x0, z0 are the position of the line source
x0 = 0;
dt1 = (t1(2)-t1(1));
dt2 = (t2(2)-t2(1));
dt3 = (t3(2)-t3(1));



% The position of source z0
z0 = -25:25/599:0;
Om = 1; % Terahertz source frequency
mu = 1; % Magnetic permeability
Et = 1; % The amplitude of the terahertz source
sig = 1*10^6; % Conductivity of the infinite wire
RR = 0.01; % Radius of the wire
theta = 1*pi/2; % Angle of incidence for the terahertz source
kt = 1; % the wave vector of the terahertz source
Dme = 1/100;
Dmu = 1;
ratio = 1/(Dme*Dmu);
for i = 1:1:length(z0)
    Ein(i) = -(1/4)*besselh(0,2,0.01);% The weight of the incident field from the infinite line source
    %%
    %Non oscillatory integration
    ai1 = 0;
    for k1 = 1:1:length(t1)
        ai1= ai1+1*exp(sqrt((2*z0(i))^2)*sinh(t1(k1)))...
            *(abs(sinh(t1(k1)))-Dmu*sqrt(cosh(t1(k1))^2-1/(Dmu*Dme)))/(abs(sinh(t1(k1)))+Dmu*sqrt(cosh(t1(k1))^2-1/(Dmu*Dme)));
    end
    AI1(i) = dt1*ai1;
    ai3 = 0;
    for k3 = 1:1:length(t3)
        ai3= ai3+1*exp(-sqrt((2*z0(i))^2)*sinh(t3(k3)))...
            *((abs(sinh(t3(k3)))-Dmu*sqrt(cosh(t3(k3))^2-1/(Dmu*Dme)))/(abs(sinh(t3(k3)))+Dmu*sqrt(cosh(t3(k3))^2-1/(Dmu*Dme))));
    end
    AI3(i) = dt3*ai3;
    %Oscillatory integration
    ai2 = 0;
    for k2 = 1:1:length(t2)
        ai2= ai2+1*exp(-1i*sqrt((2*z0(i))^2)*sin(t2(k2)))...
            *((1*(sin(t2(k2)))-Dmu*sqrt(-cos(t2(k2))^2+1/(Dmu*Dme)))/(0.00001*1i+1*(sin(t2(k2)))+Dmu*sqrt(-cos(t2(k2))^2+1/(Dmu*Dme))));
    end
    AI2(i) = (-1i)*dt2*ai2;
    Er(i) = -(1i/(4*pi))*(1*AI1(i)+1*AI2(i)+1*AI3(i));
    alpha(i) = Ein(i)+Er(i); % The proportionality weight for left side field of the infinite line source
    R(i) = sqrt(z0(i)^2);
    I(i) = 1/(1+alpha(i));  % Selfconsistent current, Isc = It+Isc*alpha
end

%%
for i = 1:1:length(z0)
    alphaa(i) = -(1/4)*(besselh(0,2,0.01)-besselh(0,2,abs(2*z0(i)))); % The proportionality weight for self-consistent current
    Ia(i) =1/(1+alphaa(i));  % Selfconsistent field,
end




% semilogx((1/(2*pi))*R,abs(Esc),'-ob',(1/(2*pi))*R,(Esca),'-or','LineWidth',0.3)


save('Data_Fixed_Current_3.mat','R',"I","Ia")
figure('DefaultAxesFontSize',24)
% 
plot((1/(2*pi))*R,imag(I),'-ob',(1/(2*pi))*R,imag(Ia),'-or','LineWidth',0.3)
% xline(1,'--')
% xline(1.5,'--')
% plot((1/(2*pi))*R,real(I-Ia),'-ob','LineWidth',0.3)
% ylim([-100,100])
%     plot(Rr,imag((ANIS)),'-og',Rr,-imag((ANRS)),'-ok','LineWidth',1)
% yline(0)
% plot(Rt,imag(ANTS),'-ob',Rr,imag(ANIS+ANRS),'-oc','LineWidth',3)
% 
% %plot(sort(R),real(ANIS),'-r')
xlabel('rk')
ylabel('Re[E[N/C]]')
% % nexttile
% % plot(R,imag(AN1),'or',R,imag(AN2),'ok')
% % xlabel('rk')
% % ylabel('Im[H_{0}^{1}]')
 legend('E','E^{PEC}','Location','southeast')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters Self_Field_Near
toc

