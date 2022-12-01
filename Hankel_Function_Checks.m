clear all
tic
%t1, t2, and t3 are integration variables
t1 = -99:99/(999):-0.000;
t2 = 0:pi/(999):pi;
t3 = 0.001:99/(999):99;
%x0, z0 are the position of the line source
dt1 = (t1(2)-t1(1));
dt2 = (t2(2)-t2(1));
dt3 = (t3(2)-t3(1));
z = -49.9:49.4/499:49.9;

%%Hankel Function Comparison 
%% Hankel Function Type-2
for i = 1:1:length(z)
    %Non oscillatory integration
    ai1 = 0;
    for k1 = 1:1:length(t1)
        ai1= ai1+1*exp(abs(z(i))*sinh(t1(k1)));
    end
    AI1(i) = dt1*ai1;
    ai3 = 0;
    for k3 = 1:1:length(t3)
        ai3= ai3+1*exp(-abs(z(i))*sinh(t3(k3)));
    end
    AI3(i) = dt3*ai3;
    %Oscillatory integration
    ai2 = 0;
    for k2 = 1:1:length(t2)
        ai2= ai2+1*exp(1i*abs(z(i))*sin(t2(k2)));
    end
    AI2(i) = (1i)*dt2*ai2;
    Bes1(i) = (-1i/4)*(1*AI1(i)+1*AI2(i)+1*AI3(i));
end

%% Hankel Function Type-2
for i = 1:1:length(z)
    %Non oscillatory integration
    ai1 = 0;
    for k1 = 1:1:length(t1)
        ai1= ai1+1*exp(abs(z(i))*sinh(t1(k1)));
    end
    AI1(i) = dt1*ai1;
    ai3 = 0;
    for k3 = 1:1:length(t3)
        ai3= ai3+1*exp(-abs(z(i))*sinh(t3(k3)));
    end
    AI3(i) = dt3*ai3;
    %Oscillatory integration
    ai2 = 0;
    for k2 = 1:1:length(t2)
        ai2= ai2+1*exp(-1i*abs(z(i))*sin(t2(k2)));
    end
    AI2(i) = (-1i)*dt2*ai2;
    Bes2(i) = (1i/4)*(1*AI1(i)+1*AI2(i)+1*AI3(i));
end

plot(z,imag(Bes2),'-ob',z,imag(besselh(0,2,abs(z))),'-or','LineWidth',0.3)

toc

