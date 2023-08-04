function [dydt] = stress(t,y,struc)
dydt=zeros(2,1);
physio=struc.physiofun;
predictive=struc.predictivefun;
r1=struc.r1;
r2 = struc.r2;
r3 = struc.r3;
sigmoid = @(x) 1/(1+exp(-50*x));
%Senesence
%Linear
% dydt(1)=-r2*max(0,sign(physio(t)-predictive(t)))+r3*max(0,sign(predictive(t)-physio(t)))*max(0,sign(y(2)-y(1)))+max(0,sign(y(1)-y(2)))*(-(t<.1)*(0)+((t>.1)&(t<.9))*((1-.6)/(.1-.9))-(t>.9)*(0));
% dydt(2)=-r1*max(0,sign(physio(t)-y(1)))-(t<.1)*(0)+((t>.1)&(t<.9))*((1-.6)/(.1-.9))-(t>.9)*(0);
%Sigmoid
% dydt(1)=-r2*max(0,sign(physio(t)-predictive(t)))+r3*max(0,sign(predictive(t)-physio(t)))*max(0,sign(y(2)-y(1)))+max(0,sign(y(1)-y(2)))*(-10*(1-.6)*exp(10*(.4-t)))/(1+exp(10*(.4-t)))^2;
% dydt(2)=-r1*max(0,sign(physio(t)-y(1)))-(10*(1-.6)*exp(10*(.4-t)))/(1+exp(10*(.4-t)))^2;
%No Senesence
%Heaviside
dydt(1)=-r2*max(0,sign(physio(t)-predictive(t)))+r3*max(0,sign(predictive(t)-physio(t)))*max(0,sign(y(2)-y(1)));
dydt(2)=-r1*max(0,sign(physio(t)-y(1)));
%Sigmoid
% dydt(1)=-r2*sigmoid(physio(t)-predictive(t))+r3*sigmoid(predictive(t)-physio(t))*max(0,sign(y(2)-y(1)));
% dydt(2)=-r1*sigmoid(physio(t)-y(1));

end


%P=2*sin(t/2)+sin(5*t)+5+STRESS
%PH=3*sin(t/2)+1+5+1
%plot(2*sin(T/2)+sin(5*T)+5+@STRESS(T,mu,sigma)
%STRESS_T=