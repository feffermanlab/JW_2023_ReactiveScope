function [dydt] = stress(t,y,struc)
dydt=zeros(2,1);
physio=struc.physiofun;
predictive=struc.predictivefun;
allstress = struc.stressfun;
r=struc.r;
%Senesence
%dydt(1)=-r*max(0,sign(physio(t)-predictive(t)))+r*max(0,sign(predictive(t)-physio(t)))*max(0,sign(y(2)-y(1)))+max(0,sign(y(1)-y(2)))*(-(t<.1)*(0)+((t>.1)&(t<.9))*((1-.6)/(.1-.9))-(t>.9)*(0));
%dydt(2)=-r*max(0,sign(allstress(t)))-(t<.1)*(0)+((t>.1)&(t<.9))*((1-.6)/(.1-.9))-(t>.9)*(0);
%No Senesence
dydt(1)=-r*sign(allstress(t))+r*(allstress(t)==0)*max(0,sign(y(2)-y(1)));
dydt(2)=-r*max(0,sign(physio(t)-y(1)));
end


%P=2*sin(t/2)+sin(5*t)+5+STRESS
%PH=3*sin(t/2)+1+5+1
%plot(2*sin(T/2)+sin(5*T)+5+@STRESS(T,mu,sigma)
%STRESS_T=