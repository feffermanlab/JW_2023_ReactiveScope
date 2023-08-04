clear


save_figures = 1; %0 for true
FigName = 'SigmoidSwitch';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRESS PARAMETERS
%define stress schedule 
numberofstressors=0;

%Exponential Stressor
%one_stressor_event=@(t, s, mu, sigma) s*exp(-(((t-mu)/sigma).^2));

%Gate Function Stressor
%one_stressor_event=@(t, s, mu, sigma) s*((t>mu)&(t<mu+sigma));

%Spliney Function Stresso
alpha=3;
one_stressor_event=@(t, s, mu, sigma) s*((mu < t & t<= mu+sigma).*((1/sigma).*(t-mu))+(mu+sigma<t & t<mu+alpha*sigma).*((1/(sigma-alpha*sigma)).*(t-mu-sigma)+1));
   %stressorinfo.s=[0,5,5]; %intensity no early stress
    %stressorinfo.s=[5,5,5]; %intensity yes early stress
     %stressorinfo.s=[5,2,6]; %intensity stress graphic
     % stressorinfo.s=[0,0,0]; %intensity psysiological mediator graphic
     stressorinfo.s= [];%[.4 .4];


%stressorinfo.sigma=[.75,.75,.75]; %duration 
 stressorinfo.sigma = [];%[.1 .1]; %for stress graphic
 %stressorinfo.mu=[14,38,57]; %time
stressorinfo.mu=[];%[0.2 .6] ;
%stressorinfo.sigma=.25*ones(1,numberofstressors);
    %come back and automate this
allstress = @(t) 0;
for i=1:numberofstressors
    allstress=@(t) allstress(t)+one_stressor_event(t,stressorinfo.s(i),stressorinfo.mu(i),stressorinfo.sigma(i));
end
%allstress=@(t) one_stressor_event(t,stressorinfo.s(1), stressorinfo.mu(1), stressorinfo.sigma(1)) +one_stressor_event(t,stressorinfo.s(2), stressorinfo.mu(2), stressorinfo.sigma(2))+one_stressor_event(t,stressorinfo.s(3), stressorinfo.mu(3), stressorinfo.sigma(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define level of physiological medaitor 
mediatorinfo.amp.seasonal=.2;%.2; %2
mediatorinfo.amp.circadian=.05;%.1; %1
years = 50/365;
mediatorinfo.per.seasonal = years*(2*pi);
mediatorinfo.per.circadian = 2*50*pi;%2*pi;
mediatorinfo.min=0.4;
mediatorinfo.active=0.05;

%define amount of decrease
r1 = .8;
r2 = .8;
r3 = .8;


 physiofun=@(t) mediatorinfo.amp.seasonal*sin(mediatorinfo.per.seasonal*(t-.5))+mediatorinfo.amp.circadian*sin(mediatorinfo.per.circadian*t)+allstress(t)+mediatorinfo.min;
 predictivescopefun=@(t) mediatorinfo.amp.seasonal*sin(mediatorinfo.per.seasonal*(t-.5))+mediatorinfo.amp.circadian+mediatorinfo.min+mediatorinfo.active;


posfun=@(x)max(0,sign(x));
pas.stressfun=allstress;
pas.physiofun=physiofun;
pas.predictivefun=predictivescopefun;
pas.posfun=posfun;
pas.r1=r1;
pas.r2 = r2;
pas.r3 = r3;

Tmax=1;
inc = Tmax/600;
Tfine=0:inc:Tmax;
opts = odeset('Events', @myEvent,'RelTol',10^-5, 'MaxStep', .01);%
options = odeset('Events',@(t,m, T_amb, p_amb,comp, M_sp, M_amb, Y_sp_inf, visk_amb, rho_amb,omega_amb_x, omega_amb_y, omega_amb_z, K) StopEvent(t,m, T_amb, p_amb,comp, M_sp, M_amb, Y_sp_inf, visk_amb, rho_amb,omega_amb_x, omega_amb_y, omega_amb_z, K),'Refine', 1);

%Initial conditions are reactive scope then max threshold
M0=1;
[T,Y]=ode45(@stress,[0 Tmax],[M0,M0],opts,pas);
figure('Units','inches','Position',[0.13,0.112857142857143,7.75,5.75],'PaperPositionMode','auto');
fill([T', Tmax,0], [Y(:,1)',[1.06*M0,1.06*M0]],[231/255, 161/255, 161/255])


 hold on

fill([T', fliplr(T')], [predictivescopefun(T'),fliplr(Y(:,1)')],[1, 1, 170/255] )
%[255, 255, 170]
fill([0,Tmax, fliplr(T')], [0,0,fliplr(predictivescopefun(T'))],[193/255, 227/255, 181/255])

 h1=plot(T(1:5:end),Y(1:5:end,1),'k.', 'LineWidth',1.5,'MarkerSize',12);

 % plot(T,Y(:,2),'k-.', 'LineWidth',1.5)
 
h2=plot(Tfine, predictivescopefun(Tfine), 'k-.','LineWidth',1.5);
 %patch([T';fliplr(T)'],[(predictivescopefun(T')); fliplr(Y(:,1)')], 'g')

h3=plot(Tfine, physiofun(Tfine),'k-','LineWidth',1.5);

h4 = plot(T,Y(:,2),'k--','LineWidth',1.5);

 legend([h4 h1 h2 h3], {'$M(t)$', '$R(t)$', '$P(t)$', '$y(t)$'},'Location','northeast','Interpreter','latex','FontSize',14)
 ylim([0,M0*1.06])

 xlabel('Time', 'FontSize',24)
 ylabel('Mediator Level', 'FontSize',24)
 ax=gca;
 ax.FontSize=24;
text(.05*Tmax,1.03*M0,{'Homeostatic Overload'},'FontSize',16,'fontweight','bold')
text(.05*Tmax,1.1*(mediatorinfo.active+mediatorinfo.min+mediatorinfo.amp.circadian),{'Reactive Homeostasis'},'FontSize',16,'fontweight','bold')
text(.05*Tmax,.03*M0,{'Predictive Homeostasis'},'FontSize',16,'fontweight','bold')
% patch([T',fliplr(T)'],[(predictivescopefun(T')), fliplr(Y(:,1)')], 'g')


if save_figures==0
    GSName = [FigName, 'RSplot_gs.eps'];
    CName = [FigName, 'RSplot_color.eps'];
    print('-deps2',GSName)
    print('-depsc2',CName)
end

%  figure('Units','inches','Position',[0.13,0.112857142857143,7.75,5.75],'PaperPositionMode','auto');
%   ylim([0,max(stressorinfo.s)*1.06])
%   hold on
%  ax=gca;
%  ax.FontSize=24;
%  ax.Legend
% box on
%   ylabel('Stress Schedule','FontSize',24)
%  xlabel('Time','FontSize',24)
%  plot(T, allstress(T),'k-', 'LineWidth',1.5)
 

if save_figures==0
    StressName = [FigName, 'Event.eps'];
    print('-deps',StressName);

    TextName = [FigName, 'Variables.txt'];
    TextID = fopen(TextName, 'w');
    fprintf(TextID,'The stressor magnitudes (si) are [');
    fprintf(TextID,'%g ',stressorinfo.s);
    fprintf(TextID,']\n');
    fprintf(TextID,'The stressor durations (sigma) are [');
    fprintf(TextID,'%g ',stressorinfo.sigma);
    fprintf(TextID,']\n');
    fprintf(TextID,'The stressor times (mu) are [');
    fprintf(TextID,'%g ',stressorinfo.mu);
    fprintf(TextID,']\n');
    fprintf(TextID,'The r values are [');
    fprintf(TextID,'%f %f %f ',r1,r2,r3);
    fprintf(TextID,']\n');
    fprintf(TextID,'%s %f','The seasonal amplitude (as) is ',mediatorinfo.amp.seasonal); fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The circadian amplitude (ac) is ',mediatorinfo.amp.circadian);fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The years is ',years);fprintf(TextID,'\n'); 
    fprintf(TextID,'%s %f','The per seaonal is ',mediatorinfo.per.seasonal); fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The per circadian is ',mediatorinfo.per.circadian); fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The ymin is ',mediatorinfo.min); fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The yact is ',mediatorinfo.active); fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The Tmax is ',Tmax); fprintf(TextID,'\n');
    fprintf(TextID,'%s %f','The M0 is ',M0);
    fclose(TextID);
end
%   figure(3)
%   ylim([0,6.5])
%   hold on
%  ylabel('Stress Schedule')
%   xlabel('Time')
%  plot(T, allstress(T),'k-', 'LineWidth',1)
 
%%
function [value,isterminal,direction] = myEvent(t,y,~)
value = y(2)-1;     % Detect height = 0
isterminal = 0;   % Stop the integration
direction = -1;   % Negative direction only

end
