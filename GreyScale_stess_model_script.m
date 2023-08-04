clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRESS PARAMETERS
%define stress schedule 
numberofstressors=2;

%Exponential Stressor
%one_stressor_event=@(t, s, mu, sigma) s*exp(-(((t-mu)/sigma).^2));

%Gate Function Stressor
one_stressor_event=@(t, s, mu, sigma) s*((t>mu)&(t<mu+sigma));

%Spliney Function Stresso
%one_stressor_event=@(t, s, mu, sigma) s*((mu < t & t<= mu+sigma).*((1/sigma).*(t-mu))+(mu+sigma<t & t<mu+3*sigma).*((1/(sigma-3*sigma)).*(t-mu-sigma)+1));
   %stressorinfo.s=[0,5,5]; %intensity no early stress
    %stressorinfo.s=[5,5,5]; %intensity yes early stress
     %stressorinfo.s=[5,2,6]; %intensity stress graphic
     % stressorinfo.s=[0,0,0]; %intensity psysiological mediator graphic
     stressorinfo.s=[.4, .4];


%stressorinfo.sigma=[.75,.75,.75]; %duration 
 stressorinfo.sigma=[.1, .1]; %for stress graphic
 %stressorinfo.mu=[14,38,57]; %time
stressorinfo.mu=[.2, .6];
%stressorinfo.sigma=.25*ones(1,numberofstressors);
    %come back and automate this
allstress = @(t) 0;
for i=1:numberofstressors
    allstress=@(t) allstress(t)+one_stressor_event(t,stressorinfo.s(i),stressorinfo.mu(i),stressorinfo.sigma(i));
end
%allstress=@(t) one_stressor_event(t,stressorinfo.s(1), stressorinfo.mu(1), stressorinfo.sigma(1)) +one_stressor_event(t,stressorinfo.s(2), stressorinfo.mu(2), stressorinfo.sigma(2))+one_stressor_event(t,stressorinfo.s(3), stressorinfo.mu(3), stressorinfo.sigma(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define level of physiological medaitor 
mediatorinfo.amp.seasonal=0.2; %2
mediatorinfo.amp.circadian=0.05; %1
years = 1;
mediatorinfo.per.seasonal=years*(2*pi);
mediatorinfo.per.circadian=years*50*2*pi;
mediatorinfo.min=.3;
mediatorinfo.active=0.01;

%define amount of decrease
r1 = .6;
r2 = 0.6;
r3 = 0.6;


 physiofun=@(t) mediatorinfo.amp.seasonal*sin(mediatorinfo.per.seasonal*t)+mediatorinfo.amp.circadian*sin(mediatorinfo.per.circadian*t)+allstress(t)+mediatorinfo.min;
 predictivescopefun=@(t) mediatorinfo.amp.seasonal*sin(mediatorinfo.per.seasonal*t)+mediatorinfo.amp.circadian+mediatorinfo.min+mediatorinfo.active;


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
figure1 = figure;
colormap(gray);

axes1 = axes('Parent',figure1);
hold(axes1,'on');


% Create multiple lines using matrix input to plot
plot1 = plot(T,Y,'LineWidth',1.5,'Color',[0 0 0],'Parent',axes1);
set(plot1(1),'DisplayName','$R(t)$','Marker','.','LineStyle','none');
set(plot1(2),'DisplayName','$M(t)$','LineStyle','--');

% Create patch
patch([T',fliplr(T')],[Y(:,1)',1.06*M0*ones(1,length(T))],[0.501960784313725 0.501960784313725 0.501960784313725]);

patch([T',fliplr(T')],[Y(:,1)',predictivescopefun(fliplr(T'))],[0.8 0.8 0.8]);

patch([T',fliplr(T')],[predictivescopefun(T'),zeros(1,length(T))],[0.941176470588235 0.941176470588235 0.941176470588235]);
 h1=plot(T,Y(:,1),'k.', 'LineWidth',1.5);

 % plot(T,Y(:,2),'k-.', 'LineWidth',1.5)
 
h2=plot(Tfine, predictivescopefun(Tfine), 'k-.','LineWidth',1.5);
 %patch([T';fliplr(T)'],[(predictivescopefun(T')); fliplr(Y(:,1)')], 'g')

h3=plot(Tfine, physiofun(Tfine),'k-','LineWidth',1.5);

h4 = plot(T,Y(:,2),'k--','LineWidth',1.5);

 legend([h4 h1 h2 h3], {'$M(t)$', '$R(t)$', '$P(t)$', '$y(t)$'},'Location','northeast','Interpreter','latex','FontSize',14)
 ylim([0,1.06])

 xlabel('Time', 'FontSize',16)
 ylabel('Mediator Level', 'FontSize',16)
 ax=gca;
 ax.FontSize=16;
text(.05*Tmax,1.03*M0,{'Homeostatic Overload'},'FontSize',14,'fontweight','bold')
text(.05*Tmax,.9*M0,{'Reactive Homeostasis'},'FontSize',14,'fontweight','bold')
text(.05*Tmax,.03*M0,{'Predictive Homeostasis'},'FontSize',14,'fontweight','bold')
% patch([T',fliplr(T)'],[(predictivescopefun(T')), fliplr(Y(:,1)')], 'g')

 figure2 = figure;
  ylim([0,1.05])
  hold on
 ax=gca;
 ax.FontSize=16;
box on
  ylabel('Stress Schedule','FontSize',16)
 xlabel('Time','FontSize',16)
 plot(T, allstress(T),'k-', 'LineWidth',1.5)
 
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
