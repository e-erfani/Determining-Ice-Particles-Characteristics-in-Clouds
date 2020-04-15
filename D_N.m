clear
clc
plot_control
grey = [0.4,0.4,0.4] ;
pink = [1.0,0.4,0.6] ;
purple = [0.7,0,0.7] ;
blu = [0,0,1] ;

%   Establish constants:
rhoi=0.917;  % density of bulk ice
rhoi2=0.50;  % ice density assumed in CAM5
beta3=3.0 ;  % m-D power assumed for ice in CAM5
alpha3=(pi/6.)*rhoi2 ;  % m-D prefactor assumed for ice in CAM5
delta3=2.0;  % A-D power assumed for ice in CAM5
gamma3=pi/4.;           % A-D prefactor assumed for ice in CAM5
nu=0.0 ;     % PSD dispersion parameter

% constants and coefficients for Heymsfield and Westbrook fallspeed      
g = 981. ;
ra = 2.867e6 ;  % gas constant for air 
p0=500. ;  %   pressure at 500 hPa  
p=1000.*p0 ; % Convert to dynes/cm2:
tk=273.-20. ;  % Temperature assumed at 500 hPa  unit: K
delH=8.0 ; %  Coeffs. for governing eqn. via Bohm (1989) for ice particles
c0H=0.35 ; %  Coeffs. for governing eqn. via Bohm (1989) for ice particles

      
D = 0.001:0.0001:0.1 ;      

%  m-D expression of Cotton et al. (2012, QJRMS):
alpha4(D<=0.007) = (pi./6).*0.700 ;   % for cgs units (g/cm3)
beta4(D<=0.007) = 3.00 ;
alpha4(D>0.007) = 0.00257 ; % for cgs units
beta4(D>0.007) = 2.00 ;
m_cot = 1.e3.*alpha4.*D.^beta4 ;  % milligram units
    
%   Select value for IWC in g/cm3:
      iwc=5.e-8 ;  % IWC = 10 mg/m3
      
%   Define constants for m-D 2nd order polynomial fit for ln(m) vs. ln(D)
%   for synoptic ice clouds, -40 < T < -20 C.  Curve fit is based on cgs units.
a0=-6.72924;
a1=1.17421;
a2=-0.15980;

% for ln(A) vs. ln(D)
ab0= -2.46356 ; %-2.12631;
ab1= 1.25892 ; %1.48745;
ab2= -0.07845; %-0.0376357;

%   Calculate beta and alpha for power law m = alpha*D**beta based on above polynomial fit
%    and an ice particle maximum dimension D of 500 mirons. Also do this for the projected area
%     power law A = gamma*D**delta:
D0=0.0500 ;  % reference dimension; cm units
beta = a1 + 2.*a2.*log(D0);  % synoptic cirrus; -40 < T < -10C
delta = ab1 + 2.*ab2.*log(D0) ;  % synoptic cirrus; -40 < T < -10C
%   Now calc. the prefactors:
mlog=a0 + a1.*log(D0) +a2.*(log(D0)).^2;
m500=exp(mlog);
alpha=m500/(D0.^beta);
Alog=ab0 + ab1.*log(D0) + ab2.*(log(D0)).^2;
A500=exp(Alog);
gama=A500/(D0.^delta);
     

slp=(nu+1.) ./ D ;
d_N=(nu+0.67) ./ slp  ;       % median number conc. dimension
d_A=(delta+nu+0.67) ./ slp  ; % median area dimension (approx.)
d_m=(beta+nu+0.67) ./ slp  ;  % median mass dimension (approx.)

%     Calc. the ice particle number concentration using methods 1 & 2:
n1 = gamma(nu+1.).*iwc.*slp.^beta ./ (alpha.*gamma(beta+nu+1.)) ;  % constant alpha and beta... 
n3 = gamma(nu+1.).*iwc.*slp.^beta3 ./ (alpha3.*gamma(beta3+nu+1.))  ; % CAM5 method
n4 = gamma(nu+1.).*iwc.*slp.^beta4 ./ (alpha4.*gamma(beta4+nu+1.))  ; % Cotton method
%       New method for N from paper:
beta2 = a1 + 2.*a2.*log(d_N) ;  % synoptic cirrus; -40 < T < -10C
mlog = a0 + a1.*log(d_N) +a2.*(log(d_N)).^2 ;  
mtmp = exp(mlog) ;
alpha2 = mtmp ./ (d_N.^beta2) ;
n2 = gamma(nu+1.).*iwc.*slp.^beta2 ./ (alpha2.*gamma(beta2+nu+1.))  ; % new method

d_N2 = d_N ;

for j=1:length(D)
    d_N(j,1) = d_N2(j) ;
end

for i = 2:5   

    %      - Mass for Method 2
    beta2 = a1 + 2.*a2.*log(d_N(:,i-1)) ; % synoptic cirrus; -40 < T < -10C
    mlog = a0 + a1.*log(d_N(:,i-1)) + a2.*log(d_N(:,i-1)).^2 ;  % ditto
    mtmp = exp(mlog) ;
    alpha2 = mtmp ./ (d_N(:,i-1).^beta2) ;

    n2_itrt = gamma(nu+1.).*iwc.*slp'.^beta2 ./ (alpha2.*gamma(beta2+nu+1.))  ; % new method

    d_N(:,i) = (nu+0.67) ./ slp(:) ;    % median mass dimension (approx.)      
      
end      
      
%%%%%%%%%%%%%      
fig_name = 'D_N';
fig_dum = figure(1);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

h1 = loglog(D(find(D==0.004):end).*1E4,n1(find(D==0.004):end).*1E3,'-k','LineWidth',2.5) ;
hold on
h2 = loglog(D(find(D==0.004):end).*1E4,n2(find(D==0.004):end).*1E3,'b','LineWidth',2.5) ;
hold on
h3 = loglog(D(find(D==0.004):end).*1E4,n3(find(D==0.004):end).*1E3,'-r','LineWidth',2.5) ;
hold on
h4 = loglog(D(find(D==0.004):end).*1E4,n4(find(D==0.004):end).*1E3,'color',purple,'LineWidth',2.5) ;
hold on

IWC_n = 'IWC = 50 mg m^-^3' ;
equation = 'synoptic ice clouds' ;
cor = '-40\circC < T \leq -20\circC' ;
str = {equation  , IWC_n , cor};
rrr=annotation('textbox', [.143 .34, .1, .1], 'String', str);
set(rrr,'Fontsize',20)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Number Concentration (liter^-^1)','fontSize',h_axis+6);
box on
%ylim([1E-6 1E0])
%xlim([1E1 1E4])

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+7);

legend('N using constant \alpha & \beta','N using variable \alpha & \beta',...
 'N using CAM5 \alpha & \beta','N using Cotton et al. (2012) \alpha & \beta',h_legend-8,'location','southwest');  
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       

n_error = 100 .* 2.* (n1 - n2) ./ (n1 + n2) ;
      
%%%%%%%%%%%%%      
fig_name = 'D_N_iterate';
fig_dum = figure(10);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

h1 = loglog(D(find(D==0.004):end).*1E4,n1(find(D==0.004):end).*1E3,'-k','LineWidth',2.5) ;
hold on
h2 = loglog(D(find(D==0.004):end).*1E4,n2_itrt(find(D==0.004):end).*1E3,'b','LineWidth',2.5) ;
hold on
h3 = loglog(D(find(D==0.004):end).*1E4,n3(find(D==0.004):end).*1E3,'-r','LineWidth',2.5) ;
hold on
h4 = loglog(D(find(D==0.004):end).*1E4,n4(find(D==0.004):end).*1E3,'color',purple,'LineWidth',2.5) ;
hold on

IWC_n = 'IWC = 50 mg m^-^3' ;
equation = 'synoptic ice clouds' ;
cor = '-40\circC < T \leq -20\circC' ;
str = {equation  , IWC_n , cor};
rrr=annotation('textbox', [.143 .34, .1, .1], 'String', str);
set(rrr,'Fontsize',20)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Number Concentration (liter^-^1)','fontSize',h_axis+6);
box on
%ylim([1E-6 1E0])
%xlim([1E1 1E4])

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+7);

legend('N using constant \alpha & \beta','N using variable \alpha & \beta',...
 'N using CAM5 \alpha & \beta','N using Cotton et al. (2012) \alpha & \beta',h_legend-8,'location','southwest');  
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%    Calc. PSD effective diameter De using methods 1 and 2:
%      - Area for Method 2
delta2 = ab1 + 2.*ab2.*log(d_A) ;   % synoptic cirrus; -40 < T < -10C
Alog = ab0 + ab1.*log(d_A) + ab2.*log(d_A).^2 ;  
Atmp = exp(Alog) ;
gamma2 = Atmp ./ (d_A .^ delta2) ;
%      - Mass for Method 2
beta2 = a1 + 2.*a2.*log(d_m) ; % synoptic cirrus; -40 < T < -10C
mlog = a0 + a1.*log(d_m) + a2.*log(d_m).^2 ;  % ditto
mtmp = exp(mlog) ;
alpha2 = mtmp ./ (d_m.^beta2) ;

d_m2 = d_m ;
d_A2 = d_A ;

for j=1:length(D)
    d_m(j,1) = d_m2(j) ;
    d_A(j,1) = d_A2(j) ;
end

for i = 2:5   

    %      - Area for Method 2    
    delta2 = ab1 + 2.*ab2.*log(d_A(:,i-1)) ;   % synoptic cirrus; -40 < T < -10C
    Alog = ab0 + ab1.*log(d_A(:,i-1)) + ab2.*log(d_A(:,i-1)).^2 ;  
    Atmp = exp(Alog) ;
    gamma2 = Atmp ./ (d_A(:,i-1) .^ delta2) ;
    %      - Mass for Method 2
    beta2 = a1 + 2.*a2.*log(d_m(:,i-1)) ; % synoptic cirrus; -40 < T < -10C
    mlog = a0 + a1.*log(d_m(:,i-1)) + a2.*log(d_m(:,i-1)).^2 ;  % ditto
    mtmp = exp(mlog) ;
    alpha2 = mtmp ./ (d_m(:,i-1).^beta2) ;

    %      Calc. of De using variable alpha, beta, gamma & delta:
    De2 = (3./2).*alpha2.*gamma(beta2+nu+1.).*(slp'.^(delta2-beta2)) ./ (rhoi.*gamma2.*gamma(delta2+nu+1.)) ;

    d_m(:,i) = (beta2(:)+nu+0.67) ./ slp(:) ;    % median mass dimension (approx.)      
    d_A(:,i) = (delta2(:)+nu+0.67) ./ slp(:) ;    % median mass dimension (approx.)      
      
end

%      Calc. of standard De:
De1 = (3./2).*alpha.*gamma(beta+nu+1.).*(slp.^(delta-beta)) ./ (rhoi.*gama.*gamma(delta+nu+1.)) ;
%      Calc. of De using CAM5 alpha, beta, gamma & delta:
De3 = (3./2).*alpha3.*gamma(beta3+nu+1.).*(slp.^(delta3-beta3)) ./ (rhoi.*gamma3.*gamma(delta3+nu+1.)) ;


%%%%%%%%
fig_name = 'D_De';
fig_dum = figure(2);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

h1 = loglog(D.*1E4,De1.*1E4,'-k','LineWidth',2.5) ;
hold on
h2 = loglog(D.*1E4,De2.*1E4,'b','LineWidth',2.5) ;
hold on
h3 = loglog(D.*1E4,De3.*1E4,'-r','LineWidth',2.5) ;
hold on

equation = 'synoptic ice clouds' ;
cor = '-40\circC < T \leq -20\circC' ;
str = {equation  , cor};
rrr=annotation('textbox', [.15 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',20)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Effective Diameter (\mum)','fontSize',h_axis+6);
box on
ylim([1E1 1E3])
%xlim([1E1 1E4])

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+7);

legend('D_e using constant \alpha, \beta, \gamma & \delta','D_e using variable \alpha, \beta, \gamma & \delta',...
 'D_e using CAM5 \alpha, \beta, \gamma & \delta','D_e using Cotton et al. (2012) \alpha & \beta',h_legend-8,'location','southeast');  
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       

De_error = 100 .* 2.* (De1' - De2) ./ (De1' + De2) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calc. the mass-weighted ice fall speed using methods 1 & 2.
%     Use schemes of Heymsfield & Westbrook (2010) 

for i = 2:5   
    
      delta2 = ab1 + 2.*ab2.*log(d_m(:,i-1)) ;   % synoptic cirrus; -40 < T < -10C
      Alog = ab0 + ab1.*log(d_m(:,i-1)) + ab2.*log(d_m(:,i-1)).^2 ;  
      Atmp = exp(Alog) ;
      gamma2 = Atmp ./ (d_m(:,i-1) .^ delta2) ;
%      - Mass for Method 2
      beta2 = a1 + 2.*a2.*log(d_m(:,i-1)) ; % synoptic cirrus; -40 < T < -10C
      mlog = a0 + a1.*log(d_m(:,i-1)) + a2.*log(d_m(:,i-1)).^2 ;  % ditto
      mtmp = exp(mlog) ;
      alpha2 = mtmp ./ (d_m(:,i-1).^beta2) ;
    
    
      area_x2 = gamma2.*d_m(:,i-1).^delta2 ;  %   best est. (variable gamma & delta)
      area_s = (pi./4.).*d_m(:,i-1).^2. ;      %   area for ice spheres
      aratio2 = area_x2 ./ area_s  ;     %   best est.

%     Calculate air density and dynamic & kinematic viscosity
      rhoa = p ./ (ra.*tk) ;
      visc = 2.48e-6 .* (tk .^ 0.754) ;  %    dynamic viscosity
      visk = visc ./ rhoa ;            %    kinematic viscosity      
      
%     Calc. the modified Best number
      xxh2 = 8.*mtmp.*g.*rhoa ./ (pi.*visc.^2.*aratio2.^0.5) ;
      
%     Calc. Reynolds number:
      ReH2 = (delH.^2./4).*((1. + (4.*(xxh2.^0.5)) ./ (delH.^2.*c0H.^0.5)).^0.5 - 1.).^2 ;

      vHW2 = ReH2.*visk./d_m(:,i-1) ; %    fallspeed best estimate (variable m-D/A-D terms) 

      d_m(:,i) = (beta2(:)+nu+0.67) ./ slp(:) ;    % median mass dimension (approx.)      
      
end

area_x1 = gama.*d_m(:,1).^delta ;    %   fixed gamma & delta (std. method)
area_x3 = gamma3.*d_m(:,1).^delta3 ;  %   fixed gamma & delta from CAM5
area_s = (pi./4.).*d_m(:,1).^2. ;      %   area for ice spheres

aratio1 = area_x1 ./ area_s ;      %   fixed gamma & delta
aratio3 = area_x3 ./ area_s ;      %   fixed gamma & delta from CAM5
vmass = alpha.*d_m(:,1).^beta  ;      %   fixed alpha & beta (std. method)
vmass3 = alpha3.*d_m(:,1).^beta3 ;    %   fixed alpha & beta from CAM5

xxh1 = 8.*vmass.*g.*rhoa ./ (pi.*visc.^2.*aratio1.^0.5) ;
xxh3 = 8.*vmass3.*g.*rhoa ./ (pi.*visc.^2.*aratio3.^0.5) ;

ReH1 = (delH.^2./4).*((1. + (4.*(xxh1.^0.5)) ./ (delH.^2.*c0H.^0.5)).^0.5 - 1.).^2 ;
ReH3 = (delH.^2./4).*((1. + (4.*(xxh3.^0.5)) ./ (delH.^2.*c0H.^0.5)).^0.5 - 1.).^2 ;

vHW1 = ReH1.*visk./d_m(:,1) ;  %    fallspeed using fixed alpha, beta, gamma, delta
vHW3 = ReH3.*visk./d_m(:,1) ; %    fallspeed using CAM5 alpha, beta, gamma, delta
      
%%%%%%%%
fig_name = 'D_Vm';
fig_dum = figure(3);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

h1 = loglog(D.*1E4,vHW1,'-k','LineWidth',2.5) ;
hold on
h2 = loglog(D.*1E4,vHW2,'b','LineWidth',2.5) ;
hold on
h3 = loglog(D.*1E4,vHW3,'-r','LineWidth',2.5) ;
hold on

equation = 'synoptic ice clouds' ;
cor = '-40\circC < T \leq -20\circC' ;
str = {equation  , cor};
rrr=annotation('textbox', [.15 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',20)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Mean Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Mass-weighted Ice Fallspeed (cm s^-^1)','fontSize',h_axis+6);
box on
%ylim([1E1 1E3])
%xlim([1E1 1E4])

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+7);

legend('V_m using constant \alpha, \beta, \gamma & \delta','V_m using variable \alpha, \beta, \gamma & \delta',...
     'V_m using CAM5 \alpha, \beta, \gamma & \delta','D_e using Cotton et al. (2012) \alpha & \beta',h_legend-8,'location','southeast');  
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       


vHW_error = 100 .* 2.* (vHW1 - vHW2) ./ (vHW1 + vHW2) ;


%%%

fig_name = 'D_Di';
fig_dum = figure(100);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

h1 = loglog(D.*1E4,d_N(:,5).*1E4,'-k','LineWidth',2.5) ;
hold on
h2 = loglog(D.*1E4,D.*1E4,'-g','LineWidth',2.5) ;
hold on
h3 = loglog(D.*1E4,d_A(:,5).*1E4,'-b','LineWidth',2.5) ;
hold on
h4 = loglog(D.*1E4,d_m(:,5).*1E4,'-r','LineWidth',2.5) ;

set(gca, 'XScale', 'log')
xlabel('Mean Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Dimension of Interest (\mum)','fontSize',h_axis+6);
box on
%ylim([1E1 1E3])
%xlim([1E1 1E4])

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+7);

legend('D_N','D_m_e_a_n',...
 'D_A','D_m',h_legend-8,'location','southeast');  
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
