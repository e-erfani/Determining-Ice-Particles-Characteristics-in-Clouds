clear
clc
plot_control
grey = [0.4,0.4,0.4] ;
pink = [1.0,0.4,0.6] ;
purple = [0.5,0,0.5] ;
phos = [0,0.7,0.7] ;
violet = [0.2,0.2,0.5] ;
banafsh = [0.5,0,1] ;
blu = [0,0.5,1] ;
grin = [0,1,0.4] ;

%%% read data  %%%
fnames = dir('nd*.txt') ;

fig_name = 'dot_fit_warm_m_D_SCPP_new';
fig_dum = figure(1);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');
    
%% read variables; change units %%%
for kk=1:size(fnames)

    id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1}.*1E-4;
    area = data{1,3};
    mass = data{1,4};
    mass(LENGTH>0.01) = 1E-3 .* 0.115 .* (area(LENGTH>0.01) .* 1E2) .^ 1.218 ; % different method for small particles
 
    LL_all(1:length(LENGTH),kk) = LENGTH(:) ;
    mm_all(1:length(LENGTH),kk) = mass(:) ;

    idx2 = find(LENGTH==0.0195) ;
    idx = find(LENGTH==0.0095) ;
    
%%%% scatterplot of observed data %%%%

    if (kk==1)
        h=scatter(LENGTH(1:idx).*1E4,mass(1:idx).*1E3,'v');
        set(h,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',0.6)
    elseif (kk==2)
        h=scatter(LENGTH(1:idx).*1E4,mass(1:idx).*1E3,'s');
        set(h,'MarkerEdgeColor',blu,'MarkerFaceColor',blu,'LineWidth',0.6)
    elseif (kk==3)
        h=scatter(LENGTH(1:idx).*1E4,mass(1:idx).*1E3,'d');
        set(h,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.6)
    elseif (kk==4)
        h=scatter(LENGTH(1:idx).*1E4,mass(1:idx).*1E3,'^');
        set(h,'MarkerEdgeColor',grin,'MarkerFaceColor',grin,'LineWidth',0.6)
    elseif (kk==5)
        h=scatter(LENGTH.*1E4,mass.*1E3);
        set(h,'MarkerEdgeColor','b','MarkerFaceColor','y','LineWidth',0.6)
    elseif (kk==6)
        h=scatter(LENGTH.*1E4,mass.*1E3);
        set(h,'MarkerEdgeColor','b','MarkerFaceColor','c','LineWidth',0.6)
    elseif (kk==7)
        h=scatter(LENGTH.*1E4,mass.*1E3);
        set(h,'MarkerEdgeColor','b','MarkerFaceColor',grey,'LineWidth',0.6)
    elseif (kk==8)
        h=scatter(LENGTH.*1E4,mass.*1E3);
        set(h,'MarkerEdgeColor','g','MarkerFaceColor','m','LineWidth',0.6)
    elseif (kk==9)
        h=scatter(LENGTH.*1E4,mass.*1E3);
        set(h,'MarkerEdgeColor','g','MarkerFaceColor',purple,'LineWidth',0.6)
    end

hold on


end

%%%%%%%%%% read ground data
cd ./SCPP
fnames = dir('m-D*.txt') ;

for kk=1:size(fnames)

    id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1}.*1E-4;
    mass = data{1,2};
   
    %%%% add ground dots to scatterplot %%%%
    h=scatter(LENGTH(1:end).*1E4,mass(1:end),'o');%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor',banafsh,'MarkerFaceColor',banafsh,'LineWidth',0.6)
    hold on
    
end

%%%%%%%%%%%
cd ..     

%%%% implement nonlinear regression model

LL = nanmean(LL_all,2) ;
mm = nanmean(mm_all,2) ;

L_tot = [LL(1:idx2) ; LENGTH(1:end)] ;
m_tot = [mm(1:idx2) ; mass(1:end)*1E-3] ;

p = polyfit(log(L_tot),log(m_tot),2);
f = polyval(p,log(L_tot));

%%% add the curve fit to plot
crvft = p(1).*log(L_max).^2 + p(2).*log(L_max) + p(3) ;
loglog(L_max.*1E4,2.71828.^(crvft).*1E3,'-k','LineWidth',2)
  
%%% calculate accuracy score
[ffit,gof2] = fit(log(L_tot),log(m_tot),'poly2') ;
%            sse: 0.8555
%        rsquare: 0.9970
%            dfe: 30
%     adjrsquare: 0.9968
%           rmse: 0.1689
          
equation=sprintf('synoptic ice clouds');
cor = '-40\circC < T \leq -20\circC' ;
str = {equation , cor};
rrr=annotation('textbox', [.15 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',19)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Ice Particle Mass (mg)','fontSize',h_axis+6);
box on
ylim([1E-6 1E0])
xlim([1E1 1E4])

hold on
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_sphere(m_sphere>0.0005) = NaN ;
loglog(L_max.*1E4,m_sphere.*1E3,'color',grey,'LineWidth',2)

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);
legend('CPI,   -25\circC < T \leq -20\circC   ' ,'CPI,   -30\circC < T \leq -25\circC   ',...
     'CPI,   -35\circC < T \leq -30\circC   ','CPI,   -40\circC < T \leq -35\circC   ',...
     'SCPP, -40\circC < T \leq -20\circC  ','CPI & 2D-S curve fit','CPI & SCPP curve fit','Ice Spheres','-50','-55','-60'...
     ,'fontsize',h_legend-4,'location','southeast');
set(gca,'fontsize',h_axis+3,'LineWidth',2);
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
