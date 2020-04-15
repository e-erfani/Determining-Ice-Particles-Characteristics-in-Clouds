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

%% read initial file
fnames = dir('nd*.txt') ;

fig_name = 'dot_fit_warm_A_D_gap';
fig_dum = figure(2);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

%%% read variables
for kk=1:size(fnames)

    id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1}.*1E-4;
    area = data{1,3};
   
    if (kk == 1) 
        LL = LENGTH ;
        mm = area ;
    else
        LL = [LL ; LENGTH] ;
        mm = [mm ; area] ;
    end
    
%%%% scatterplot %%%%

    if (kk==1)
        h=scatter(LENGTH.*1E4,area.*1E2,'o');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',0.6)
    elseif (kk==2)
        h=scatter(LENGTH.*1E4,area.*1E2,'s');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor',blu,'MarkerFaceColor',blu,'LineWidth',0.6)
    elseif (kk==3)
        h=scatter(LENGTH.*1E4,area.*1E2,'d');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.6)
    elseif (kk==4)
        h=scatter(LENGTH.*1E4,area.*1E2,'^');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor',grin,'MarkerFaceColor',grin,'LineWidth',0.6)
    elseif (kk==5)
        h=scatter(LENGTH.*1E4,area.*1E2,'o');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','b','MarkerFaceColor','k','LineWidth',0.6)
    elseif (kk==6)
        h=scatter(LENGTH.*1E4,area.*1E2,'s');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',0.6)
    elseif (kk==7)
        h=scatter(LENGTH.*1E4,area.*1E2,'d');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','b','MarkerFaceColor','r','LineWidth',0.6)
    elseif (kk==8)
        h=scatter(LENGTH.*1E4,area.*1E2,'o');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','g','MarkerFaceColor','k','LineWidth',0.6)
    elseif (kk==9)
        h=scatter(LENGTH.*1E4,area.*1E2,'s');%,'black','LineWidth',3);
        set(h,'MarkerEdgeColor','g','MarkerFaceColor','b','LineWidth',0.6)
    end

    hold on

end


%%% Implement nonlinear regression model

L_max = 25*1E-4:1E-4:3000*1E-4 ; % ; max(LL) ;

p = polyfit(log(LL(isnan(mm)==0)),log(mm(isnan(mm)==0)),2);
f = polyval(p,log(LL));

%%% Add model curve fit to the plot
crvft = p(1).*log(L_max).^2 + p(2).*log(L_max) + p(3) ;
loglog(L_max.*1E4,2.71828.^(crvft).*1E2,'-k','LineWidth',2)
  
%%% Calculate score

[ffit,gof2] = fit(log(LL(isnan(mm)==0)),log(mm(isnan(mm)==0)),'poly2') ;
%            sse: 1.8274
%        rsquare: 0.9978
%            dfe: 206
%     adjrsquare: 0.9978
%           rmse: 0.0942
equation=sprintf('SPARTICUS synoptic ice clouds');
cor = '-40\circC < T \leq -20\circC' ;
str = {equation , cor};
rrr=annotation('textbox', [.15 .81, .1, .1], 'String', str);
set(rrr,'Fontsize',18)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('ice particle projected area (mm^2)','fontSize',h_axis+6);
box on
ylim([1E-4 2E0])
xlim([1E1 1E4])

hold on
m_sphere = pi .* L_max .^ 2 ./ 4 ;
m_sphere(m_sphere>0.01) = NaN ;
loglog(L_max.*1E4,m_sphere.*1E2,'color',grey,'LineWidth',2)

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);
legend('-25\circC < T \leq -20\circC   ' ,'-30\circC < T \leq -25\circC   ',...
     '-35\circC < T \leq -30\circC   ','-40\circC < T \leq -35\circC   ',...
     'curve fit','Ice Spheres','-50','-55','-60'...
     ,'fontsize',h_legend-4,'location','southeast');
set(gca,'fontsize',h_axis+6,'LineWidth',2);
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
