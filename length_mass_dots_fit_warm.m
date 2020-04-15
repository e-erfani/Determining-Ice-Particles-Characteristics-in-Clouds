clear
clc
plot_control
  grey = [0.4,0.4,0.4] ;
  pink = [1.0,0.4,0.6] ;
  purple = [0.5,0,0.5] ;

    fnames = dir('nd*.txt') ;

        fig_name = 'dot_fit_warm_m_D';
        fig_dum = figure(1);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');
    

for kk=1:4%size(fnames)

   id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1}.*1E-4;
    mass = data{1,4};
    area = data{1,3};
    area(LENGTH > 0.021 & LENGTH < 0.022) = 2 .* area(LENGTH > 0.023 & LENGTH < 0.024) - area(LENGTH > 0.025 & LENGTH < 0.026) ;
   % mass(LENGTH > 0.01 & LENGTH < 0.02) = NaN ;
    mass(LENGTH>0.01) = 1E-3 .* 0.115 .* (area(LENGTH>0.01) .* 1E2) .^ 1.218 ;
 
    if (kk == 1) 
        LL = LENGTH ;
        mm = mass ;
    else
        LL = [LL ; LENGTH] ;
        mm = [mm ; mass] ;
    end
    

%%%% plot %%%%


if (kk==1)
    h=scatter(LENGTH.*1E4,mass.*1E3,'o');%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','k','LineWidth',0.6)
elseif (kk==2)
    h=scatter(LENGTH.*1E4,mass.*1E3,'s');%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','b','LineWidth',0.6)
elseif (kk==3)
    h=scatter(LENGTH.*1E4,mass.*1E3,'d');%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.6)
elseif (kk==4)
    h=scatter(LENGTH.*1E4,mass.*1E3,'^');%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','g','LineWidth',0.6)
elseif (kk==5)
    h=scatter(LENGTH.*1E4,mass.*1E3);%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','b','MarkerFaceColor','y','LineWidth',0.6)
elseif (kk==6)
    h=scatter(LENGTH.*1E4,mass.*1E3);%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','b','MarkerFaceColor','c','LineWidth',0.6)
elseif (kk==7)
    h=scatter(LENGTH.*1E4,mass.*1E3);%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','b','MarkerFaceColor',grey,'LineWidth',0.6)
elseif (kk==8)
    h=scatter(LENGTH.*1E4,mass.*1E3);%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','g','MarkerFaceColor','m','LineWidth',0.6)
elseif (kk==9)
    h=scatter(LENGTH.*1E4,mass.*1E3);%,'black','LineWidth',3);
    set(h,'MarkerEdgeColor','g','MarkerFaceColor',purple,'LineWidth',0.6)
end

hold on


end


L_max = 25*1E-4:1E-4:6000*1E-4 ; %max(LL) ;

p = polyfit(log(LL(isnan(mm)==0)),log(mm(isnan(mm)==0)),2);
f = polyval(p,log(LL));

crvft = p(1).*log(L_max).^2 + p(2).*log(L_max) + p(3) ;
loglog(L_max.*1E4,2.71828.^(crvft).*1E3,'-k','LineWidth',2)
  

%  hold on
% crvft2 = -0.120320.*log(L_max).^2 + 1.35772.*log(L_max) - 6.59802 ;
% loglog(L_max.*1E4,2.71828.^(crvft2).*1E3,'-b','LineWidth',2)

[ffit,gof2] = fit(log(LL(isnan(mm)==0)),log(mm(isnan(mm)==0)),'poly2') ;

cor=sprintf('R^2=%g',gof2.rsquare);
equation=sprintf('ln(y)=%gln(x)^2+(%g)ln(x)+(%g)',p(1),p(2),p(3));
str = {equation , cor};
rrr=annotation('textbox', [.15 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',15)


set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('Ice Particle Mass (mg)','fontSize',h_axis+6);
%title(fnames(kk).name,'fontsize',h_title+6,'fontweight','bold');
box on
ylim([1E-6 1E0])
xlim([1E1 1E4])

hold on
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_sphere(m_sphere>0.0001) = NaN ;
loglog(L_max.*1E4,m_sphere.*1E3,'--k','LineWidth',2)

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);
     legend('-20','-25','-30','-35','curve fit','Ice Spheres','-50','-55','-60'...
         ,'fontsize',h_legend-4,'location','southeast');
     set(gca,'fontsize',h_axis+6,'LineWidth',2);
    eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
