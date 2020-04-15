clear
clc
plot_control
  grey = [0.4,0.4,0.4] ;
  pink = [1.0,0.4,0.6] ;
  purple = [0.5,0,0.5] ;

    fnames = dir('nd*.txt') ;

        fig_name = 'dot_A_D_2nd';
        fig_dum = figure(2);
      set(fig_dum, 'name', fig_name,'numbertitle','on');
      set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
      set(fig_dum,'paperpositionmode','auto');
    

for kk=1:size(fnames)

   id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1};
    area = data{1,3} .* 1E2;
    area(LENGTH > 210 & LENGTH < 220) = 2 .* area(LENGTH > 230 & LENGTH < 240) - area(LENGTH > 250 & LENGTH < 260) ;

    ratio(kk) = area(LENGTH==215) ./ area(LENGTH==195) ;
%%%% plot %%%%
h=scatter(LENGTH,area);%,'black','LineWidth',3);


if (kk==1)
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','k','LineWidth',0.6)
elseif (kk==2)
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','b','LineWidth',0.6)
elseif (kk==3)
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.6)
elseif (kk==4)
    set(h,'MarkerEdgeColor','r','MarkerFaceColor','g','LineWidth',0.6)
elseif (kk==5)
    set(h,'MarkerEdgeColor','b','MarkerFaceColor','y','LineWidth',0.6)
elseif (kk==6)
    set(h,'MarkerEdgeColor','b','MarkerFaceColor','c','LineWidth',0.6)
elseif (kk==7)
    set(h,'MarkerEdgeColor','b','MarkerFaceColor',grey,'LineWidth',0.6)
elseif (kk==8)
    set(h,'MarkerEdgeColor','g','MarkerFaceColor','m','LineWidth',0.6)
elseif (kk==9)
    set(h,'MarkerEdgeColor','g','MarkerFaceColor',purple,'LineWidth',0.6)
end

hold on


end
    
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('ice particle size (\mum)','fontSize',h_axis+6);
ylabel('ice particle projected area (mm^2)','fontSize',h_axis+6);
%title(fnames(kk).name,'fontsize',h_title+6,'fontweight','bold');
box on
ylim([1E-4 1E1])


set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
  set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);
     legend('-20','-25','-30','-35','-40','-45','-50','-55','-60'...
         ,'fontsize',h_legend-4,'location','northwest');
     set(gca,'fontsize',h_axis+6,'LineWidth',2);
    eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
