clear
clc
plot_control
grey = [0.4,0.4,0.4] ;
pink = [1.0,0.4,0.6] ;
purple = [0.5,0,0.5] ;


fig_name = 'warm_ratio_m_A_2DS';
fig_dum = figure(1);
set(fig_dum, 'name', fig_name,'numbertitle','on');
set(fig_dum,'units','inches','position',[0.3,0.3,8.8,8.8]);
set(fig_dum,'paperpositionmode','auto');

%%%%%%%
cd ../2DS  
  
fnames = dir('nd*.txt') ;    

for kk=1:size(fnames)

    id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1}.*1E-4;
    area = data{1,3};
    mass = data{1,4};
 
    if (kk == 1) 
        LL = LENGTH ;
        AA = area ;
        mm = mass ;
    else
        LL = [LL ; LENGTH] ;
        AA = [AA ; area] ;
        mm = [mm ; mass] ;
    end
   
end


L_max = 15*1E-4:1E-4:max(LL) ;

p = polyfit(log(LL(isnan(AA)==0)),log(AA(isnan(AA)==0)),2);
f = polyval(p,log(LL));

crvft = p(1).*(log(L_max)).^2 + p(2).*log(L_max) + p(3) ;
A_crv = 2.71828 .^ crvft ;

p2 = polyfit(log(LL(isnan(mm)==0)),log(mm(isnan(mm)==0)),2);
f2 = polyval(p2,log(LL));

crvft2 = p2(1).*(log(L_max)).^2 + p2(2).*log(L_max) + p2(3) ;
m_crv = 2.71828 .^ crvft2 ;

m_to_A = m_crv ./ A_crv ;  
  

A_sphere = pi .* L_max .^ 2 ./ 4 ;
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_to_A_sp = m_sphere ./ A_sphere ;
L_max = L_max ( L_max >=0.0015) ;
h2 = semilogx(L_max.*1E4,m_to_A ./ m_to_A_sp,'-k','LineWidth',3) ;
hold on

clear m_to_A m_crv A_crv crvft2 crvft AA mm p p2 f f2 area mass LENGTH LL fnames  
  
  
%%%%%%%%
cd ../Synoptic_100_200_mass_baker_lawson
fnames = dir('nd*.txt') ;
    

for kk=1:size(fnames)

    id = fopen(fnames(kk).name);
    data = textscan(id,'%f %f %f %f');    
    LENGTH = data{1,1}.*1E-4;
    area = data{1,3};
    mass = data{1,4};
    mass(LENGTH>0.01) = 1E-3 .* 0.115 .* (area(LENGTH>0.01) .* 1E2) .^ 1.218 ;  

    if (kk == 1) 
        LL = LENGTH ;
        AA = area ;
        mm = mass ;
    else
        LL = [LL ; LENGTH] ;
        AA = [AA ; area] ;
        mm = [mm ; mass] ;
    end
   
end

L_max = 15*1E-4:1E-4:max(LL) ;

p = polyfit(log(LL(isnan(AA)==0)),log(AA(isnan(AA)==0)),2);
f = polyval(p,log(LL));

crvft = p(1).*(log(L_max)).^2 + p(2).*log(L_max) + p(3) ;

A_crv = 2.71828 .^ crvft ;

p2 = polyfit(log(LL(isnan(mm)==0)),log(mm(isnan(mm)==0)),2);
f2 = polyval(p2,log(LL));

crvft2 = p2(1).*(log(L_max)).^2 + p2(2).*log(L_max) + p2(3) ;
m_crv = 2.71828 .^ crvft2 ;

m_to_A = m_crv ./ A_crv ;

[ffit,gof2] = fit(log(LL),log(AA),'poly2') ;
cor=sprintf('R^2=%g',gof2.rsquare);
equation=sprintf('ln(y)=%gln(x)^2+(%g)ln(x)+(%g)',p(1),p(2),p(3));
equation = 'SPARTICUS synoptic ice clouds' ;
cor = '-40\circC < T \leq -20\circC' ;
str = {equation , cor};
rrr=annotation('textbox', [.35 .8, .1, .1], 'String', str);
set(rrr,'Fontsize',20)


set(gca, 'XScale', 'log')
xlabel('Ice Particle Size (\mum)','fontSize',h_axis+6);
ylabel('(m / A)_i_c_e / (m / A)_s_p_h_e_r_e','fontSize',h_axis+6);
box on

A_sphere = pi .* L_max .^ 2 ./ 4 ;
m_sphere = 0.917 .* pi .* L_max .^ 3 ./ 6 ;
m_to_A_sp = m_sphere ./ A_sphere ;

h2 = semilogx(L_max.*1E4,m_to_A ./ m_to_A_sp,'-b','LineWidth',3) ;

ylim([0.0 1.201])
xlim([1E1 1E4])

set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'Fontsize',25,'linewidth',1.5)
set(gca,'XMinorTick','on','YMinorTick','on','fontsize',h_tick+4);

legend('2DS','CPI & 2DS','fontsize',h_legend-4,'location','southwest');  
set(gca,'fontsize',h_axis+6,'LineWidth',2);
eval(['print -r600 -djpeg ', fig_name,'.jpg']);       
