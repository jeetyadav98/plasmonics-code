function vis_s1(params,t_al, string)


sz= size(t_al);
max_dips= zeros(size(t_al));
fwhms= zeros(size(t_al));
sens= zeros(size(t_al));
FOM= zeros(size(t_al));

for i=1:sz(2)
   max_dips(i)= 1-params(2, 1, i);
   fwhms(i)= params(2, 3, i);
   sens(i)= (params(3,2,i)-params(1,2,i))/(0.02);
end

%%%%%%%%%%%%%%%%%%%%%%%% Calculate FOM

w1= 1; w2=0.5; w3= 0.5;
FOM= w1.*normalize(max_dips,'range') ...
        + w2.*normalize(sens, 'range')...
        + w3.*normalize(1./fwhms, 'range');
    
%%%%%%%%%%%%%%%%%%%%%%%% Optionally plot graphs

if(strcmp(string,'plotA') | strcmp(string,'plotS'))
    figure;
    plot(t_al.*10^9,max_dips);
    xlabel('thickness of aluminium');
    ylabel('max dips');

    figure;
    plot(t_al.*10^9,sens);
    xlabel('thickness of aluminium');
    ylabel('Sensitivity');

    figure;
    plot(t_al.*10^9,fwhms);
    xlabel('thickness of aluminium');
    ylabel('Full width at half maxima');
end

if(strcmp(string,'plotA') |strcmp(string,'plotC'))
    figure;
    hold on;
    plot(t_al.*10^9,normalize(max_dips,'range'));
    plot(t_al.*10^9,normalize(sens,'range'));
    plot(t_al.*10^9,normalize(fwhms,'range'));
    legend('Max dips', 'Sensitivity', 'FWHM');
    hold off;
end

if(strcmp(string,'plotA') |strcmp(string,'plotFOM'))
    figure;
    plot(t_al.*10^9, FOM);
end

end