function vis_s2(params,t_g, t_bg, t_m, string)

[temp,n_g]= size(t_g);
[temp,n_bg]= size(t_bg);
[temp,n_m]= size(t_m);

dim= n_g*n_bg*n_m;

% (1,dim) dimensional arrays populated with thickness for plotting
t_g1= zeros(1,dim);
t_bg1= zeros(1,dim);
t_m1= zeros(1,dim);

% Initialize arrays for parameters of interest
max_dips= zeros(1,dim);
fwhms= zeros(1,dim);
sens= zeros(1,dim);
FOM= zeros(1,dim);

t_g1= repmat(t_g,1,(dim/n_g));

for i= 1:(n_bg)
   temp(((i-1)*n_g +1) : (i*n_g))= repmat(t_bg(i), 1,n_g);
end
t_bg1= repmat(temp,1,n_m);

for i= 1:(n_m)
   t_m1(((i-1)*n_g*n_bg +1) : (i*n_g*n_bg))= repmat(t_m(i),1, n_g*n_bg); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:dim
   max_dips(i)= 1- mean(params(:, 1, i));
   sens(i)= (params(3,2,i)-params(1,2,i))/(0.02);
   fwhms(i)= mean(params(:, 3, i));
end

%%%%%%%%%%%%%%%%%%%%%%%% Calculate FOM

w1= 2; w2=0.5; w3= 0.5;
FOM= w1.*normalize(max_dips,'range') ...
        + w2.*normalize(sens, 'range')...
        + w3.*normalize(1./fwhms, 'range');

%%%%%%%%%%%%%%%%%%%%%%%% Get some params :P

% figure;
% plot(FOM);

% [maxx,index]= max(max_dips);
% index=index-11
% t_g1(index)
% t_bg1(index)
% t_m1(index)
% max_dips(index)
% sens(index)
% fwhms(index)
% FOM(index)

%%%%%%%%%%%%%%%%%%%%%%%% 4D colour plot (params vs FOM)
if(strcmp(string,'plotA') | strcmp(string,'plot4'))
    figure;
    scatter3(t_g1.*10^9,t_bg1.*10^9,t_m1.*(10^9), 100, FOM,'MarkerFaceColor', 'Flat');
    colorbar();
    caxis([min(FOM), max(FOM)]);
    xlabel('Tg (nm)');
    ylabel('Tbuff-Tg (nm)');
    zlabel('Tm (nm)');
end
    
%%%%%%%%%%%%%%%%%%%%%%%% Optionally plot graphs

if(strcmp(string,'plotA') | strcmp(string,'plotS'))
    figure;
    plot(max_dips);
    xlabel('change in counter');
    ylabel('max dips');

    figure;
    plot(sens);
    xlabel('change in counter');
    ylabel('Sensitivity');

    figure;
    plot(fwhms);
    xlabel('change in counter');
    ylabel('Full width at half maxima');
end

end