%%%%%%%%%% for tm binary grating glass metal  metalgrating dielectric
%%%%%%%%%% wavelength modulation
clc; clear all;

wbar= waitbar(0, 'Progress Bar');

c=3e8;
e=1.62e-19;
h=6.62e-34;
% omegal=2*3.14*650.07e12;
% gama=15.92*2*3.14e12;
% gamal=104.86*2*3.14e12;
% omegap=2113.6*2*3.14e12;
% deltaepsilon=1.09;

count=1;


for kkk=1:1:1
  
    n_omega= 2827;
  for iii=1:1:n_omega
    
    omega=9.42e14+(iii-1)*0.005e14;
    lambda0=(2*3.14*c)./omega;
    lambda(iii)=lambda0;
    lamdac=2.4511e-5;
    lamdap=1.0657e-7;

    
    epsilonreal=1-((lambda0.^2.*(lamdac)^2)./(lamdap^2.*(lambda0.^2+lamdac^2)));
    epsilonim=((lambda0.^3.*(lamdac))./(lamdap^2.*(lambda0.^2+lamdac^2)));
    
    n2=sqrt((sqrt(epsilonreal.^2+epsilonim.^2)+epsilonreal)./2);
    k2=sqrt((sqrt(epsilonreal.^2+epsilonim.^2)-epsilonreal)./2);

    n22(iii)=n2;
    k22(iii)=k2;

    Numords=101;    %%%%%%%%% number of diffractive orders maintained
    nc=1.64;        %%%%%%%%% region 1 cover refractive index
    ns=1;           %%%%%%%%% region 3 substrate refractive index
    Ngrat=3;        %%%%%%%%% number of grating slices
    period=1100e-9; %%%%%%%%% grating period in microns

    nd=1.33;
    nm=n2-1i*k2;

    depth=[15e-9,40e-9,1000e-9];  %%%% Height for each grating
    j=sqrt(-1);

    nr=[nm,nm,nd];                %%%%%%%%%% Ridge refractive index for each grating
    ng=[nm,nd,nd];                %%%%%%%%%% index for ridge each grating
    Filfac=[.5 .5 .5 ];           %%%%%%%%%% fill factor for ridges
    Disp=[0 0 0 ];                %%%%%%%%%% ridge displacement in a frac                                                                                                                                                                                                                                              tion of period

    theta0=0;                     %%%%%%%%%% angle of incidence
    phi0=0;                       %%%%%%%%%% azimuthal angle of incidence
    deg=pi/180; 

    Nmax=(Numords-1)/2;           %%%%%%%%%% highest order number retained
    I=(-Nmax:Nmax)';              %%%%%%%%%% I is the order index
    p=(Numords+1)/2;              %%%%%%%%%% index of zeroth order

    theta0=theta0*deg;
    phi0=phi0*deg;                %%%%%%%%%% converting in radians

    epsc=nc^2;                    %%%%%%%%%% relative permittivity
    epss=ns^2;
    k0=2*pi/lambda0;              %%%%%%%%%% free space vector
    K=2*pi/period;                %%%%%%%%%% grating vector %%%%%%%%%%%%
    kc=k0*nc;

    kx=kc*sin(theta0)*cos(phi0)-I*K;  %%%%%%%% region1 wave vector components
%     kspp(iii)=k0*sqrt((epsilonreal*nd)/(epsilonreal+nd));   % Surface plasmon wave vector
%     kx1(iii)=kx(p);
%     kx2(iii)=kx(p-1);
%     kx3(iii)=kx(p-2);
%     kx4(iii)=kx(p+1);
%     kx5(iii)=kx(p+2);

    ky=kc*sin(theta0)*sin(phi0)*ones(size(kx)); %y direction
    
    kzc=sqrt(kc^2-kx.^2-ky.^2);
    bad_indices=find((real(kzc)-imag(kzc))<0);  
    kzc(bad_indices)=-kzc(bad_indices);              
    ks=k0*ns;
    kzs=sqrt(ks^2-kx.^2-ky.^2);
    bad_indices=find((real(kzs)-imag(kzs))<0);   %%%%%%%%%%%region3 wavevector
    kzs(bad_indices)=-kzs(bad_indices); 

    %%%%%%%%%%define some  matrices and vectors
    Zv=zeros(Numords,1);
    Dv=Zv;
    Dv(p)=1;
    Zv2=[Zv;Zv];
    Eye=eye(Numords);  %%%%%%% identity matrix
    Kx=diag(kx/k0);
    Kxsq=Kx.^2;
    Kzc=diag(kzc/k0);
    Kzcsq= Kzc.^2;
    Kzs=diag(kzs/k0);
    Kzssq= Kzs.^2;
    M=Numords-1;
    temp1=Eye/ns;
    fmat=Eye;
    gmat=j*Kzs/ns^2;
       
    % tc=-atan2(kx,kzc);
    % ts=-atan2(kx,kzs);
    % tc=tc.*(180/pi);
    % ts=ts.*(180/pi);

    for ii=Ngrat:-1:1
      epsg=ng(ii).^2;                             %%%%%%%%% groove permittivity
      epsr=nr(ii).^2;                             %%%%%%%%% ridge permittivity

      epsG=(1-Filfac(ii))*epsg+Filfac(ii)*epsr;   %%%%%%%%% average grating
      iepsG=(1-Filfac(ii))/epsg+Filfac(ii)/epsr; 
      Sinc=sin(pi*Filfac(ii)*[1:M])./(pi*[1:M]);

      vm=(epsr-epsg)*fliplr(Sinc);
      v0=epsG;
      vp=(epsr-epsg)*Sinc;
      v=[vm v0 vp].*exp(+j*2*pi*Disp(ii)*[-M:M]);
      ivm=(1/epsr-1/epsg)*fliplr(Sinc);
      iv0=iepsG;
      ivp=(1/epsr-1/epsg)*Sinc;
         
      iv=[ivm iv0 ivp].*exp(+j*2*pi*Disp(ii)*[-M:M]);
      Epsilon=toeplitz(fliplr(v(1:Numords)),v(Numords:2*Numords-1));
      Alpha=toeplitz(fliplr(iv(1:Numords)),iv(Numords:2*Numords-1));
      
      clear Sinc v  vm  v0  vp 

      B=Kx*(Epsilon\Kx)-Eye;  %%%%%%%%%%% cofficient matrix
      [W,V]=eig(Alpha\B);     %%%%%%%%%%% W is the eigen vector and V are the eigen values
      Q=sqrt(V);
      M0=Alpha*W*Q;
      E=expm(-k0*Q*depth(ii));
      v=[W W;M0,-M0]\[fmat;gmat];

      temp2=v(1:Numords,:)\E;
      temp3=E*v(Numords+1:2*Numords,:)*temp2;
      temp1=temp1*temp2;
      fmat=W+W*temp3;
      gmat=M0-M0*temp3;
    end

    gfi=gmat/fmat;
    RHS=-gfi(:,p);
    RHS(p)=RHS(p)+j*kzc(p)/k0/epsc;
    Rs=[gfi+j*Kzc/nc^2]\RHS;
    Ts=(temp1/fmat)*(Rs+Dv)*nc;
     

    IR1=(abs(Rs).^2).*real(kzc./kzc(p));
    IT1=(abs(Ts).^2).*real(kzs./kzc(p));


    e=sum(IT1);
    f=sum(IR1);
    g=1-e-f;
    IT12(count)=e;
    IR12(count)=f; 
    loss(count)=g;
    count=count+1;
    
    waitbar(count/n_omega,wbar, sprintf('Progress: %.2f %%', (count*100/n_omega)));
    
  end

  lambda=lambda.*10^9;
  plot(lambda,IR12);
  hold all

  count=1;
  close(wbar);
end