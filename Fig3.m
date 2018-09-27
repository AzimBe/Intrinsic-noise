%   *****************************************************************
%   ***************************** Cooerative   ***********************
clear
clc
figure()
S0=0.135;  % S_0
Sa=1.8;    % S_a
K=1.;

cd fnM
xd = fig3ac ; % SSA results


p= [1.  -(S0+Sa) 1. -S0];
r=roots(p);  % roots of first derivative of potential 
xa=r(3);  % concentration at stable point A
xb=r(2);  % concentration at unstable point B
xc=r(1);  % concentration at stable point C
ua= -(S0+Sa)*xa+Sa*atan(xa)+xa^2/2.; % potential at point A 
ub= -(S0+Sa)*xb+Sa*atan(xb)+xb^2/2.; % potential at point B 
uc= -(S0+Sa)*xc+Sa*atan(xc)+xc^2/2.; % potential at point C 
uza= -2.*Sa*xa/(1+xa^2)^2+1 ; % second derivative of potential at point A 
uzb= -2.*Sa*xb/(1+xb^2)^2+1 ; % second derivative of potential at point B 
uzc= -2.*Sa*xc/(1+xc^2)^2+1 ; % second derivative of potential at point C
tr= [20:20:100]; 
tr2= [20 60 140 200 300]; 
tc2= [1 3 7 10 11];

    % escape time A to C
    for i=1:5
        pt=20.*i ;
        V = pt/xa ;
        tac(i) = xd(i,1) ; % SSA, escape time A to C at T(i), translation compensated
        eac(i) = xd(i,2)  ; % SSA, error bar 
        gp=log(2.); % gamma_p
        gt= 1. ;  % 1. / T
        d=gp*sqrt(W2d(xa,S0,Sa,V,gt,K)*abs(V2d(xb,S0,Sa,V,gt,K))) ;
        ex=exp( Vd(xb,S0,Sa,V,gt,K) - Vd(xa,S0,Sa,V,gt,K) - log(ttd1(xa,V,gt,K))) ;  
        tta(i)= 2*pi*ex/d;  % Theoretical escape time A to C , translation compensated
        
        % constant theta
        kt=xa*2./V  ;           % theta, transC
        d=gp*sqrt(uza*abs(uzb)) ;  
        ex=exp((ub-ua)/kt) ; 
        ttac(i)= 2*pi*ex/d;
    end
    % escape time C to A
    for i=1:5
        pt=20*i ;
        V = pt/xa ;
        tca(i) = xd(i+5,1) ; % SSA, escape time A to C at T(i), translation compensated
        eca(i) = xd(i+5,2)  ; % SSA, error bar 
        gp=log(2.); % gamma_p
        gt= 1. ;  % 1. / T
        d=gp*sqrt(W2d(xc,S0,Sa,V,gt,K)*abs(V2d(xb,S0,Sa,V,gt,K))) ;
        ex=exp( Vd(xb,S0,Sa,V,gt,K) - Vd(xc,S0,Sa,V,gt,K) - log(ttd1(xc,V,gt,K))) ;  
        ttc(i)= 2*pi*ex/d;  % Theoretical escape time C to A , translation compensated
        
        % constant theta
        kt=xc*2./V  ;           % theta, transcription compensated
        d=gp*sqrt(uzc*abs(uzb)) ;  
        ex=exp((ub-uc)/kt) ; 
        ttcc(i)= 2*pi*ex/d;
    end
    % plots
    subplot(2,2,1)
    hold on  
    loglog(tr,tta, 'k') ;
    loglog(tr,ttac, 'r') ;
    h1 = errorbar(tr,tac,eac,'kx');
    set(get(h1,'Parent'), 'YScale', 'log')
    set(get(h1,'Parent'), 'XScale', 'log')
    xlabel('P_{T, low}','FontSize',16)
    ylabel('\Gamma_{ low \rightarrow high}','FontSize',16)
    xlim([10 110])
    
    subplot(2,2,3)
    hold on
    loglog(tr,ttc, 'k')   ;
    loglog(tr,ttcc, 'r')   ;
    h2 = errorbar(tr,tca,eca,'kx');
    set(get(h2,'Parent'), 'YScale', 'log')
    set(get(h2,'Parent'), 'XScale', 'log')  
    xlabel('P_{T, low}','FontSize',16)
    ylabel('\Gamma_{ high \rightarrow low}','FontSize',16)
    xlim([10 110])
cd ..

%   *****************************************************************
%   ***************************** Dimer   ***********************
S0=0.11;
Sa=3.3;
al=1.0;
gpr= 10 ; 
K=1.;
K1=2.*K/(al.*gpr);
bt= 1./gpr-1.0;
KK=K1/K;
Kr=K/K1;

cd fnD
xl = fig3b;
xh = fig3d;


pp=[-al -1. S0+Sa-al -1. S0];
r = roots(pp); % calculation of roots of first derivative of potential 
xa= r(4)+2.*(r(4))^2*Kr;  % concentration at stable point A
xb= r(3)+2.*(r(3))^2*Kr; % concentration at unstable point B
xc= r(2)+2.*(r(2))^2*Kr;  % concentration at stable point C

    for i=1:5
        % SSA
        tac(i) = xl(tc2(i),1) ; % SSA, escape time A to C at T(i), transC
        eac(i) = xl(tc2(i),2)  ; % SSA, error bar 
        pt=tr2(i) ;
        V = pt/xa ;
        gp=gpr*log(2.)  ;  % gamma_p1
        gt= 1. ;
        % FP integration
        d=gp*sqrt(dwdi(xa,S0,Sa,V,gt,bt,K,K1)*abs(dvdi(xb,S0,Sa,V,gt,bt,K,K1))) ;
        ex=exp( vddm(xb,S0,Sa,V,gt,bt,K,K1) - vddm(xa,S0,Sa,V,gt,bt,K,K1) - log(ttdi(xa,V,gt,bt,K,K1))) ;  % numerator of tau formula 
        tta(i)= 2*pi*ex/d;  % escape time C to A , translation compensated
        % Constant theta
        vwd = integral(@(p) fund(p,S0,Sa,bt,KK,Kr),xa,xb); 
        ei=exp(vwd/ttdi(xa,V,gt,bt,K,K1)) ;  % numerator of tau formula 
        d=gp*sqrt(dvd(xa,S0,Sa,bt,KK,Kr)*abs(dvd(xb,S0,Sa,bt,KK,Kr))) ;
        ttac(i)= 2.*pi*ei/d; 
    end
    for i=1:5
        % SSA
        tca(i) = xh(i,1) ; % SSA, escape time A to C at T(i), transC
        eca(i) = xh(i,2)  ; % SSA, error bar  
        pt=20*i ;
        V = pt/xa ;
        gp=gpr*log(2.)  ;  % gamma_p1
        gt= 1. ;
        % FP integration
        d=gp*sqrt(dwdi(xc,S0,Sa,V,gt,bt,K,K1)*abs(dvdi(xb,S0,Sa,V,gt,bt,K,K1))) ;
        ex=exp( vddm(xb,S0,Sa,V,gt,bt,K,K1) - vddm(xc,S0,Sa,V,gt,bt,K,K1) - log(ttdi(xc,V,gt,bt,K,K1))) ;  % numerator of tau formula 
        ttc(i)= 2*pi*ex/d;  % escape time C to A 
        % Constant theta
        vwd = integral(@(p) fund(p,S0,Sa,bt,KK,Kr),xc,xb); 
        ei=exp(vwd/ttdi(xc,V,gt,bt,K,K1)) ;  % numerator of tau formula 
        d=gp*sqrt(dvd(xc,S0,Sa,bt,KK,Kr)*abs(dvd(xb,S0,Sa,bt,KK,Kr))) ;
        ttcc(i)= 2.*pi*ei/d;
    end
    % plots
    subplot(2,2,2)
    hold on  
    loglog(tr2,tta, 'b') ;
    loglog(tr2,ttac, 'r') ;
    h1 = errorbar(tr2,tac,eac,'bx');
    set(get(h1,'Parent'), 'YScale', 'log')
    set(get(h1,'Parent'), 'XScale', 'log')
    xlabel('P_{T, low}','FontSize',16)
    ylabel('\Gamma_{ low \rightarrow high}','FontSize',16)
    xlim([10 310])
    
    subplot(2,2,4)
    hold on
    loglog(tr,ttc, 'b')   ;
    loglog(tr,ttcc, 'r')   ;
    h2 = errorbar(tr,tca,eca,'bx');
    set(get(h2,'Parent'), 'YScale', 'log')
    set(get(h2,'Parent'), 'XScale', 'log')  
    xlabel('P_{T, low}','FontSize',16)
    ylabel('\Gamma_{ high \rightarrow low}','FontSize',16)
    xlim([10 110])
cd ..
h = legend('Theoretical, interpolated \theta','Theoretical, constant \theta','SSA');
