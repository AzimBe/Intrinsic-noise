%   *****************************************************************
%   ***************************** Cooerative   ***********************
clear
clc
figure()
S0=0.135;  % S_0
Sa= 1.8;    % S_a
K=1.;

cd fnM

p= [1  -(S0+Sa) 1 -S0];
r=roots(p);  % roots of first derivative of potential 
xa=r(3);  % concentration at stable point A
xb=r(2);  % concentration at unstable point B
xc=r(1);  % concentration at stable point C
tr = [  1.; sqrt(2);2 ;2.*sqrt(2.);4];
 
    y(1:5,1:2) = MtransC  ; % SSA result transC
    y(1:5,3:4) = MtransL  ; % SSA result transL

    V= 40./xa ;
    for i=1:5
        tac(i) = y(i,1) ; 
        eac(i) = y(i,2)  ; 
        tacr(i) = y(i,3) ; 
        eacr(i) = y(i,4)  ; 
        gp=log(2.)/tr(i)  ; 
        gt= 1./tr(i) ;  
        d=gp*sqrt(W2d(xc,S0,Sa,V,gt,K)*abs(V2d(xb,S0,Sa,V,gt,K))) ;
        dr=gp*sqrt(W2r(xc,S0,Sa,V,K)*abs(V2r(xb,S0,Sa,V,K))) ;
        ex=exp( Vd(xb,S0,Sa,V,gt,K) - Vd(xc,S0,Sa,V,gt,K) - log(ttd1(xc,V,gt,K))) ;  
        exr=exp(Vr(xb,S0,Sa,V,K) - Vr(xc,S0,Sa,V,K) - log(ttr1(xc,V,K))) ;
        tta(i)= 2*pi*ex/d;  % Theoretical escape time C to A , translation compensated
        ttar(i)= 2*pi*exr/dr; % Theoretical escape time C to A , transcription compensated
    end
    % plots
    subplot(1,2,1)
    hold on
    
    loglog(tr,tta(1:5), 'k');
    errorbar(tr,tac,eac,'kx');
    loglog(tr,ttar(1:5), 'r');
    h1 = errorbar(tr,tacr,eacr,'rx');
    set(get(h1,'Parent'), 'YScale', 'log')
    set(get(h1,'Parent'), 'XScale', 'log')

    xlabel('generation time','FontSize',16)
    ylabel('\Gamma_{ high \rightarrow low}','FontSize',16)
%     title('cooperative, S_0 = 0.135, S_a = 1.8','FontSize',16)
%     axis square
%     xlim([1 4.2])
cd ..

%   *****************************************************************
%   ***************************** Dimer   ***********************
al= 1.;
S0=0.11;
Sa=3.3;
gpr = 1. ;
K=1.;
K1=2.*K/(al.*gpr);
bt= 1./gpr - 1.0 ;
KK=K1/K;
Kr=K/K1;


cd fnD

pp=[-al -1 S0+Sa-al -1 S0];
r = roots(pp); % calculation of roots of first derivative of potential 
xa= r(4)+2.*(r(4))^2*Kr;  % concentration at stable point A
xb= r(3)+2.*(r(3))^2*Kr; % concentration at unstable point B
xc= r(2)+2.*(r(2))^2*Kr;  % concentration at stable point C

    x(1:5,1:2) = DtransC  ; % SSA result transC
    x(1:5,3:4) = DtransL  ; % SSA result transL
    V= 40./xa ;

    for i=1:5
        tca(i) = x(i,1) ; % log of SSA, 2 comp. escape time C to A at T(i)
        eca(i) = x(i,2)  ; % SSA, 2 comp. error bar at T(i)
        tcar(i) = x(i,3) ; % log of SSA, 2 comp. escape time C to A at T(i)
        ecar(i) = x(i,4)  ; % SSA, 2 comp. error bar at T(i)
        gp=gpr*log(2.)/tr(i)  ;  % gamma_p
        gt= 1./tr(i) ;
        d=gp*sqrt(dwdi(xc,S0,Sa,V,gt,bt,K,K1)*abs(dvdi(xb,S0,Sa,V,gt,bt,K,K1))) ;
        dr=gp*sqrt(dwri(xc,S0,Sa,V,gt,bt,K,K1)*abs(dvri(xb,S0,Sa,V,gt,bt,K,K1))) ;
        ex=exp( vddm(xb,S0,Sa,V,gt,bt,K,K1) - vddm(xc,S0,Sa,V,gt,bt,K,K1) - log(ttdi(xc,V,gt,bt,K,K1))) ;  % numerator of tau formula , Eq. (11)
        exr=exp( vdrm(xb,S0,Sa,V,gt,bt,K,K1) - vdrm(xc,S0,Sa,V,gt,bt,K,K1) - log(ttri(xc,V,bt,K,K1))) ;
        ttc(i)= 2.*pi*ex/d;  % Theoretical escape time C to A , Eq. (11)
        ttcr(i)= 2.*pi*exr/dr;
        % numerical integration
%         vwd = integral(@(p) fund(p,S0,Sa,V,gt,bt,K,K1),xc,xb); 
%         vwr = integral(@(p) funr(p,S0,Sa,V,gt,bt,K,K1),xc,xb);
%         ei=exp(vwd-log(ttdi(xc,V,gt,bt,K,K1))) ;  % numerator of tau formula , Eq. (11)
%         eir=exp(vwr-log(ttri(xc,V,bt,K,K1))) ;
%         d=gp*sqrt(dwdi(xc,S0,Sa,V,gt,bt,K,K1)*abs(dvdi(xb,S0,Sa,V,gt,bt,K,K1))) ;
%         dr=gp*sqrt(dwri(xc,S0,Sa,V,gt,bt,K,K1)*abs(dvri(xb,S0,Sa,V,gt,bt,K,K1))) ;
%         tic(i)= 2.*pi*ei/d;  % Theoretical escape time C to A , Eq. (11)
%         ticr(i)= 2.*pi*eir/dr;
    end
    subplot(1,2,2)
    hold on
    
    loglog(tr,ttc(1:5), 'k')
%     loglog(tr,tic(1:5), '-.k')    
    errorbar(tr,tca,eca,'kx');
    loglog(tr,ttcr(1:5), 'r')
%     loglog(tr,ticr(1:5), '-.r')
    h2 = errorbar(tr,tcar,ecar,'rx');
    set(get(h2,'Parent'), 'YScale', 'log')
    set(get(h2,'Parent'), 'XScale', 'log')
    

    xlabel('generation time','FontSize',16)
    ylabel('\Gamma_{ high \rightarrow low}','FontSize',16)
%     title('dimer, S_0 = 0.11, S_a = 3.3, \alpha = 1., \gamma_p = \gamma_p_2, K=1.','FontSize',16)
%     axis square
% h = legend('SSA, transl.','Theoretical, transl.','Num int, transl.','constant \theta, transl.','SSA, transc.','Theoretical, transc.','Num int, transc.','constant \theta, transc.');
cd ..