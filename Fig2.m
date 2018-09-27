clear
clc
figure()
%   *****************************************************************
%   ***************************** bifurcation curves   ***********************
subplot(2,2,1)  % 
h=[2 2 2];  % hill function coefficient
al=[1. 0.5 0.0];  % alpha 
CArray = {'Color'}; VCArray = {'b','r','k'}';
LArray = {'LineStyle'}; VLArray = {'-','-','-'}';
S0V=[0.11  0.125 0.135];
SaV=[3.3  2.55 1.8];
tc=['b+';'r+';'k+'];

for i=1:3
    n=h(i);
    a=al(i);
    syms t
    Sas = (1.+2.*a*t)*(1+t^n)^2/(n*t^(n-1));  % Stationary S0
    S0s = t+a*t^2-t*(1.+2.*a*t)*(1+t^n)/n ;   % Stationary Sa
    h1 = ezplot(S0s , Sas) ;
    set(h1, LArray,VLArray(i),CArray,VCArray(i))
    axis square
    hold on
    S0=S0V(i); Sa=SaV(i); plot(S0, Sa, tc(i,:),'MarkerSize',10)
end

for i=1:1300      
    x(i)=0.135*i/1000. ;
    y(i)=1.8*i/1000. ;
end
plot(x,y,'k-.')  % constant fold change line

plot(0.1289, 1.7188 , '+m','MarkerSize',10)   % doubling rate 1.04
plot(0.1430, 1.9060 , '+g','MarkerSize',10)   % doubling rate 0.94

xlim([0.0 0.25])
ylim([1. 4.])
xlabel('S_0','FontSize',12)
ylabel('S_a','FontSize',12)
title( '         ','FontSize',16 )
h = legend('\alpha = 1.0','fold change 31.','\alpha = 0.5' ,'fold change 21.4','\alpha = 0.0', 'fold change 14.3','constant fold change', 'doubling rate 1.04', 'doubling rate 0.94');

%   *****************************************************************
%   ***************************** alpha = 0   ***********************
subplot(2,2,2)
S0=S0V(3); Sa=SaV(3);V=187.44;
p= [1  -(S0+Sa) 1 -S0];
r=roots(p);  % roots of first derivative of potential 
xa=r(3);  % concentration at stable point A
xb=r(2);  % concentration at unstable point B
xc=r(1);  % concentration at stable point C
% pb = round(xb*40./xa);
for i=1:500
    phi(i)=( -(S0+Sa)*(i/V)+Sa*atan(i/V)+(i/V)^2/2. ) -  ( -(S0+Sa)*(xb)+Sa*atan(xb)+(xb)^2/2. ) ;
%     phi(i)=S0+Sa*(i/V)^2/(1.+(i/V)^2)-(i/V);
end
plot(phi,'k')
axis square 
xlabel('P_T','FontSize',12)
ylabel('\int \Phi','FontSize',12)
ylim([-0.1 0.02])
hold on
subplot(2,2,3)
for i=1:800
    tt(i)=2.*(i/V)/V ;
end
plot(tt,'k')
axis square 
xlabel('P_T','FontSize',12)
ylabel('\theta','FontSize',12)
ylim([0. 0.08])
hold on
subplot(2,2,4)
gt = 1. ;
for i=1:600
    psi(i)= ((-S0*log(i/V)-Sa*log(1+(i/V)^2)/2.+(i/V))*V/2)- ((-S0*log(xb)-Sa*log(1+(xb)^2)/2.+(xb))*V/2) ;
end
plot(psi,'k')
axis square
xlabel('P_T','FontSize',12)
ylabel('\Psi','FontSize',12)
ylim([-3. 0.5])
hold on

cd fnD

%   *******************************************************************
%   ***************************** alpha = 0.5   ***********************
subplot(2,2,2)
S0=S0V(2); Sa=SaV(2);a = al(2); 
gpr= 10 ; 
K=1.;
K1=2.*K/(a*gpr);
KK=K1/K;
Kr=K/K1;
bt= 1./gpr-1.0;
pp=[-a -1. S0+Sa-a -1. S0];
r = roots(pp); % calculation of roots of first derivative of potential 
xa= r(4)+2.*(r(4))^2*Kr;
xb= r(3)+2.*(r(3))^2*Kr; % concentration at unstable point B
xc= r(2)+2.*(r(2))^2*Kr;  % concentration at stable point C
V = 40./xa ;
for i=1:800
    phi(i)= integral(@(p) fund(p,S0,Sa,bt,KK,Kr),xb,i/V);
end
plot(phi,'r')
subplot(2,2,3)
for i=1:800
    tt(i)= ttdi(i/V,V,gt,bt,K,K1) ; 
end
plot(tt,'r')
subplot(2,2,4)
for i=1:800
    psi(i) = vddm(i/V,S0,Sa,V,gt,bt,K,K1) - vddm(xb,S0,Sa,V,gt,bt,K,K1) ; 
end
plot(psi,'r')

%   *******************************************************************
%   ***************************** alpha = 1.0   ***********************
subplot(2,2,2)
S0=S0V(1); Sa=SaV(1);a = al(1);
gpr= 10 ; 
K=1.;
K1=2.*K/(a*gpr);
KK=K1/K;
Kr=K/K1;
bt= 1./gpr-1.0;
pp=[-a -1. S0+Sa-a -1. S0];
r = roots(pp); % calculation of roots of first derivative of potential 
xa= r(4)+2.*(r(4))^2*Kr;
xb= r(3)+2.*(r(3))^2*Kr; % concentration at unstable point B
xc= r(2)+2.*(r(2))^2*Kr;  % concentration at stable point C
V = 40./xa ;
for i=1:800
    phi(i)= integral(@(p) fund(p,S0,Sa,bt,KK,Kr),xb,i/V);
end
plot(phi,'b')
subplot(2,2,3)
for i=1:800
    tt(i)= ttdi(i/V,V,gt,bt,K,K1) ; 
end
plot(tt,'b')
subplot(2,2,4)
gt = 1. ;
for i=1:800
    psi(i) = vddm(i/V,S0,Sa,V,gt,bt,K,K1) - vddm(xb,S0,Sa,V,gt,bt,K,K1) ; 
end
plot(psi,'b')
cd ..

