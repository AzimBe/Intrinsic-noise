clear
clc
% S0=0.135;
% fn=1.0/(1.+1.8/0.135);   % inverse of f
% x= [0.001:0.0001:1.0];  % concentration
% for i=1:9991
%     t=x(i);
% %     g=(fn+t)/(1+t);
% %     gp=(1-fn)/(1+t)^2; % first derivative of g
% %     bt(i)=2*t/(2*t*gp-g);  % beta
% %     al(i)=((bt(i)*g-2*t)/(sqrt(t)))^2;  % alhpa
%     f(i)=(-sqrt(t)*(1+t)^2/2.+S0+2.*S0*t)/(S0*t^2);
%     al(i)=2.*t*(1+t)/(S0*(1+t*f(i))-sqrt(t)*(1+t));
% end
% 
% 
% loglog(al,f)
% hold on
% % axis([0,11000,0,300])
% % plot(1000, 64, '+k')
% % plot(10000,200,'+r')
% % plot(0.0, 0.0, 'ob')
% xlabel('\alpha','FontSize',18)
% ylabel('\beta','FontSize',18)
% title('stability phase plot','FontSize',18)
% % h = legend('stability regipn','\alpha = 1000, \beta = 64','\alpha = 10000, \beta = 200',3);
% z=[0.01 0.001 0.0001];
% subplot(2,2,4)
hold on
for i=1:1
    a=20.0;
    syms t
    Sa = (1.+2.*a*t)*(1+t^2)^2/(2.*t);
    S0 = t+a*t^2-t*(1.+2.*a*t)*(1+t^2)/2. ;
    ezplot(S0 , Sa)
    axis square
    hold on
end
% xlim([0.0 0.4])
% ylim([0.9 2.4])
% xlim([0.0 0.2])
% ylim([7.5 9.])
% xlim([0.0 0.4])
% ylim([5.0 10])
xlim([0.0 0.1])
ylim([24. 28.])
xlabel('S_0','FontSize',12)
ylabel('S_a','FontSize',12)
title( '\alpha = 20.0 ','FontSize',16 )

t0=0.03;
ta=26.5;
plot(t0, ta, '+r')
% plot(0.065, 10., '+r')
V=93.72;
K1=0.1;
K=1.0;
p=[-a -1 t0+ta-a -1 t0];
r = roots(p);
P1= r*K*V
P2= r1.^2/(K1*V)
Pt= r1+2*r2
b=V*K1/8.