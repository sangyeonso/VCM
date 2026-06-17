%% codesign_beta.m - 발열(J<=10) + d_c=0.5(스펙)에서 B_g 배율 beta 스윕
%  "gamma>=1 되려면 B_g를 몇 배(beta) 키워야 하나?" -> Halbach 필요수준 도출
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; P.rho_e=1.68e-8; P.Jmax=10.2; dc=0.5e-3;  % d_c=0.5,I=2A -> J=10.2 (스펙 열한계)
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);
DB={'Mg',45,200,1740;'Al6061',69,276,2700;'Al7075',72,503,2810;'Ti64',114,880,4430;
    'BeCu',128,1100,8250;'17-4PH',197,1170,7800;'Spring',200,1200,7850;'Maraging',190,2000,8100};

fprintf('=== d_c=0.5(스펙) + J<=10(발열) : B_g 배율 beta 스윕 ===\n');
fprintf('%-8s%10s%14s\n','beta','best gamma','best 재료');
fprintf('%s\n',repmat('-',1,34));
for beta=[1.0 1.2 1.4 1.6 1.8 2.0]
  P.beta=beta; bestg=-Inf; bn='';
  for k=1:8
    m.E=DB{k,2}*1e9;m.Sy=DB{k,3}*1e6;m.rho=DB{k,4};
    g=solveleaf(m,dc,P,opts,12);
    if g>bestg, bestg=g; bn=DB{k,1}; end
  end
  vd=''; if bestg>=1.0, vd=' <= 충족!'; end
  fprintf('%-8.1f%10.3f%14s%s\n', beta, bestg, bn, vd);
end
fprintf('\n주: beta=1.0 은 현재 N52 단일자석. Halbach/자속집속이 beta 몇까지 주는지가 관건.\n');

%% ===== 함수 =====
function g=solveleaf(m,dc,P,opts,ns)
  lb=[0.3,15,5,15,5,2,1.2,0.1]; ub=[2,60,20,45,8,7,2.5,3]; x0=[0.6,35,10,25,7,4,1.5,0.7];
  bg=-Inf; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conL(x,m,dc,P),opts);
      if ef>0,c=conL(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;end,end
    catch
    end
  end
  g=bg;
end
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=P.beta*(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;   % beta 배율
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3; wc=tgM-2*c_c; dcs=rp*dc; hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs)); L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n=nsym*n_half; n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half, n);
  v.Kf=n_eff*Bg*L_m; v.mcoil=8960*n*(pi/4)*dc^2*L_m;
  A_g=2*pi*(rmM+tgM/2)*hyM; v.By1=Bg*A_g/(pi*rmM^2); v.VCMd=2*(rmM+tgM+ty3); v.wc=wc; v.dcs=dcs;
end
function r=perfL(tm,Lm,bm,am,rm,hy,tg,E,rho,dc,P)
  t=tm/1e3;L=Lm/1e3;b=bm/1e3;a=am/1e3; v=vcm(rm,hy,tg,dc,P); db=40e-3;hi=5e-3;
  Keff=2*E*b*t^3/L^3; Mflex=rho*a^2*b+0.25*(2*rho*db*hi*b)+(1/3)*(8*rho*L*t*b);
  Meff=Mflex+v.mcoil; xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.x=xs;r.f1=f1;r.By1=v.By1;r.VCMd=v.VCMd;r.dcs=v.dcs;r.wc=v.wc;r.t=t;r.L=L;r.foot=a+2*(2*L+hi);
end
function [c,ceq]=conL(x,m,dc,P)
  r=perfL(x(1),x(2),x(3),x(4),x(5),x(6),x(7),m.E,m.rho,dc,P); g=x(8); xop=g*1e-3;
  peak=P.KT*3*m.E*r.t*(xop/2)/r.L^2; J=P.I/((pi/4)*dc^2)/1e6;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2; J-P.Jmax];
  ceq=[];
end
