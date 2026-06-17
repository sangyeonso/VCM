%% codesign_current.m - 전류 I를 올려 전 KPI 충족 + 발열 대가 정량화
%  단일 VCM, DCP+Al7075, d_c=0.5. I 스윕. J<=10 제약은 풀고 J·P를 보고.
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.SIG=0.2; P.KT=1.7; P.mi=0.5; P.rho_e=1.68e-8;
m.E=72e9; m.Sy=503e6; m.rho=2810;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',300,...
     'MaxFunctionEvaluations',3000,'ConstraintTolerance',1e-9);

fprintf('=== 전류 I 상향에 따른 γ + 발열 (단일 VCM, DCP+Al7075, d_c=0.5) ===\n');
fprintf('%-7s%8s%11s%9s%9s%8s%9s%8s\n','I[A]','gamma','stroke[mm]','f1[Hz]','J[A/mm2]','R[Ω]','P=I2R[W]','판정');
fprintf('%s\n',repmat('-',1,68));
for I=[2.0 2.3 2.5 2.7 3.0]
  P.I=I;
  [g,r]=solveDCP(m,P,opts,18);
  J=I/((pi/4)*(0.5e-3)^2)/1e6; Pw=I^2*r.R;
  vd='미충족'; if g>=1.0, vd='충족 OK'; end
  fprintf('%-7.1f%8.3f%11.3f%9.1f%9.1f%8.3f%9.2f%8s\n', I, g, r.x*1e3, r.f1, J, r.R, Pw, vd);
end
fprintf('\n주: J<=10(연속 공랭 경험치) 기준 — I>2A는 J 초과라 냉각 강화 필요. 단 총 P는 낮음(저저항 코일).\n');
fprintf('    원 스펙 "최대 2A"를 푸는 것과 동치. 듀얼 VCM과 비교: 앰프 1개(고전류) vs 앰프 2개.\n');

%% ===== 함수 =====
function K=hingeK(L,t,b,th,s,E)
  A=b*t;Iz=b*t^3/12; Ch=[L/(E*A),0,0;0,L^3/(3*E*Iz),L^2/(2*E*Iz);0,L^2/(2*E*Iz),L/(E*Iz)];
  c=cos(th);sn=sin(th);R=[c,-sn,0;sn,c,0;0,0,1];T=[1,0,-s(2);0,1,s(1);0,0,1]; K=inv(T*R*Ch*R'*T');
end
function [Kp,info]=flexA(L,t,b,W,H,E)
  dy=H/2-t/2; sL={[W/2;-dy],[W/2;dy],[-W/2;-dy],[-W/2;dy]}; th=[0,0,pi,pi];
  Kt=zeros(3);Ki=cell(4,1); for i=1:4,Ki{i}=hingeK(L,t,b,th(i),sL{i},E);Kt=Kt+Ki{i};end
  Kp=Kt(2,2); info.Ki={Ki};info.sL={sL};info.th=th;info.L=L;info.t=t;info.b=b;info.Iz=b*t^3/12;
end
function sig=flexStress(info,defl,KT)
  u=[0;defl;0];Ki=info.Ki{1};sL=info.sL{1};th=info.th;L=info.L;t=info.t;b=info.b;Iz=info.Iz; sig=0;
  for i=1:4
    fb=Ki{i}*u;s=sL{i};Tc=[1,0,0;0,1,0;-s(2),s(1),1];fh=Tc*fb;
    R=[cos(th(i)),-sin(th(i)),0;sin(th(i)),cos(th(i)),0;0,0,1];fl=R'*fh;
    Mf=abs(fl(2)*L+fl(3));sg=Mf*(t/2)/Iz+abs(fl(1))/(b*t); if sg>sig,sig=sg;end
  end
  sig=KT*sig;
end
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3;wc=tgM-2*c_c;dcs=rp*dc;hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs));L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n=nsym*n_half;n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half,n);
  v.Kf=n_eff*Bg*L_m;v.mcoil=8960*n*(pi/4)*dc^2*L_m;
  v.R=P.rho_e*n*L_m/((pi/4)*dc^2);
  A_g=2*pi*(rmM+tgM/2)*hyM;v.By1=Bg*A_g/(pi*rmM^2);v.VCMd=2*(rmM+tgM+ty3);v.wc=wc;v.dcs=dcs;
end
function r=perfDCP(x,m,P)
  L=x(1)/1e3;t=x(2)/1e3;b=x(3)/1e3;a=x(4)/1e3; dc=0.5e-3; v=vcm(x(5),x(6),x(7),dc,P);
  [Kp,info]=flexA(L,t,b,a,a,m.E); Keff=Kp;
  Mout=m.rho*a*b*a+v.mcoil; Mi=P.mi*m.rho*a*b*a;
  K2=[4*Kp,-2*Kp;-2*Kp,2*Kp]; M2=diag([Mi,Mout]);
  fn=sort(sqrt(abs(eig(K2,M2)))/(2*pi)); f1=fn(1);
  xs=v.Kf*P.I/Keff; g=x(8);
  r.Keff=Keff;r.f1=f1;r.x=xs;r.Kf=v.Kf;r.mcoil=v.mcoil;r.R=v.R;r.By1=v.By1;r.VCMd=v.VCMd;
  r.dcs=v.dcs;r.wc=v.wc;r.info=info;r.foot=a+4*L; r.sig=flexStress(info,(g*1e-3)/2,P.KT);
end
function [c,ceq]=conDCP(x,m,P)
  r=perfDCP(x,m,P); g=x(8);   % J 제약 제거(전류 상향 탐색)
  c=[g*1e-3-r.x; g*100-r.f1; r.sig-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end
function [g,r]=solveDCP(m,P,opts,ns)
  lb=[10,0.2,5,15,5,2,1.2,0.1]; ub=[70,1.5,20,50,8,7,2.5,3]; x0=[20,0.3,8,25,7,4,1.5,0.6];
  bg=-Inf;bx=x0; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conDCP(x,m,P),opts);
      if ef>0,c=conDCP(xo,m,P); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;bx=xo;end,end
    catch
    end
  end
  g=bg; r=perfDCP(bx,m,P);
end
