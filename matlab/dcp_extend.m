%% dcp_extend.m - 단일 평행판 -> DCP(이중복합) 확장, d_c=0.5 마진 비교
%  행렬법(검증됨)으로 단일 K_p 산출. DCP는 2자유도(중간단+출력) 모달.
%   K_eff(DCP)=K_p (직렬x병렬 상쇄), 각 단 x/2 변형 -> 응력 절반, 중간단 질량 추가
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; P.Jmax=10.2; dc=0.5e-3;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);
DB={'Mg',45,200,1740;'Al6061',69,276,2700;'Al7075',72,503,2810;'Ti64',114,880,4430;
    'BeCu',128,1100,8250;'17-4PH',197,1170,7800;'Spring',200,1200,7850;'Maraging',190,2000,8100};

fprintf('=== d_c=0.5 : 단일 평행판 vs DCP(이중복합) gamma 비교 ===\n');
fprintf('%-9s%14s%14s\n','재료','gamma(단일)','gamma(DCP)');
fprintf('%s\n',repmat('-',1,38));
for k=1:8
  m.E=DB{k,2}*1e9;m.Sy=DB{k,3}*1e6;m.rho=DB{k,4};
  gS=solveX(m,dc,P,opts,14,@conSingle);
  gD=solveX(m,dc,P,opts,14,@conDCP);
  fprintf('%-9s%10.3f%s%10.3f%s\n',DB{k,1}, gS,tag(gS), gD,tag(gD));
end

%% ===== 함수 =====
function s=tag(g), if g>=1.0,s=' OK';else,s='   ';end, end
function K=hingeK(L,t,b,th,s,E)
  A=b*t;Iz=b*t^3/12; Ch=[L/(E*A),0,0;0,L^3/(3*E*Iz),L^2/(2*E*Iz);0,L^2/(2*E*Iz),L/(E*Iz)];
  c=cos(th);sn=sin(th);R=[c,-sn,0;sn,c,0;0,0,1];T=[1,0,-s(2);0,1,s(1);0,0,1]; K=inv(T*R*Ch*R'*T');
end
function [Kp,f1s,info]=flexA(L,t,b,W,H,mcoil,E,rho)
  dy=H/2-t/2; sL={[W/2;-dy],[W/2;dy],[-W/2;-dy],[-W/2;dy]}; th=[0,0,pi,pi];
  Kt=zeros(3);Ki=cell(4,1); for i=1:4,Ki{i}=hingeK(L,t,b,th(i),sL{i},E);Kt=Kt+Ki{i};end
  mB=rho*W*b*H; M=diag([mB+mcoil,mB+mcoil,mB*(W^2+H^2)/12]);
  Kp=Kt(2,2); fn=sort(sqrt(abs(eig(Kt,M)))/(2*pi)); f1s=fn(1);
  info.Ki={Ki};info.sL={sL};info.th=th;info.L=L;info.t=t;info.b=b;info.Iz=b*t^3/12;
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
  A_g=2*pi*(rmM+tgM/2)*hyM;v.By1=Bg*A_g/(pi*rmM^2);v.VCMd=2*(rmM+tgM+ty3);v.wc=wc;v.dcs=dcs;
end
% ---- 단일 평행판 ----
function [c,ceq]=conSingle(x,m,dc,P)
  L=x(1)/1e3;t=x(2)/1e3;b=x(3)/1e3;a=x(4)/1e3; v=vcm(x(5),x(6),x(7),dc,P);
  [Kp,f1,info]=flexA(L,t,b,a,a,v.mcoil,m.E,m.rho);
  xs=v.Kf*P.I/Kp; g=x(8); J=P.I/((pi/4)*dc^2)/1e6; sig=flexStress(info,g*1e-3,P.KT);
  c=[g*1e-3-xs; g*100-f1; sig-P.SIG*m.Sy; v.VCMd-20e-3; (a+2*L)-200e-3; v.dcs-v.wc; v.By1-2.2; J-P.Jmax];
  ceq=[];
end
% ---- DCP (이중복합, 2자유도) ----
function [c,ceq]=conDCP(x,m,dc,P)
  L=x(1)/1e3;t=x(2)/1e3;b=x(3)/1e3;a=x(4)/1e3; v=vcm(x(5),x(6),x(7),dc,P);
  [Kp,~,info]=flexA(L,t,b,a,a,0,m.E,m.rho);
  Keff=Kp; Mout=m.rho*a*b*a+v.mcoil; Mi=m.rho*a*b*a;   % 중간단 총질량 ~ 출력 구조질량
  K2=[4*Kp,-2*Kp;-2*Kp,2*Kp]; M2=diag([Mi,Mout]);
  fn=sort(sqrt(abs(eig(K2,M2)))/(2*pi)); f1=fn(1);
  xs=v.Kf*P.I/Keff; g=x(8); J=P.I/((pi/4)*dc^2)/1e6;
  sig=flexStress(info,(g*1e-3)/2,P.KT);              % 각 단 x/2 -> 응력 절반
  c=[g*1e-3-xs; g*100-f1; sig-P.SIG*m.Sy; v.VCMd-20e-3; (a+4*L)-200e-3; v.dcs-v.wc; v.By1-2.2; J-P.Jmax];
  ceq=[];
end
function g=solveX(m,dc,P,opts,ns,cf)
  lb=[10,0.2,5,15,5,2,1.2,0.1]; ub=[70,1.5,20,50,8,7,2.5,3]; x0=[30,0.4,10,25,7,4,1.5,0.6];
  bg=-Inf; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)cf(x,m,dc,P),opts);
      if ef>0,c=cf(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;end,end
    catch
    end
  end
  g=bg;
end
