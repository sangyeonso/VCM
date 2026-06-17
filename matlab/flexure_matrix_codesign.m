%% flexure_matrix_codesign.m
%  유연기구를 강의 정석(컴플라이언스 행렬법, homework1-3)으로 재구축 + VCM 공동최적+발열
%  (A) homework2 재현 검증 (6043Hz, 16.69um)
%  (B) 행렬법 유연기구로 d_c=0.5(스펙)·J<=10.2·8재료 공동최적 -> gamma 재판정
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; P.Jmax=10.2;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);

%% (A) 검증: homework2 (L4 t0.3 b5, body 6x5, Al)
[Ke,f1,Me,~]=flexA(4e-3,0.3e-3,5e-3,6e-3,5e-3,0,71e9,2770);
fprintf('=== (A) 행렬법 검증 (homework2) ===\n');
fprintf('Kyy=%.3e N/m (정답 5.99e5), f1=%.0f Hz (정답 6043), dy@10N=%.2f um (16.69)\n\n',...
   Ke, f1, 10/Ke*1e6);

%% (B) 행렬법 유연기구 공동최적 (d_c=0.5, J<=10.2)
DB={'Mg',45,200,1740;'Al6061',69,276,2700;'Al7075',72,503,2810;'Ti64',114,880,4430;
    'BeCu',128,1100,8250;'17-4PH',197,1170,7800;'Spring',200,1200,7850;'Maraging',190,2000,8100};
fprintf('=== (B) 행렬법 유연기구 + VCM 공동최적 (d_c=0.5, J<=10.2, 피로0.2Sy) ===\n');
fprintf('%-9s%9s%11s%9s%9s%10s\n','재료','gamma','stroke[mm]','f1[Hz]','Kf','Keff[N/mm]');
fprintf('%s\n',repmat('-',1,60));
best=-Inf; bn='';
for k=1:8
  m.E=DB{k,2}*1e9;m.Sy=DB{k,3}*1e6;m.rho=DB{k,4};
  [bx,r]=solveCD(m,0.5e-3,P,opts,14);
  if bx(8)>best,best=bx(8);bn=DB{k,1};end
  vd=''; if bx(8)>=1.0,vd=' <=충족';end
  fprintf('%-9s%9.3f%11.3f%9.1f%9.2f%10.1f%s\n',DB{k,1},bx(8),r.x*1e3,r.f1,r.Kf,r.Keff/1e3,vd);
end
fprintf('%s\n',repmat('-',1,60));
fprintf('d_c=0.5 최고: %s gamma=%.3f -> %s\n', bn, best, tern(best>=1,'충족!','미충족'));
fprintf('(이전 단순식 24EI/L^3 결과는 0.92 였음 -> 행렬법으로 바뀌나 비교)\n');

%% ===== 함수 =====
function s=tern(c,a,b), if c,s=a;else,s=b;end, end
function K=hingeK(L,t,b,th,s,E)
  A=b*t; Iz=b*t^3/12;
  Ch=[L/(E*A),0,0; 0,L^3/(3*E*Iz),L^2/(2*E*Iz); 0,L^2/(2*E*Iz),L/(E*Iz)];
  c=cos(th);sn=sin(th); R=[c,-sn,0;sn,c,0;0,0,1]; T=[1,0,-s(2);0,1,s(1);0,0,1];
  K=inv(T*R*Ch*R'*T');
end
function [Keff,f1,Meff,info]=flexA(L,t,b,W,H,mcoil,E,rho)
  % 4-힌지 평행판 (homework2 layout), 운동방향=y
  dy=H/2-t/2; sL={[W/2;-dy],[W/2;dy],[-W/2;-dy],[-W/2;dy]}; th=[0,0,pi,pi];
  Kt=zeros(3); Ki=cell(4,1);
  for i=1:4, Ki{i}=hingeK(L,t,b,th(i),sL{i},E); Kt=Kt+Ki{i}; end
  mB=rho*W*b*H; mtot=mB+mcoil; Im=mB*(W^2+H^2)/12;
  M=diag([mtot,mtot,Im]);
  Keff=Kt(2,2); Meff=mtot;
  fn=sort(sqrt(abs(eig(Kt,M)))/(2*pi)); f1=fn(1);
  info.Kt=Kt; info.Ki={Ki}; info.sL={sL}; info.th=th; info.L=L; info.t=t; info.b=b; info.Iz=b*t^3/12;
end
function sig=flexStress(info,x,KT)
  % 운동방향(y) stroke x 에서 힌지별 고정단 굽힘응력 max (homework3 방식) * KT
  u=[0;x;0]; Ki=info.Ki{1}; sL=info.sL{1}; th=info.th; L=info.L; t=info.t; b=info.b; Iz=info.Iz;
  sig=0;
  for i=1:4
    fb=Ki{i}*u; s=sL{i}; Tc=[1,0,0;0,1,0;-s(2),s(1),1]; fh=Tc*fb;
    R=[cos(th(i)),-sin(th(i)),0;sin(th(i)),cos(th(i)),0;0,0,1]; fl=R'*fh;
    Mf=abs(fl(2)*L+fl(3)); sg=Mf*(t/2)/Iz+abs(fl(1))/(b*t);
    if sg>sig,sig=sg; end
  end
  sig=KT*sig;
end
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3; wc=tgM-2*c_c; dcs=rp*dc; hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs)); L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n=nsym*n_half; n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half,n);
  v.Kf=n_eff*Bg*L_m; v.mcoil=8960*n*(pi/4)*dc^2*L_m;
  A_g=2*pi*(rmM+tgM/2)*hyM; v.By1=Bg*A_g/(pi*rmM^2); v.VCMd=2*(rmM+tgM+ty3); v.wc=wc; v.dcs=dcs;
end
function r=perfCD(x,m,dc,P)
  L=x(1)/1e3;t=x(2)/1e3;b=x(3)/1e3;a=x(4)/1e3; v=vcm(x(5),x(6),x(7),dc,P);
  [Keff,f1,~,info]=flexA(L,t,b,a,a,v.mcoil,m.E,m.rho);
  xs=v.Kf*P.I/Keff;
  r.Keff=Keff;r.f1=f1;r.x=xs;r.Kf=v.Kf;r.mcoil=v.mcoil;r.By1=v.By1;r.VCMd=v.VCMd;
  r.dcs=v.dcs;r.wc=v.wc;r.info=info;r.foot=a+2*L;
end
function [c,ceq]=conCD(x,m,dc,P)
  r=perfCD(x,m,dc,P); g=x(8); J=P.I/((pi/4)*dc^2)/1e6;
  sig=flexStress(r.info, g*1e-3, P.KT);
  c=[g*1e-3-r.x; g*100-r.f1; sig-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2; J-P.Jmax];
  ceq=[];
end
function [bx,r]=solveCD(m,dc,P,opts,ns)
  lb=[10,0.2,5,15,5,2,1.2,0.1]; ub=[70,1.5,20,50,8,7,2.5,3]; x0=[30,0.4,10,25,7,4,1.5,0.6];
  bg=-Inf;bx=x0; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conCD(x,m,dc,P),opts);
      if ef>0,c=conCD(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;bx=xo;end,end
    catch
    end
  end
  r=perfCD(bx,m,dc,P);
end
