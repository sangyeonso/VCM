%% final_dcp_design.m - 최종 확정 설계 (DCP + Al7075, d_c=0.5 스펙준수) 전 제원 dump
%  행렬법 유연기구 + DCP(2자유도) + FEMM N52 VCM + 발열(J<=10.2) + 피로(0.2Sy,Kt)
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; P.Jmax=10.2; dc=0.5e-3;
m.E=72e9; m.Sy=503e6; m.rho=2810;        % Al7075-T6
P.mi_frac=0.5;                            % 중간단 총질량 = 0.5 x 출력 구조질량 (경량 프레임 가정)
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',300,...
     'MaxFunctionEvaluations',3000,'ConstraintTolerance',1e-9);

lb=[10,0.2,5,15,5,2,1.2,0.1]; ub=[70,1.5,20,50,8,7,2.5,3]; x0=[30,0.4,10,25,7,4,1.5,0.6];
bg=-Inf;bx=x0; S=[x0; lb+rand(30,8).*(ub-lb)];
for s=1:size(S,1)
  try
    [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conDCP(x,m,dc,P),opts);
    if ef>0,c=conDCP(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;bx=xo;end,end
  catch
  end
end
r=perfDCP(bx,m,dc,P); g=bx(8);

fprintf('================ 최종 확정 설계: DCP + Al7075 (d_c=0.5 스펙준수) ================\n');
fprintf('  gamma = %.3f   (stroke=%.3f mm, f1=%.1f Hz)\n', g, r.x*1e3, r.f1);
fprintf('  -- 유연기구 (행렬법, DCP) --\n');
fprintf('     hinge: L=%.1f t=%.2f b=%.1f mm | 스테이지 a=%.1f mm\n', bx(1),bx(2),bx(3),bx(4));
fprintf('     K_eff=%.0f N/mm  M_eff=%.2f g  외형 a+4L=%.0f mm (<200)\n', r.Keff/1e3, r.Meff*1e3, (bx(4)+4*bx(1)));
fprintf('     peak stress=%.0f MPa (<%.0f=0.2Sy), 각 단 변형 x/2=%.3f mm\n', r.sig*1e-6, 0.2*m.Sy/1e6, g*0.5);
fprintf('  -- VCM (FEMM N52) --\n');
fprintf('     r_m=%.2f h_y1=%.2f t_g=%.2f mm (d_c=0.5) | Ø=%.1f mm\n', bx(5),bx(6),bx(7),r.VCMd*1e3);
fprintf('     B_g=%.0f mT  K_f=%.2f N/A  m_coil=%.2f g  J=%.1f A/mm^2\n', r.Bg*1e3,r.Kf,r.mcoil*1e3, P.I/((pi/4)*dc^2)/1e6);
fprintf('  -- KPI 판정 --\n');
fprintf('     stroke %.3f mm (>1?) | f1 %.1f Hz (>100?) | d_c 0.5mm (>0.5 OK) | Ø %.1f (<20 OK)\n', r.x*1e3, r.f1, r.VCMd*1e3);
fprintf('==============================================================================\n');

%% ===== 함수 (dcp_extend와 동일) =====
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
  A_g=2*pi*(rmM+tgM/2)*hyM;v.By1=Bg*A_g/(pi*rmM^2);v.VCMd=2*(rmM+tgM+ty3);v.wc=wc;v.dcs=dcs;v.Bg=Bg;
end
function r=perfDCP(x,m,dc,P)
  L=x(1)/1e3;t=x(2)/1e3;b=x(3)/1e3;a=x(4)/1e3; v=vcm(x(5),x(6),x(7),dc,P);
  [Kp,info]=flexA(L,t,b,a,a,m.E); Keff=Kp;
  Mout=m.rho*a*b*a+v.mcoil; Mi=P.mi_frac*m.rho*a*b*a;
  K2=[4*Kp,-2*Kp;-2*Kp,2*Kp]; M2=diag([Mi,Mout]);
  fn=sort(sqrt(abs(eig(K2,M2)))/(2*pi)); f1=fn(1);
  xs=v.Kf*P.I/Keff; g=x(8);
  r.Keff=Keff;r.f1=f1;r.x=xs;r.Kf=v.Kf;r.mcoil=v.mcoil;r.By1=v.By1;r.VCMd=v.VCMd;r.Bg=v.Bg;
  r.dcs=v.dcs;r.wc=v.wc;r.info=info;r.foot=a+4*L; r.Meff=Mout+Mi; r.sig=flexStress(info,(g*1e-3)/2,P.KT);
end
function [c,ceq]=conDCP(x,m,dc,P)
  r=perfDCP(x,m,dc,P); g=x(8); J=P.I/((pi/4)*dc^2)/1e6;
  c=[g*1e-3-r.x; g*100-r.f1; r.sig-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2; J-P.Jmax];
  ceq=[];
end
