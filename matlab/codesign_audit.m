%% codesign_audit.m - 전수 감사: KPI 독립 재계산 + 결측 물리(발열/전류밀도) 점검
%  (A) 권장설계 2종(Al7075 d_c0.3, Mg d_c0.5) 전 항목 dump + 발열지표
%  (B) 전류밀도 J<=10A/mm^2 제약 추가 -> d_c 자유화 -> 진짜 최적 d_c
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; P.rho_e=1.68e-8;   % 구리 전기저항률
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);

%% (A) 두 설계 전 항목 dump + 발열
fprintf('=== (A) KPI 독립 재계산 + 발열지표 ===\n');
cases={'Al7075 d_c0.3',72e9,503e6,2810,0.3e-3; 'Mg d_c0.5',45e9,200e6,1740,0.5e-3};
for c=1:2
  m.E=cases{c,2};m.Sy=cases{c,3};m.rho=cases{c,4}; dc=cases{c,5};
  [bx,r]=solveleaf(m,dc,P,opts,16);
  J=P.I/((pi/4)*dc^2)/1e6;       % A/mm^2
  fprintf('\n[%s]  gamma=%.3f\n', cases{c,1}, bx(8));
  fprintf('  유연: t=%.3f L=%.1f b=%.1f a=%.1f mm | K_eff=%.1f N/mm\n',bx(1),bx(2),bx(3),bx(4),r.Keff/1e3);
  fprintf('  VCM : r_m=%.2f h_y1=%.2f t_g=%.2f mm | B_g=%.0fmT K_f=%.2f N/A n=%.0f turns\n',...
     bx(5),bx(6),bx(7),r.Bg*1e3,r.Kf,r.n);
  fprintf('  질량: M_flex=%.2fg m_coil=%.2fg M_eff=%.2fg\n', (r.Meff-r.mcoil)*1e3, r.mcoil*1e3, r.Meff*1e3);
  fprintf('  KPI : stroke=%.3fmm (>1?) f1=%.1fHz (>100?) peak=%.0fMPa (<%.0f=0.2Sy?)\n',...
     r.x*1e3, r.f1, P.KT*3*m.E*r.t*((bx(8)*1e-3)/2)/r.L^2/1e6, 0.2*m.Sy/1e6);
  fprintf('  >> 발열: 전류밀도 J=%.1f A/mm^2, R=%.2f ohm, V=I*R=%.1fV, P=I^2R=%.2fW\n',...
     J, r.R, P.I*r.R, P.I^2*r.R);
end

%% (B) J<=10 제약 추가, d_c 자유화 -> 진짜 최적 d_c (Al7075)
fprintf('\n=== (B) 전류밀도 J<=10 A/mm^2 제약 추가, d_c 자유화 (Al7075) ===\n');
m.E=72e9;m.Sy=503e6;m.rho=2810;
for Jmax=[Inf 10 7]
  % vars [t,L,b,a, rm,hy,tg, dc, g]
  lb=[0.3,15,5,15,5,2,1.2,0.2,0.1]; ub=[2,60,20,45,8,7,2.5,1.0,3]; x0=[0.6,35,10,25,7,4,1.5,0.4,0.7];
  bg=-Inf;bx=x0; S=[x0; lb+rand(24,9).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(9),S(s,:),[],[],[],[],lb,ub,@(x)conJ(x,m,P,Jmax),opts);
      if ef>0,c=conJ(xo,m,P,Jmax); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;bx=xo;end,end
    catch
    end
  end
  Js=P.I/((pi/4)*(bx(8)*1e-3)^2)/1e6;
  jt=' (무시)'; if isfinite(Jmax), jt=sprintf(' (<=%g)',Jmax); end
  fprintf('  J_max=%-6s%s : gamma=%.3f, 최적 d_c=%.3fmm, J=%.1f A/mm^2\n',...
     num2str(Jmax), jt, bx(9), bx(8), Js);
end
fprintf('\n해석: J 무시하면 d_c가 가늘게 몰림(유리). J 제약 넣으면 최적 d_c가 굵어짐.\n');
fprintf('      -> d_c>0.5 스펙이 사실 "발열 한계"의 proxy인지 이 결과로 판단.\n');

%% ===== 함수 =====
function [bx,r]=solveleaf(m,dc,P,opts,ns)
  lb=[0.3,15,5,15,5,2,1.2,0.1]; ub=[2,60,20,45,8,7,2.5,3]; x0=[0.6,35,10,25,7,4,1.5,0.7];
  bg=-Inf;bx=x0; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conL(x,m,dc,P),opts);
      if ef>0,c=conL(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg,bg=-fv;bx=xo;end,end
    catch
    end
  end
  r=perfL(bx(1),bx(2),bx(3),bx(4),bx(5),bx(6),bx(7),m.E,m.rho,dc,P);
end
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3; wc=tgM-2*c_c; dcs=rp*dc; hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs)); L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n=nsym*n_half; n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half, n);
  v.Kf=n_eff*Bg*L_m; v.mcoil=8960*n*(pi/4)*dc^2*L_m; v.n=n;
  v.R=P.rho_e*n*L_m/((pi/4)*dc^2);            % 코일 저항
  A_g=2*pi*(rmM+tgM/2)*hyM; v.By1=Bg*A_g/(pi*rmM^2); v.VCMd=2*(rmM+tgM+ty3); v.wc=wc; v.dcs=dcs; v.Bg=Bg;
end
function r=perfL(tm,Lm,bm,am,rm,hy,tg,E,rho,dc,P)
  t=tm/1e3;L=Lm/1e3;b=bm/1e3;a=am/1e3; v=vcm(rm,hy,tg,dc,P); db=40e-3;hi=5e-3;
  Keff=2*E*b*t^3/L^3; Mflex=rho*a^2*b+0.25*(2*rho*db*hi*b)+(1/3)*(8*rho*L*t*b);
  Meff=Mflex+v.mcoil; xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.Keff=Keff;r.Meff=Meff;r.x=xs;r.f1=f1;r.Kf=v.Kf;r.mcoil=v.mcoil;r.n=v.n;r.R=v.R;r.Bg=v.Bg;
  r.By1=v.By1;r.VCMd=v.VCMd;r.dcs=v.dcs;r.wc=v.wc;r.t=t;r.L=L;r.foot=a+2*(2*L+hi);
end
function [c,ceq]=conL(x,m,dc,P)
  r=perfL(x(1),x(2),x(3),x(4),x(5),x(6),x(7),m.E,m.rho,dc,P); g=x(8); xop=g*1e-3;
  peak=P.KT*3*m.E*r.t*(xop/2)/r.L^2;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end
function [c,ceq]=conJ(x,m,P,Jmax)   % d_c 자유 + J 제약
  dc=x(8)*1e-3;
  r=perfL(x(1),x(2),x(3),x(4),x(5),x(6),x(7),m.E,m.rho,dc,P); g=x(9); xop=g*1e-3;
  peak=P.KT*3*m.E*r.t*(xop/2)/r.L^2; J=P.I/((pi/4)*dc^2)/1e6;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  if isfinite(Jmax), c=[c; J-Jmax]; end
  ceq=[];
end
