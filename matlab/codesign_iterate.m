%% codesign_iterate.m - 반복(block-coordinate) 공동최적화 + 일괄과 교차검증
%  트랙: leaf + Mg + d_c=0.3  (일괄 결과 gamma=1.002 와 비교)
%  [A] 유연기구(t,L,b,a) 최적화 | Kf,mcoil 고정  -> Keff,Mflex
%  [B] VCM(rm,hy,tg) 최적화      | Keff,Mflex 고정 -> Kf,mcoil
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; dc=0.3e-3;
m.E=45e9; m.Sy=200e6; m.rho=1740;     % Mg
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200,...
     'MaxFunctionEvaluations',2000,'ConstraintTolerance',1e-9);

% 초기 VCM — 약한/가벼운 쪽에서 출발 (basin 테스트)
Kf=1.3; mcoil=0.8e-3;
fprintf('=== 반복 공동최적화 (leaf+Mg+d_c0.3) ===\n');
fprintf('%-5s%9s%9s%10s%9s%9s\n','iter','gamma','Kf[N/A]','mcoil[g]','Keff','f1[Hz]');
fprintf('%s\n',repmat('-',1,52));
gprev=-1;
for it=1:8
  % ---- [A] 유연기구: vars [t,L,b,a,g], Kf/mcoil 고정 ----
  lbA=[0.3,15,5,15,0.1]; ubA=[2.0,60,20,45,3]; x0A=[0.6,35,10,25,0.7];
  [xA]=solveMS(@(x)conA(x,m,Kf,mcoil,P),lbA,ubA,x0A,opts,12);
  rA=flexL(xA(1),xA(2),xA(3),xA(4),m); Keff=rA.Keff; Mflex=rA.Mflex;
  % ---- [B] VCM: vars [rm,hy,tg,g], Keff/Mflex 고정 ----
  lbB=[5,2,1.2,0.1]; ubB=[8,7,2.5,3]; x0B=[7,4,1.5,0.7];
  [xB]=solveMS(@(x)conB(x,Keff,Mflex,dc,P),lbB,ubB,x0B,opts,12);
  v=vcm(xB(1),xB(2),xB(3),dc,P); Kf=v.Kf; mcoil=v.mcoil;
  % ---- 시스템 gamma ----
  xs=Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/(Mflex+mcoil));
  g=min(xs/1e-3, f1/100);
  fprintf('%-5d%9.4f%9.2f%10.2f%9.1f%9.1f\n', it, g, Kf, mcoil*1e3, Keff/1e3, f1);
  if abs(g-gprev)<1e-3, fprintf('  -> 수렴 (Δγ<1e-3)\n'); break; end
  gprev=g;
end
fprintf('\n일괄(monolithic) 결과: gamma=1.002.  반복 수렴값과 비교 -> 일치하면 교차검증 OK\n');

%% ===== solver (max last var) =====
function bx=solveMS(cf,lb,ub,x0,opts,ns)
  n=numel(lb); bg=-Inf; bx=x0; S=[x0; lb+rand(ns,n).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(end),S(s,:),[],[],[],[],lb,ub,cf,opts);
      if ef>0, c=cf(xo); if all(c<=1e-6)&&(-fv)>bg, bg=-fv; bx=xo; end, end
    catch
    end
  end
end

%% ===== leaf flexure =====
function r=flexL(tm,Lm,bm,am,m)
  t=tm/1e3;L=Lm/1e3;b=bm/1e3;a=am/1e3; db=40e-3;hi=5e-3;
  r.Keff=2*m.E*b*t^3/L^3;
  r.Mflex=m.rho*a^2*b + 0.25*(2*m.rho*db*hi*b) + (1/3)*(8*m.rho*L*t*b);
  r.t=t;r.L=L;r.foot=a+2*(2*L+hi);
end
function [c,ceq]=conA(x,m,Kf,mcoil,P)   % flexure block, Kf/mcoil fixed
  r=flexL(x(1),x(2),x(3),x(4),m); g=x(5);
  xs=Kf*P.I/r.Keff; f1=(1/(2*pi))*sqrt(r.Keff/(r.Mflex+mcoil));
  peak=P.KT*3*m.E*r.t*((g*1e-3)/2)/r.L^2;
  c=[g*1e-3-xs; g*100-f1; peak-P.SIG*m.Sy; r.foot-200e-3];
  ceq=[];
end

%% ===== VCM (N52 surrogate) =====
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3; wc=tgM-2*c_c; dcs=rp*dc; hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs)); L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half, nsym*n_half);
  v.Kf=n_eff*Bg*L_m; v.mcoil=8960*(nsym*n_half)*(pi/4)*dc^2*L_m;
  A_g=2*pi*(rmM+tgM/2)*hyM; v.By1=Bg*A_g/(pi*rmM^2); v.VCMd=2*(rmM+tgM+ty3); v.wc=wc; v.dcs=dcs;
end
function [c,ceq]=conB(x,Keff,Mflex,dc,P)   % VCM block, Keff/Mflex fixed
  v=vcm(x(1),x(2),x(3),dc,P); g=x(4);
  xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/(Mflex+v.mcoil));
  c=[g*1e-3-xs; g*100-f1; v.VCMd-20e-3; v.dcs-v.wc; v.By1-2.2];
  ceq=[];
end
