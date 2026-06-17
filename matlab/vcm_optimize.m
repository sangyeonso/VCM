%% vcm_optimize.m
%  VCM 모터 자체 최적화 (surrogate 아님 — 완전 해석 퍼미언스 모델, 전체 형상 6변수)
%  목적: K_f/m_coil 최대화  (= d_c 고정시 B_g 최대화, n_eff/n=1 클리핑)
%  변수 x=[r_m, h_m, h_y1, h_y2, t_g, t_y3] (mm)
%  고정: d_c=0.5mm, N52(Br=1.43), FeCo 요크(mu_s 높음, By1<=2.2T 포화)
%  제약: Ø=2(r_m+t_g+t_y3)<=20, 코일끼움 w_c>=d_c*, 포화 By1<=2.2
%  ※ 해석모델은 절대 B_g를 ~1.5x 과대평가 -> 최적 '형상' 탐색용, 절대값은 FEMM 검증
clear; clc; rng(3);
mu0=4*pi*1e-7;

P.Br=1.43; P.Hc=1.43/(mu0*1.07); P.mu_s=5000; P.dc=0.5e-3; P.rp=1.15;
P.cc=0.3e-3; P.tb=0.3e-3; P.nsym=2; P.I=2; P.By1max=2.2;

% 변수 경계 [r_m h_m h_y1 h_y2 t_g t_y3] (mm)
lb=[2, 10, 2, 1, 1.2, 0.5];
ub=[8, 100,16, 8, 4.0, 3.0];
x0=[6, 40, 5, 3, 2.0, 1.0];
opts=optimoptions('fmincon','Display','off','Algorithm','sqp',...
     'MaxIterations',400,'MaxFunctionEvaluations',4000,'ConstraintTolerance',1e-9);

% 다중시작 최적화 (B_g 최대화)
best=-Inf; bx=x0;
starts=[x0; lb+rand(40,6).*(ub-lb)];
for s=1:size(starts,1)
  try
    [xo,fv,ef]=fmincon(@(x)-vcm(x,P,mu0).Bg, starts(s,:),[],[],[],[],lb,ub,@(x)con(x,P,mu0),opts);
    if ef>0
      c=con(xo,P,mu0);
      if all(c<=1e-6) && (-fv)>best, best=-fv; bx=xo; end
    end
  catch
  end
end
r=vcm(bx,P,mu0);

% --- 기준선(현 설계) 비교 ---
xb=[6.5, 25, 5.8, 3, 2.5, 1.0];  rb=vcm(xb,P,mu0);

fprintf('=== VCM 모터 최적화 (완전 해석모델, 전체 6변수, max K_f/m_coil) ===\n\n');
fprintf('%-18s%12s%12s\n','항목','기준설계','최적설계');
fprintf('%s\n',repmat('-',1,42));
fprintf('%-18s%12s%12s\n','--- 형상[mm] ---','','');
nm={'r_m','h_m','h_y1','h_y2','t_g','t_y3'};
for i=1:6, fprintf('%-18s%12.2f%12.2f\n', nm{i}, xb(i), bx(i)); end
fprintf('%-18s%12s%12s\n','--- 성능 ---','','');
fprintf('%-18s%12.0f%12.0f\n','B_g [mT] (해석)', rb.Bg*1e3, r.Bg*1e3);
fprintf('%-18s%12.2f%12.2f\n','K_f [N/A]', rb.Kf, r.Kf);
fprintf('%-18s%12.2f%12.2f\n','m_coil [g]', rb.mc*1e3, r.mc*1e3);
fprintf('%-18s%12.0f%12.0f\n','K_f/m_coil [N/A/kg... ]', rb.Kf/rb.mc, r.Kf/r.mc);
fprintf('%-18s%12.2f%12.2f\n','By1 [T] (<=2.2)', rb.By1, r.By1);
fprintf('%-18s%12.1f%12.1f\n','VCM Ø [mm]', rb.D*1e3, r.D*1e3);
fprintf('%-18s%12.1f%12.1f\n','추력 @2A [N]', rb.Kf*P.I, r.Kf*P.I);
fprintf('\nK_f/m_coil 개선: %.0f -> %.0f (%.1f%% UP)\n', rb.Kf/rb.mc, r.Kf/r.mc, ...
        100*(r.Kf/r.mc - rb.Kf/rb.mc)/(rb.Kf/rb.mc));
fprintf('-> 최적 형상: 검증용 FEMM 파라미터로 사용 (다음 단계)\n');

%% ===== 함수 =====
function r=vcm(x,P,mu0)
  r_m=x(1)*1e-3; h_m=x(2)*1e-3; h_y1=x(3)*1e-3; h_y2=x(4)*1e-3; t_g=x(5)*1e-3; t_y3=x(6)*1e-3;
  mu_r=P.Br/(mu0*P.Hc);
  P_m =pi*mu0*mu_r*r_m^2/h_m;
  P_y1=pi*mu0*P.mu_s*r_m^2/h_y1;
  P_y2=2*pi*mu0*P.mu_s*h_y2/(log(r_m+t_g+t_y3/2)-log(r_m/2));
  P_y3=pi*mu0*P.mu_s*((r_m+t_g+t_y3)^2-(r_m+t_g)^2)/(h_y1+h_m);
  P_g =2*pi*mu0*h_y1/(log(r_m+t_g)-log(r_m));
  P_l =0.52*(2*pi*mu0*r_m);
  Phi_r=P.Br*pi*r_m^2;
  P_main=1/(1/P_g+1/P_y1+1/P_y3+1/P_y2);
  Phi_g=Phi_r*P_main/(P_m+P_l+P_main);
  A_g=2*pi*(r_m+t_g/2)*h_y1; Bg=Phi_g/A_g; By1=Phi_g/(pi*r_m^2);
  % 코일
  dcs=P.rp*P.dc; w_c=t_g-2*P.cc; h_c=h_y1;
  n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs));
  L_m=2*pi*(r_m+P.cc+P.tb+w_c/2);
  n_eff=P.nsym*(h_y1+0.08*h_m)/h_c*n_half; n_eff=min(n_eff,P.nsym*n_half);
  Kf=n_eff*Bg*L_m; mc=8960*(P.nsym*n_half)*(pi/4)*P.dc^2*L_m;
  r.Bg=Bg; r.By1=By1; r.Kf=Kf; r.mc=mc; r.D=2*(r_m+t_g+t_y3); r.wc=w_c; r.dcs=dcs;
end

function [c,ceq]=con(x,P,mu0)
  r=vcm(x,P,mu0);
  c(1)=r.D-20e-3;          % VCM 직경
  c(2)=r.dcs-r.wc;         % 코일 끼움
  c(3)=r.By1-P.By1max;     % 요크 포화
  ceq=[];
end
