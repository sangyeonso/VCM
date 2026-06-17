%% flexure_notch_lagrange.m
%  강의(제07강) 방법으로 유연기구 재모델링: 노치힌지(Paros-Weisbord) + 라그랑지
%   k_th=2Eb t^2.5/(9pi sqrt(R)),  K_eff=4k_th/l^2
%   M_eff=M+1/2 m+1/2 m_l+2J/l^2  (라그랑지)
%   theta_max=(4K/Kt)(R/Et)sig_max,  x_max=2 l theta_max  (응력한계 travel, Kt 포함)
%  (A) 강의 예제 재현 검증   (B) 우리 스펙(1mm,100Hz) 적용/최적화
clear; clc; rng(5);

%% ===== (A) 강의 예제 재현 (Al6061, R2 b12 t0.6 l15, sig=0.1Sy) =====
matAl.E=72e9; matAl.Sy=270e6; matAl.rho=2770;
v=notchDCP(matAl,2,12,0.6,15,40,0.1,0,3,3);
fprintf('=== (A) 강의 예제 재현 검증 ===\n');
fprintf('K=%.4f (강의 0.3355), Kt=%.4f (1.0732)\n', v.K, v.Kt);
fprintf('theta_max=%.4f rad (0.0016), x_max=%.1f um (46.9)\n', v.thmax, v.xmax*1e6);
fprintf('K_eff=%.0f N/mm (210), M_eff=%.4f kg (~0.06), f1=%.0f Hz (300)\n\n', v.Keff/1e3, v.Meff, v.f1);

%% ===== (B) 우리 스펙 적용: x_max>=1mm & f1>=100Hz 동시충족 가능? =====
% 재료 3종, gamma 최대화: x_max>=g*1mm, f1>=g*100Hz
mats={'Al6061',72e9,270e6,2770; 'Spring-steel',200e9,1200e6,7850; 'BeCu',128e9,1100e6,8250};
SIG=0.3;        % sigma_max=0.3*Sy (공격적; Kt는 이미 포함). 피로면 더 낮춤
MCOIL=4e-3;     % VCM 가동코일 [kg]
KF=3.39; II=2;  % 최적 VCM(d_c=0.3) 힘상수[N/A], 전류[A] -> F=6.78N
%      [ R    t    b    l    a_plat  gamma]
lb=[0.3, 0.2, 5,   10,  15,  0.1];
ub=[5.0, 1.5, 20,  75,  50,  3.0];
x0=[2.0, 0.5, 12,  40,  30,  0.5];
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',300,...
     'MaxFunctionEvaluations',3000,'ConstraintTolerance',1e-9);

fprintf('=== (B) 노치힌지로 가능? VCM 구동력 제약 포함 (x_force=K_f I/K_eff>=1mm) ===\n');
fprintf('%-13s%8s%10s%11s%9s%8s%9s\n','재료','gamma','xmax[mm]','x_force[mm]','f1[Hz]','l[mm]','Keff[N/mm]');
fprintf('%s\n',repmat('-',1,72));
for k=1:3
  m.E=mats{k,2}; m.Sy=mats{k,3}; m.rho=mats{k,4};
  best=-Inf; bx=x0;
  starts=[x0; lb+rand(20,6).*(ub-lb)];
  for s=1:size(starts,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(6),starts(s,:),[],[],[],[],lb,ub,@(x)con(x,m,SIG,MCOIL,KF,II),opts);
      if ef>0, c=con(xo,m,SIG,MCOIL,KF,II); if all(c<=1e-6)&&(-fv)>best, best=-fv; bx=xo; end, end
    catch
    end
  end
  r=notchDCP(m,bx(1),bx(3),bx(2),bx(4),bx(5),SIG,MCOIL,3,3);
  xforce=KF*II/r.Keff;
  vd=''; if bx(6)>=1.0, vd=' <= 충족'; end
  fprintf('%-13s%8.3f%10.3f%11.3f%9.1f%8.1f%9.0f%s\n',...
     mats{k,1}, bx(6), r.xmax*1e3, xforce*1e3, r.f1, bx(4), r.Keff/1e3, vd);
end
fprintf('\n주: 이제 gamma는 x_max(응력)·x_force(VCM구동)·f1 셋 다 >= g*목표. x_force가 진짜 병목.\n');

%% ===== 함수 =====
function r=notchDCP(mat,R,b,t,l,a,sigfac,mx,wb,wl)
  Rm=R/1e3; bm=b/1e3; tm=t/1e3; lm=l/1e3; am=a/1e3; wbm=wb/1e3; wlm=wl/1e3;
  E=mat.E; Sy=mat.Sy; rho=mat.rho;
  K =0.565*(t/R)+0.166;
  Kt=(2.7*t+5.4*R)/(8*R+t)+0.325;
  k_th=2*E*bm*tm^2.5/(9*pi*sqrt(Rm));          % N*m/rad
  Keff=4*k_th/lm^2;                            % N/m
  thmax=(4*K/Kt)*(Rm/(E*tm))*(sigfac*Sy);      % rad
  xmax=2*lm*thmax;                             % m
  M=am^2*bm*rho + mx;                          % 플랫폼(+코일)
  m_bar=am*wbm*bm*rho;                         % 중간바 1개
  m_lnk=lm*wlm*bm*rho;                         % 링크 1개
  Meff=M + 0.5*m_bar + (2/3)*m_lnk;            % 라그랑지 (2J/l^2 -> m_lnk/6)
  f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.K=K; r.Kt=Kt; r.k_th=k_th; r.Keff=Keff; r.thmax=thmax; r.xmax=xmax;
  r.Meff=Meff; r.f1=f1; r.M=M;
end

function [c,ceq]=con(x,mat,SIG,MCOIL,KF,II)
  R=x(1); t=x(2); b=x(3); l=x(4); a=x(5); g=x(6);
  r=notchDCP(mat,R,b,t,l,a,SIG,MCOIL,3,3);   % (mat,R,b,t,l,a)
  xforce=KF*II/r.Keff;
  c(1)=g*1e-3 - r.xmax;          % 응력한계 travel >= g*1mm
  c(2)=g*100  - r.f1;            % f1>=g*100Hz
  c(3)=g*1e-3 - xforce;          % VCM 구동 변위 >= g*1mm  (병목)
  c(4)=t - R;                    % R>t
  c(5)=R - 5*t;                  % R<5t
  c(6)=(a+2*l) - 180;            % 외형 여유(<200)
  ceq=[];
end
