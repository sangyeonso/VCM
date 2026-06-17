%% fea_frame_modal.m
%  유연기구 FEA 검증 — 2D 프레임 유한요소(Euler-Bernoulli 빔) 모델 (툴박스 불필요)
%  목적: 해석모델의 (1) 강성 k_x=24EI/L^3 와 (2) 공진 f1(유효질량 1/4,1/3 가정)을
%        분산질량 consistent-mass 고유치 해석으로 교차검증.
%  Double Compound Parallelogram: blade 8개(thin) + bar 3개(stage/inter, thick),
%  상하 ground 고정. 운동축 X.
clear; clc;

%% ---- 형상/물성 (Al6061, d_c=0.3 최적해 근방) ----
E=71e9; rho=2770; nu=0.33;
t=0.56e-3;  b=8e-3;  L=29e-3;        % blade 두께/면외깊이/길이
xb=15e-3;                            % blade 두 열의 x=±xb
h_stage=12e-3; h_inter=6e-3;         % stage/inter bar 높이(Y)
nseg=4;                              % 멤버당 요소 분할

% blade 단면: 운동방향(X) 굽힘 I=b*t^3/12, A=t*b
Ab=t*b;          Ib=b*t^3/12;
% bar 단면(축=X): A=h*b, I=b*h^3/12
As=h_stage*b;    Is=b*h_stage^3/12;
Ai=h_inter*b;    Aii=b*h_inter^3/12;

%% ---- 멤버 정의: [x1 y1 x2 y2 A I] ----
yL=[0 L 2*L 3*L 4*L];
M=[];
% blades (수직) : 양 열(x=±xb), 4구간
for sx=[-xb xb]
  for q=1:4
    M=[M; sx yL(q) sx yL(q+1) Ab Ib];
  end
end
% bars (수평): y=L(inter), 2L(stage), 3L(inter)
M=[M; -xb L  xb L  Ai Aii];
M=[M; -xb 2*L xb 2*L As Is];
M=[M; -xb 3*L xb 3*L Ai Aii];

%% ---- 노드 생성(분할) & 요소 조립 ----
keyf=@(x,y)sprintf('%.6f_%.6f',round(x,9),round(y,9));
nodeMap=containers.Map('KeyType','char','ValueType','double');
X=[]; elems=[];
getNode=@(x,y)0;
for i=1:size(M,1)
    x1=M(i,1);y1=M(i,2);x2=M(i,3);y2=M(i,4);A=M(i,5);I=M(i,6);
    xs=linspace(x1,x2,nseg+1); ys=linspace(y1,y2,nseg+1);
    ids=zeros(1,nseg+1);
    for k=1:nseg+1
        kk=keyf(xs(k),ys(k));
        if ~isKey(nodeMap,kk)
            X=[X; xs(k) ys(k)]; nodeMap(kk)=size(X,1);
        end
        ids(k)=nodeMap(kk);
    end
    for k=1:nseg
        elems=[elems; ids(k) ids(k+1) A I];
    end
end
nN=size(X,1); nDOF=3*nN;
K=zeros(nDOF); Mass=zeros(nDOF);

for e=1:size(elems,1)
    n1=elems(e,1);n2=elems(e,2);A=elems(e,3);I=elems(e,4);
    x1=X(n1,1);y1=X(n1,2);x2=X(n2,1);y2=X(n2,2);
    Le=hypot(x2-x1,y2-y1); c=(x2-x1)/Le; s=(y2-y1)/Le;
    [ke,me]=beam2d(E,A,I,rho,Le);
    T=blkdiag([c s 0;-s c 0;0 0 1],[c s 0;-s c 0;0 0 1]);
    keg=T'*ke*T; meg=T'*me*T;
    d=[3*n1-2 3*n1-1 3*n1 3*n2-2 3*n2-1 3*n2];
    K(d,d)=K(d,d)+keg; Mass(d,d)=Mass(d,d)+meg;
end

%% ---- 경계조건: 상하 ground(±xb, y=0 및 4L) 고정 ----
fixed=[];
for p=[ -xb 0; xb 0; -xb 4*L; xb 4*L ]'
    kk=keyf(p(1),p(2)); nid=nodeMap(kk);
    fixed=[fixed 3*nid-2 3*nid-1 3*nid];
end
free=setdiff(1:nDOF,fixed);

%% ---- 모달 해석 ----
[V,D]=eig(K(free,free),Mass(free,free));
w=sqrt(abs(diag(D))); f=sort(w)/(2*pi);

fprintf('=== 유연기구 FEA (2D 프레임 빔요소, %d 노드, %d 요소) ===\n',nN,size(elems,1));
fprintf('형상: t=%.2fmm L=%.0fmm b=%.0fmm, stage bar h=%.0fmm, inter h=%.0fmm\n\n',...
    t*1e3,L*1e3,b*1e3,h_stage*1e3,h_inter*1e3);
fprintf('FEA 하위 4개 모드 [Hz]: %.1f, %.1f, %.1f, %.1f\n', f(1),f(2),f(3),f(4));

%% ---- 해석모델 비교 ----
k_x=24*E*Ib/L^3;
m_st = rho*(2*xb)*h_stage*b;       % stage bar
m_in = rho*(2*xb)*h_inter*b;       % inter bar (1개)
m_bl = rho*L*t*b;                  % blade 1개
m_eff= m_st + 0.25*(2*m_in) + (1/3)*(8*m_bl);
f_an = (1/(2*pi))*sqrt(k_x/m_eff);

fprintf('\n--- 해석모델 ---\n');
fprintf('k_x(폐형식 24EI/L^3) = %.0f N/m\n',k_x);
fprintf('m_eff = stage %.2fg + 1/4*inter %.2fg + 1/3*blade %.2fg = %.2f g\n',...
    m_st*1e3, 0.25*2*m_in*1e3, (1/3)*8*m_bl*1e3, m_eff*1e3);
fprintf('해석 f1 = %.1f Hz\n', f_an);
fprintf('\n>>> FEA f1 = %.1f Hz vs 해석 f1 = %.1f Hz  (오차 %.1f%%)\n',...
    f(1), f_an, 100*(f(1)-f_an)/f_an);

%% ---- 1차 모드형상 그림 저장 (보고서용) ----
[~,ord]=sort(w);
phi=zeros(nDOF,1); phi(free)=V(:,ord(1));
ux=phi(1:3:end); uy=phi(2:3:end);
sc=8e-3/max(abs([ux;uy]));
figure('Position',[100 100 420 470]); hold on; axis equal off;
for e=1:size(elems,1)
    n1=elems(e,1); n2=elems(e,2);
    plot([X(n1,1) X(n2,1)]*1e3,[X(n1,2) X(n2,2)]*1e3,'-','Color',[0.75 0.75 0.75],'LineWidth',1);
    plot(([X(n1,1) X(n2,1)]+sc*[ux(n1) ux(n2)])*1e3,([X(n1,2) X(n2,2)]+sc*[uy(n1) uy(n2)])*1e3,'-b','LineWidth',1.8);
end
plot(nan,nan,'-','Color',[0.75 0.75 0.75]); plot(nan,nan,'-b');
legend({'undeformed','1st mode (scaled)'},'Location','southoutside','Orientation','horizontal');
title(sprintf('FEA 1st mode = %.1f Hz  (X-translation, DCP)',f(1)));
exportgraphics(gcf,fullfile('..','figure','fig7_fea_modeshape.png'),'Resolution',200);

%% ---- 빔 요소 행렬 ----
function [ke,me]=beam2d(E,A,I,rho,L)
    EA=E*A; EI=E*I;
    ke=[ EA/L,0,0,-EA/L,0,0;
         0,12*EI/L^3,6*EI/L^2,0,-12*EI/L^3,6*EI/L^2;
         0,6*EI/L^2,4*EI/L,0,-6*EI/L^2,2*EI/L;
         -EA/L,0,0,EA/L,0,0;
         0,-12*EI/L^3,-6*EI/L^2,0,12*EI/L^3,-6*EI/L^2;
         0,6*EI/L^2,2*EI/L,0,-6*EI/L^2,4*EI/L];
    m=rho*A*L;
    me=(m/420)*[140,0,0,70,0,0;
        0,156,22*L,0,54,-13*L;
        0,22*L,4*L^2,0,13*L,-3*L^2;
        70,0,0,140,0,0;
        0,54,13*L,0,156,-22*L;
        0,-13*L,-3*L^2,0,-22*L,4*L^2];
end
