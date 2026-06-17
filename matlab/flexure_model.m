%% flexure_model.m
%  Double Compound Parallelogram (DCP) 1-DOF 유연기구 설계 모델
%  - 운동방향 강성 k_x = 24*E*I/L^3 (제06강 8-(2)-4, Double Compound Flexure)
%  - 변위 x = F/k_x, 1차 공진 f1 = (1/2pi)*sqrt(k_x/m_eff)
%  - blade 응력 sigma = 3*E*t*(x/2)/L^2  (compound라 단 stage당 행정 절반)
%  - 재료 3종(Al6061 / 스프링강 / BeCu) 비교
%  - 검증: 빔 강성행렬 응축으로 guided blade = 12EI/L^3 확인 + 토폴로지 직병렬 산술
%
%  좌표: X=운동방향, Y=blade 길이방향, Z=면외(판두께, <=20mm)
%  blade 단면2차모멘트 I = b*t^3/12  (t=두께(X), b=높이(Z))
clear; clc;

%% ===== 설계변수 (nominal) =====
t = 0.8e-3;     % blade 두께 [m]   (가공성 >0.3mm)
L = 40e-3;      % blade 길이 [m]
b = 15e-3;      % blade 높이(면외) [m]  (<=20mm)

%% ===== 고정 형상 (질량 산정용) =====
a_s   = 20e-3;  % moving stage 한 변 (X,Y 정사각) [m]
d_b   = 40e-3;  % parallelogram 두 blade 간격(X) [m] = 중간stage 길이
h_i   = 5e-3;   % 중간 stage 두께(Y) [m]
n_blade = 8;    % 총 blade 수 (Double Compound = 4 levels x 2)
m_coil  = 5e-3; % VCM 가동코일 질량 [kg] (임시값; VCM 모델과 결합 예정)
m_pay   = 0;    % payload [kg] (별도 명시 없음)

x_req = 1e-3;   % 요구 최대변위 [m]
f_req = 100;    % 요구 1차 공진 [Hz]

%% ===== 재료 3종 =====
%        이름            E[Pa]    rho[kg/m3]  Sy[Pa]
mats = { 'Al6061-T6',    71e9,    2770,       276e6;
         'Spring-steel', 200e9,   7850,       1200e6;
         'BeCu C17200',  128e9,   8250,       1100e6 };

I = b*t^3/12;   % 단면 2차 모멘트 [m^4]

fprintf('=== 유연기구 모델 (Double Compound Parallelogram) ===\n');
fprintf('설계변수: t=%.2f mm, L=%.1f mm, b=%.1f mm,  I=%.4e m^4\n\n', ...
        t*1e3, L*1e3, b*1e3, I);

%% ===== (검증 1) guided blade 강성 = 12EI/L^3 (빔 강성행렬 응축) =====
%  2D 빔요소 4x4 (절점당 v,theta). 절점1 고정, 절점2 회전구속(guided) -> k = K(3,3)
E0 = 1; % 단위계수 (형상검증용)
Kb = (E0*I/L^3) * [ 12,    6*L,  -12,    6*L;
                    6*L, 4*L^2, -6*L,  2*L^2;
                   -12,   -6*L,   12,   -6*L;
                    6*L, 2*L^2, -6*L,  4*L^2 ];
k_guided_matrix = Kb(3,3);          % v2-v2 (theta2=0 구속)  -> 12EI/L^3
k_guided_closed = 12*E0*I/L^3;
fprintf('[검증1] guided blade 강성  행렬응축=%.4e,  폐형식 12EI/L^3=%.4e,  오차=%.2e\n', ...
        k_guided_matrix, k_guided_closed, abs(k_guided_matrix-k_guided_closed));

%% ===== (검증 2) 토폴로지 직병렬 산술 -> DCP = 24EI/L^3 =====
E_chk = 71e9;
k1   = 12*E_chk*I/L^3;          % 단일 guided blade
k_para = 2*k1;                  % parallelogram (blade 2개 병렬)       = 24EI/L^3
k_comp = 1/(1/k_para + 1/k_para); % compound (parallelogram 2개 직렬)  = 12EI/L^3
k_dcp  = 2*k_comp;             % double compound (compound 2개 병렬)  = 24EI/L^3
k_dcp_closed = 24*E_chk*I/L^3;
fprintf('[검증2] DCP 직병렬=%.4e,  폐형식 24EI/L^3=%.4e,  오차=%.2e\n\n', ...
        k_dcp, k_dcp_closed, abs(k_dcp-k_dcp_closed));

%% ===== 재료별 성능 =====
fprintf('%-14s %8s %8s %9s %10s %9s %8s %7s\n', ...
        'Material','k_x[N/m]','m_eff[g]','f1[Hz]','F@1mm[N]','sig[MPa]','SF','판정');
fprintf('%s\n', repmat('-',1,84));

results = struct();
for i = 1:size(mats,1)
    name = mats{i,1};  E = mats{i,2};  rho = mats{i,3};  Sy = mats{i,4};

    % 운동방향 강성
    k_x = 24*E*I/L^3;

    % 질량 (운동방향 모드 유효질량)
    m_ms    = rho * a_s^2 * b;            % moving stage (solid block)
    m_int_t = 2 * (rho * d_b * h_i * b);  % 중간 stage 2개 총질량
    m_bld_t = n_blade * (rho * L*t*b);    % blade 총질량
    m_eff   = m_ms + m_coil + m_pay ...
              + 0.25*m_int_t ...          % 중간stage: v/2 -> 유효 1/4
              + (1/3)*m_bld_t;            % blade: 근사 1/3

    % 공진
    f1 = (1/(2*pi)) * sqrt(k_x/m_eff);

    % 1mm 구동력 / 응력 / 안전율
    F_1mm  = k_x * x_req;
    sigma  = 3*E*t*(x_req/2)/L^2;         % per-stage 행정 x/2
    SF     = Sy / sigma;

    pass = (f1 > f_req) && (sigma <= 0.3*Sy);
    verdict = '미충족'; if pass, verdict = 'OK'; end

    fprintf('%-14s %8.0f %8.1f %9.1f %10.1f %9.1f %8.2f %7s\n', ...
            name, k_x, m_eff*1e3, f1, F_1mm, sigma*1e-6, SF, verdict);

    results(i).name=name; results(i).k_x=k_x; results(i).m_eff=m_eff;
    results(i).f1=f1; results(i).F_1mm=F_1mm; results(i).sigma=sigma; results(i).SF=SF;
end

fprintf('\n주: F@1mm = 1mm 변위에 필요한 VCM 추력. sig=blade 최대응력(0.3Sy 이하 권장). SF=Sy/sig.\n');
fprintf('   m_eff = stage + coil + (1/4)*중간stage + (1/3)*blade.\n');
