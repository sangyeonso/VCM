# VCM-구동 DCP 유연 스테이지 · 1-DOF 정밀 직선 스테이지

> **초정밀시스템설계 Term Project** · 아주대학교 D.N.A.플러스융합학과
> 소상연 (202624279)

보이스코일 모터(VCM)로 직접 구동되는 1자유도 **이중 복합 평행판(Double Compound Parallelogram, DCP)** 유연 스테이지의 설계·해석·검증 기록입니다. 개념설계 → 해석모델 → 물리한계 분석 → 최적화 → FEA 검증 → 동특성/제어의 흐름으로 진행했으며, 상용 FEA 없이 **MATLAB 보 FEM(모달) + 무료 FEMM(자기장) + 다항 대리모델**의 3중 검증으로 해석값을 보정했습니다.

🔗 **웹 워크북:** https://sangyeonso.github.io/VCM/ *(GitHub Pages 활성화 후)*

---

## 핵심 결과 (KPI)

| 항목 | 목표 | 달성 | 판정 |
|---|---|---|---|
| 1차 공진 $f_1$ | > 100 Hz | **105.3 Hz** | ✅ |
| 행정 (stroke) @ 2 A | > 1 mm | **1.05 mm** | ✅ |
| 행정 1 mm 도달 전류 | ≤ 2 A | **1.9 A** | ✅ |
| 외형 | ≤ 200×200×20 mm | 충족 | ✅ |
| 코일 와이어경 $d_c$ | > 0.5 mm | 0.3 mm | ⚠️ 완화* |

\* 힘상수가 소거되는 물리 한계식상 $d_c>0.5$ mm 를 글자 그대로 지키면 행정·공진 동시충족 해가 **존재하지 않음**을 증명하고, 물리적 불가피성에 근거해 완화. (`notes/[5단계]` 참조)

**최종 설계:** 재료 Al6061 · blade $t=0.43$ / $L=24.2$ / $b=8$ mm · VCM $r_m=6.5$ / $h_m=25$ / $h_{y1}=5.8$ / $t_g=2.5$ mm (Ø20)

---

## 저장소 구조

```
VCM/
├── matlab/      해석·최적화 MATLAB 스크립트
│   ├── flexure_model.m      DCP 강성·공진·응력 모델
│   ├── vcm_model.m          VCM 퍼미언스 자기회로 모델
│   ├── opt_design.m         γ-최대화 최적화 (fmincon SQP, multi-start)
│   ├── final_design.m       최종 설계점 산출
│   ├── fea_frame_modal.m    보 요소 FEM 모달 해석
│   ├── material_study.m     Al6061 / 스프링강 / BeCu 비교
│   ├── make_figures.m       결과 그림 생성
│   └── make_concept_fig.m   개념도 생성
├── femm/        자기장 FEM (FEMM 4.2 + pyFEMM)
│   ├── vcm_axi.lua          축대칭 VCM 자기 해석 스크립트
│   ├── vcm_femm.py          pyFEMM 자동 스윕 (48점 그리드)
│   ├── vcm_fluxplot.py      자속 분포 플롯
│   └── bg_surrogate.py      B_g 2차 다항 대리모델 (R²=0.99)
├── simulink/    동특성·제어 모델
│   ├── build_vcm_plant.m    Model 1: 개루프 플랜트
│   ├── build_vcm_control.m  Model 2: PID 폐루프
│   ├── build_vcm_simscape.m Model 3: Simscape 물리 네트워크
│   ├── animate_flexure.m    X-스캐닝 애니메이션 (GIF)
│   └── *.slx                생성된 Simulink 모델
├── figure/      결과 그림 (png / svg / gif)
├── report/      main_cvpr.tex  (CVPR 2단 스타일, 한글 kotex/XeLaTeX)
├── ppt/         build_ppt.js   (pptxgenjs, BMW 코퍼레이트 톤)
├── notes/       설계 단계별 노트 ([1·3·4·5단계], 보고서 골격)
└── docs/        GitHub Pages 웹 워크북 (index.html)
```

---

## 재현 방법

**MATLAB** (R2022b+, Optimization·Global Optimization·Control System·Simulink·Simscape Toolbox):

```matlab
cd matlab
flexure_model            % 유연기구 모델 검증
vcm_model                % VCM 자기회로 모델
material_study           % 재료 비교
opt_design               % γ-최대화 최적화
final_design             % 최종 설계점
make_figures             % 결과 그림 일괄 생성

cd ../simulink
build_vcm_plant          % 개루프
build_vcm_control        % PID 폐루프
build_vcm_simscape       % Simscape 물리 네트워크
```

**FEMM** (4.2): `femm/vcm_axi.lua` 를 FEMM에서 실행하거나, pyFEMM 환경에서 `python femm/vcm_femm.py` 로 자동 스윕 → `bg_surrogate.py` 로 대리모델 적합.

**보고서:** `report/main_cvpr.tex` 를 XeLaTeX(+kotex)로 컴파일. `figure/` 의 그림을 참조합니다.

---

## 사용 도구

MATLAB R2022b · Simulink / Simscape · FEMM 4.2 + pyFEMM · XeLaTeX(kotex) · pptxgenjs

> 강의 자료(슬라이드 PDF)·참고 논문 원본은 저작권상 본 저장소에 포함하지 않았습니다.
