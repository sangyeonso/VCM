# 보고서 그림 인덱스 (figure/)

생성 코드: `matlab/make_figures.m` (fig1~6), `matlab/fea_frame_modal.m` (fig7), 개념도는 직접 작성(fig0).

| 파일 | 내용 | 보고서 위치(제안) | 캡션(제안) |
|---|---|---|---|
| **fig_architecture.svg** | 전체 방법론 아키텍처 (모델→최적화→FEA→surrogate→최종) | 서론/방법론 개요 | 그림 0. 설계·해석 방법론 아키텍처: 해석모델 → SQP 최적화 → 구조·자기 FEA 검증 → FEA-정보 surrogate 재최적화 |
| **fig0_concept.svg** | Double Compound Parallelogram + VCM 평면 배치 개념도 | 1. 설계 개요 | 그림 1. VCM 구동 1자유도 대칭 Double Compound Parallelogram 유연기구 개념도 (X축 직동, moving-coil) |
| **fig7_fea_modeshape.png** | FEA 1차 모드형상 (133.8Hz, X-병진) | 4. FEA 검증 | 그림 2. 빔 FEM 1차 모드형상 — blade S자 변형으로 stage가 X 병진(가이드 운동) |
| **fig3_bg_parity.png** | $B_g$ surrogate vs FEMM 패리티 ($R^2$=0.99) | 4. FEA / 5. 최적화 | 그림 3. FEMM 48점으로 적합한 $B_g$ surrogate 정확도 |
| **fig4_bg_tg_lever.png** | $B_g$ vs 공극 $t_g$ (설계 레버) | 4·5 | 그림 4. 공극을 줄이면 $B_g$ 급증 — surrogate가 활용한 설계 레버 |
| **fig5_tradeoff.png** | 스트로크–공진 trade-off, 스펙영역 | 2. 타당성 / 5 | 그림 5. 변위·공진 상충과 권장설계 위치 |
| **fig2_dc_sensitivity.png** | $\gamma$ vs 코일선경 $d_c$ (임계 0.32mm) | 5. 결과 | 그림 6. 스펙 충족률–선경 관계: $d_c\lesssim0.32$mm에서 전 스펙 충족 |
| **fig1_material_gamma.png** | 재료 3종 $\gamma$ 비교 (Al 0.75) | 5. 재료비교 | 그림 7. 유연기구 재료별 스펙 충족률 (Al6061 최선) |
| **fig6_final_spec.png** | 최종설계($d_c$=0.3) 스펙 대비 달성 | 5. 최종설계 | 그림 8. 권장 최종설계의 요구 대비 달성 (스트로크·공진 충족, 전류 여유) |
| **fig8_femm_flux.png** | FEMM 자속밀도 분포 + 자속선 (권장 VCM) | 4. FEA / VCM | 그림 9. 권장 VCM의 FEMM 자속밀도(|B|)·자속선 — 자석→중앙요크→공극→측면요크 경로, 요크 비포화 |

## 메모
- 라벨은 영어, 캡션은 한글로(폰트 안정성). PNG는 200dpi.
- fig0은 SVG(벡터) — Word/LaTeX 직접 삽입 가능. PNG가 필요하면 말해주세요(렌더 가능).
- (선택) FEMM 자속분포 그림이 필요하면 pyFEMM `mo_savebitmap`로 추가 생성 가능.
