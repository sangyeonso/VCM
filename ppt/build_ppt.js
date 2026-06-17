// build_ppt.js  — 초정밀설계 프로젝트 발표 (BMW 코퍼레이트 톤)
const pptxgen = require("pptxgenjs");
const P = new pptxgen();
P.layout = "LAYOUT_16x9";              // 10 x 5.625 in
P.author = "소상연";
P.title  = "VCM 구동 DCP 유연기구 스테이지";

const FIG = "../figure/";
// ---- BMW 코퍼레이트 토큰 ----
const C = {
  canvas:"FFFFFF", soft:"F7F7F7", card:"FAFAFA", dark:"1A2129", darkEl:"262E38",
  primary:"1C69D4", primaryA:"0653B6",
  ink:"262626", body:"3C3C3C", bodyStrong:"1A1A1A", muted:"6B6B6B", mutedSoft:"9A9A9A",
  onDark:"FFFFFF", onDarkSoft:"BBBBBB", hair:"E6E6E6", hairStrong:"CCCCCC",
  mBlueL:"0066B1", mBlueD:"1C69D4", mRed:"E22718", success:"22C55E",
};
const DISP="Inter", LIGHT="Inter Light";   // BMW Type Next 대체

function bg(s,dark){ s.background={color: dark?C.dark:C.canvas}; }
// UPPERCASE 트래킹 라벨 (BMW 시그니처)
function label(s,t,x,y,dark){ s.addText(t,{x,y,w:6,h:0.3,margin:0,fontFace:DISP,bold:true,
  fontSize:11,charSpacing:3,color: dark?C.onDarkSoft:C.primary, align:"left"}); }
function footer(s,n,dark){
  s.addText("초정밀시스템설계 프로젝트 · 소상연",{x:0.6,y:5.28,w:6,h:0.25,margin:0,
    fontFace:LIGHT,fontSize:9,color: dark?C.onDarkSoft:C.mutedSoft,align:"left"});
  s.addText(String(n).padStart(2,"0"),{x:9.0,y:5.28,w:0.4,h:0.25,margin:0,
    fontFace:DISP,bold:true,fontSize:9,color: dark?C.onDarkSoft:C.muted,align:"right"});
}
const img=(s,f,x,y,w,h)=> s.addImage({path:FIG+f,x,y,w,h,sizing:{type:"contain",w,h}});

/* ============ S1 — Title (DARK) ============ */
let s=P.addSlide(); bg(s,true);
label(s,"PRECISION MECHATRONICS  ·  TERM PROJECT",0.7,0.7,true);
s.addText([
 {text:"VCM 구동 1자유도",options:{breakLine:true}},
 {text:"Double Compound Parallelogram",options:{breakLine:true,color:C.primary}},
 {text:"유연기구 스테이지",options:{}},
],{x:0.65,y:1.5,w:6.2,h:2.2,margin:0,fontFace:DISP,bold:true,fontSize:34,color:C.onDark,lineSpacing:42,align:"left"});
s.addText("설계 · 최적화 · FEA 검증",{x:0.7,y:3.7,w:6,h:0.4,margin:0,fontFace:LIGHT,fontSize:18,color:C.onDarkSoft});
s.addText("소상연  ·  202624279  ·  아주대학교 D.N.A.플러스융합학과",
  {x:0.7,y:4.95,w:7,h:0.3,margin:0,fontFace:LIGHT,fontSize:12,color:C.onDarkSoft});
// 블루 스탯 패널
s.addShape(P.shapes.RECTANGLE,{x:7.35,y:1.5,w:2.0,h:2.6,fill:{color:C.primary}});
const stat=(t,l,yy)=>{ s.addText(t,{x:7.5,y:yy,w:1.75,h:0.5,margin:0,fontFace:DISP,bold:true,fontSize:26,color:C.onDark});
  s.addText(l,{x:7.5,y:yy+0.5,w:1.75,h:0.3,margin:0,fontFace:LIGHT,fontSize:11,color:"DCE8FB"}); };
stat("105 Hz","1차 공진",1.7); stat("1.05 mm","최대변위 @2A",2.55); stat("Al6061","최적 재료",3.4);

/* ============ S2 — 배경·요구사양 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"OVERVIEW",0.7,0.55,false);
s.addText("초정밀 위치결정 스테이지",{x:0.65,y:0.9,w:8.5,h:0.7,margin:0,fontFace:DISP,bold:true,fontSize:32,color:C.ink});
s.addText([
 {text:"원자현미경·광학 스캐너 등은 나노미터급 분해능의 위치결정을 요구한다.",options:{breakLine:true}},
 {text:"유연기구(flexure)는 마모·백래시가 없고 단일체 가공으로 조립오차가 없으며,",options:{breakLine:true}},
 {text:"변형이 탄성영역에서 반복·예측 가능해 초정밀 응용에 적합하다.",options:{}},
],{x:0.7,y:1.85,w:4.8,h:1.6,margin:0,fontFace:LIGHT,fontSize:15,color:C.body,lineSpacing:24});
s.addText("목표 — VCM으로 구동되는 1자유도 유연기구 스테이지를 요구사양 충족하도록 설계·최적화",
  {x:0.7,y:3.5,w:4.8,h:0.8,margin:0,fontFace:DISP,bold:true,fontSize:14,color:C.bodyStrong,lineSpacing:20});
// 요구사양 표
s.addText("REQUIREMENTS",{x:5.9,y:1.0,w:3,h:0.3,margin:0,fontFace:DISP,bold:true,fontSize:11,charSpacing:2,color:C.primary});
const reqRows=[["크기 (구동기 포함)","≤ 200×200×20 mm"],["최대 변위","> 1 mm"],
  ["1차 공진주파수","> 100 Hz"],["최대 전류","= 2 A"],["코일 최소 직경","> 0.5 mm"]];
s.addTable(reqRows.map(r=>[
  {text:r[0],options:{fontFace:LIGHT,fontSize:13,color:C.body,align:"left",valign:"middle"}},
  {text:r[1],options:{fontFace:DISP,bold:true,fontSize:13,color:C.ink,align:"right",valign:"middle"}}]),
  {x:5.9,y:1.45,w:3.5,colW:[2.1,1.4],rowH:0.5,border:{type:"solid",pt:1,color:C.hair},fill:{color:C.canvas}});
footer(s,2,false);

/* ============ S3 — 방법론 아키텍처 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"METHODOLOGY",0.7,0.5,false);
s.addText("설계 · 해석 방법론 아키텍처",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:28,color:C.ink});
img(s,"fig_architecture_final.png",0.5,1.65,9.0,3.4);
s.addText("해석모델  →  SQP 최적화  →  구조·자기 FEA 검증  →  FEA-정보 surrogate 재최적화",
  {x:0.7,y:5.0,w:9,h:0.3,margin:0,fontFace:LIGHT,fontSize:12,color:C.muted,align:"center"});
footer(s,3,false);

/* ============ S4 — 개념설계 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"CONCEPT DESIGN",0.7,0.5,false);
s.addText("대칭 Double Compound Parallelogram + VCM",{x:0.65,y:0.85,w:8.8,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:26,color:C.ink});
s.addText([
 {text:"유연기구",options:{bold:true,color:C.primary,breakLine:true,fontFace:DISP}},
 {text:"대칭형 직선가이드 DCP — 운동방향 강성 kx = 24EI/L³,",options:{breakLine:true}},
 {text:"자유도 1로 구속·기생운동/열오차 보상.",options:{breakLine:true}},
 {text:" ",options:{breakLine:true,fontSize:8}},
 {text:"구동기",options:{bold:true,color:C.primary,breakLine:true,fontFace:DISP}},
 {text:"원통형 대칭 moving-coil VCM — 이동질량 최소화,",options:{breakLine:true}},
 {text:"구동축을 운동축 X에 정확히 정렬.",options:{breakLine:true}},
 {text:" ",options:{breakLine:true,fontSize:8}},
 {text:"설계변수: 유연기구 {t, L, b}, VCM {r_m, h_m, h_y1, t_g, d_c}",options:{color:C.muted}},
],{x:0.7,y:1.7,w:4.4,h:3.1,margin:0,fontFace:LIGHT,fontSize:14,color:C.body,lineSpacing:21});
s.addShape(P.shapes.RECTANGLE,{x:5.45,y:1.6,w:4.0,h:3.3,fill:{color:C.card}});
img(s,"fig0_concept.png",5.55,1.7,3.8,3.1);
footer(s,4,false);

/* ============ S5 — 해석모델 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"ANALYTICAL MODELS",0.7,0.5,false);
s.addText("해석 모델 — 결합 시스템",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:28,color:C.ink});
// 카드 1: 유연기구
s.addShape(P.shapes.RECTANGLE,{x:0.65,y:1.7,w:4.3,h:2.5,fill:{color:C.card},line:{color:C.hair,width:1}});
s.addText("유연기구 모델",{x:0.9,y:1.95,w:3.8,h:0.4,margin:0,fontFace:DISP,bold:true,fontSize:18,color:C.ink});
s.addText([
 {text:"k_x = 24 E I / L³,   x = F / k_x",options:{breakLine:true,fontFace:DISP,bold:true,color:C.primary}},
 {text:"f₁ = (1/2π)·√(k_x / m_eff)",options:{breakLine:true,fontFace:DISP,bold:true,color:C.primary}},
 {text:"검증: 행렬법 = 12EI/L³, 직병렬 = 24EI/L³ (오차 ~10⁻²³)",options:{color:C.body}},
],{x:0.9,y:2.5,w:3.85,h:1.5,margin:0,fontFace:LIGHT,fontSize:13.5,color:C.body,lineSpacing:24});
// 카드 2: VCM
s.addShape(P.shapes.RECTANGLE,{x:5.15,y:1.7,w:4.3,h:2.5,fill:{color:C.card},line:{color:C.hair,width:1}});
s.addText("VCM Permeance 모델",{x:5.4,y:1.95,w:3.8,h:0.4,margin:0,fontFace:DISP,bold:true,fontSize:18,color:C.ink});
s.addText([
 {text:"자기등가회로 → 공극자속 B_g = Φ_g / A_g",options:{breakLine:true,fontFace:DISP,bold:true,color:C.primary}},
 {text:"추력 F = n_eff · B_g · L_m · I",options:{breakLine:true,fontFace:DISP,bold:true,color:C.primary}},
 {text:"검증: 2023 논문 재현 → 힘상수 Kf 오차 −6%",options:{color:C.body}},
],{x:5.4,y:2.5,w:3.85,h:1.5,margin:0,fontFace:LIGHT,fontSize:13.5,color:C.body,lineSpacing:24});
s.addText("시스템 결합:  코일질량 → m_eff (공진)   ·   추력 Kf·I → 변위",
  {x:0.65,y:4.5,w:9,h:0.4,margin:0,fontFace:DISP,bold:true,fontSize:14,color:C.bodyStrong,align:"center"});
footer(s,5,false);

/* ============ S6 — 근본 물리 한계 (DARK, key insight) ============ */
s=P.addSlide(); bg(s,true); label(s,"KEY INSIGHT  ·  FUNDAMENTAL LIMIT",0.7,0.6,true);
s.addText("스트로크와 공진을 결합하면 힘상수가 소거된다",{x:0.65,y:1.1,w:8.8,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:26,color:C.onDark});
s.addShape(P.shapes.RECTANGLE,{x:1.3,y:2.1,w:7.4,h:1.2,fill:{color:C.darkEl}});
s.addText("ρ_cu · (π/4) · d_c²  /  B_g    ≤    I_max  /  [ x_req · (2π f_req)² ]",
  {x:1.3,y:2.1,w:7.4,h:1.2,margin:0,fontFace:DISP,bold:true,fontSize:18,color:C.primary,align:"center",valign:"middle"});
s.addText([
 {text:"VCM을 아무리 강화해도 충족 여부는 ",options:{}},
 {text:"코일 선경 d_c 와 B_g 만으로 결정",options:{bold:true,color:C.onDark,fontFace:DISP}},
 {text:"된다.",options:{}},
],{x:1.3,y:3.6,w:7.4,h:0.4,margin:0,fontFace:LIGHT,fontSize:16,color:C.onDarkSoft,align:"center"});
s.addText("→  d_c > 0.5 mm 강제 시 변위·공진 동시 충족은 본질적으로 불가  (임계 d_c ≈ 0.32 mm)",
  {x:1.3,y:4.25,w:7.4,h:0.4,margin:0,fontFace:LIGHT,fontSize:14,color:C.onDarkSoft,align:"center"});
footer(s,6,true);

/* ============ S7 — 최적화 결과 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"OPTIMIZATION  ·  fmincon SQP",0.7,0.5,false);
s.addText("재료 비교 & 코일 선경 민감도",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:28,color:C.ink});
img(s,"fig1_material_gamma.png",0.6,1.7,4.3,2.9);
img(s,"fig2_dc_sensitivity.png",5.1,1.7,4.4,2.9);
s.addText("Al6061 최선 (γ=0.75)   ·   d_c ≲ 0.32 mm 에서 전 스펙 충족",
  {x:0.65,y:4.7,w:9,h:0.3,margin:0,fontFace:LIGHT,fontSize:12,color:C.muted,align:"center"});
footer(s,7,false);

/* ============ S8 — 구조 FEA (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"FEA VERIFICATION  ·  STRUCTURAL",0.7,0.5,false);
s.addText("구조 FEA — 빔 프레임 모달해석",{x:0.65,y:0.85,w:8.8,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:28,color:C.ink});
img(s,"fig7_fea_modeshape.png",5.7,1.5,3.6,3.5);
s.addText([
 {text:"2D 프레임 빔 유한요소(consistent mass)로",options:{breakLine:true}},
 {text:"Double Compound Parallelogram을 모달해석.",options:{breakLine:true}},
 {text:" ",options:{breakLine:true,fontSize:8}},
 {text:"1차 모드는 blade의 S자 변형을 통한",options:{breakLine:true}},
 {text:"stage의 X 병진(가이드 운동).",options:{}},
],{x:0.7,y:1.9,w:4.6,h:2.0,margin:0,fontFace:LIGHT,fontSize:15,color:C.body,lineSpacing:23});
s.addShape(P.shapes.RECTANGLE,{x:0.7,y:3.95,w:4.6,h:0.9,fill:{color:C.dark}});
s.addText([
 {text:"FEA 133.8 Hz  vs  해석 137.7 Hz  ",options:{color:C.onDark,fontFace:LIGHT}},
 {text:"→  −2.8%",options:{bold:true,color:C.primary,fontFace:DISP}},
],{x:0.7,y:3.95,w:4.6,h:0.9,margin:0,fontSize:16,align:"center",valign:"middle"});
footer(s,8,false);

/* ============ S9 — 자기 FEA & surrogate (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"FEA VERIFICATION  ·  MAGNETIC & SURROGATE",0.7,0.5,false);
s.addText("자기 FEA (FEMM) 와 FEA-정보 surrogate",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:25,color:C.ink});
img(s,"fig8_femm_flux.png",0.7,1.6,2.6,3.4);
img(s,"fig3_bg_parity.png",3.7,1.7,3.0,3.0);
s.addText([
 {text:"FEMM 축대칭 해석으로 해석모델이 B_g를",options:{breakLine:true}},
 {text:"과대평가함을 발견 ",options:{}},
 {text:"(450 → 291 mT).",options:{bold:true,color:C.ink,fontFace:DISP,breakLine:true}},
 {text:" ",options:{breakLine:true,fontSize:8}},
 {text:"pyFEMM 48점 격자로 B_g(r_m,h_y1,t_g)",options:{breakLine:true}},
 {text:"surrogate 적합 ",options:{}},
 {text:"(R² = 0.99)",options:{bold:true,color:C.primary,fontFace:DISP,breakLine:true}},
 {text:"→ 최적화에 대체, FEA-정확 재최적화.",options:{}},
],{x:6.95,y:1.8,w:2.6,h:3.0,margin:0,fontFace:LIGHT,fontSize:13,color:C.body,lineSpacing:20});
footer(s,9,false);

/* ============ S10 — 동특성·제어 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"DYNAMICS & CONTROL  ·  SIMULINK",0.7,0.5,false);
s.addText("폐루프 위치제어 (PID + 2A 전류포화)",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:27,color:C.ink});
img(s,"fig11_step.png",0.6,1.7,4.3,2.9);
img(s,"fig12_scan.png",5.1,1.7,4.4,2.9);
s.addText("1 mm 스텝: 정착 ~64 ms (전류한계 근처)   ·   10 Hz 스캐닝 추종 RMS 82 μm",
  {x:0.65,y:4.7,w:9,h:0.3,margin:0,fontFace:LIGHT,fontSize:12,color:C.muted,align:"center"});
footer(s,10,false);

/* ============ S11 — Simscape 물리모델 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"PHYSICAL MODEL  ·  SIMSCAPE & ANIMATION",0.7,0.5,false);
s.addText("Simscape 물리네트워크 & 모션 시각화",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:26,color:C.ink});
s.addShape(P.shapes.RECTANGLE,{x:5.5,y:1.55,w:4.0,h:3.4,fill:{color:C.card}});
img(s,"fig13_frame.png",5.6,1.6,3.8,3.3);
s.addText([
 {text:"전기(R–L–제어전압원) – 기계(질량–스프링–댐퍼)를",options:{breakLine:true}},
 {text:"전류·운동센서로 전자기 결합한 물리네트워크.",options:{breakLine:true}},
 {text:" ",options:{breakLine:true,fontSize:8}},
 {text:"검증: 1V 스텝 정상상태 변위 ",options:{}},
 {text:"0.204 mm",options:{bold:true,color:C.primary,fontFace:DISP}},
 {text:" — 이론·블록선도 모델과 일치, 속도 기반",options:{breakLine:true}},
 {text:"역기전력으로 안정 감쇠.",options:{breakLine:true}},
 {text:" ",options:{breakLine:true,fontSize:8}},
 {text:"애니메이션: DCP가 X로 스캐닝하며 blade가 S자로 휨.",options:{color:C.muted}},
],{x:0.7,y:1.75,w:4.5,h:3.1,margin:0,fontFace:LIGHT,fontSize:14,color:C.body,lineSpacing:21});
footer(s,11,false);

/* ============ S12 — 최종 설계안 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"FINAL DESIGN",0.7,0.5,false);
s.addText("최종 설계안 (Al6061, FEA-정확)",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:28,color:C.ink});
const hd={fontFace:DISP,bold:true,fontSize:12,color:C.onDark,fill:{color:C.dark},align:"center",valign:"middle"};
const cel=(t,b)=>({text:t,options:{fontFace:b?DISP:LIGHT,bold:!!b,fontSize:13,color:b?C.primary:C.body,align:"center",valign:"middle"}});
s.addTable([
 [{text:"d_c [mm]",options:hd},{text:"f₁ [Hz]",options:hd},{text:"변위@2A [mm]",options:hd},{text:"1mm 전류 [A]",options:hd},{text:"판정",options:hd}],
 [cel("0.50"),cel("75.3"),cel("0.75"),cel("2.65"),cel("미달")],
 [cel("0.40"),cel("87.4"),cel("0.87"),cel("2.29"),cel("~87%")],
 [cel("0.30",1),cel("105.3",1),cel("1.05",1),cel("1.90",1),cel("충족",1)],
],{x:0.7,y:1.7,colW:[0.9,0.85,1.35,1.3,0.8],rowH:0.5,border:{type:"solid",pt:1,color:C.hair}});
s.addText([
 {text:"권장 설계  d_c = 0.30 mm",options:{breakLine:true,bold:true,fontFace:DISP,fontSize:16,color:C.ink}},
 {text:"유연기구 t=0.43, L=24.2, b=8 mm",options:{breakLine:true,color:C.body}},
 {text:"VCM r_m=6.5, h_m=25, h_y1=5.8, t_g=2.5 mm (Ø20)",options:{breakLine:true,color:C.body}},
 {text:"응력 안전율 3.33 · 전 요구사양 충족",options:{color:C.muted}},
],{x:0.7,y:3.6,w:5.4,h:1.2,margin:0,fontFace:LIGHT,fontSize:13,lineSpacing:20});
img(s,"fig6_final_spec.png",6.3,1.9,3.2,2.6);
footer(s,12,false);

/* ============ S13 — 대안·한계 (LIGHT) ============ */
s=P.addSlide(); bg(s,false); label(s,"ALTERNATIVES & LIMITS",0.7,0.5,false);
s.addText("미충족 시 대안 & 한계",{x:0.65,y:0.85,w:9,h:0.6,margin:0,fontFace:DISP,bold:true,fontSize:28,color:C.ink});
const alt=(k,t,yy)=>{ s.addShape(P.shapes.RECTANGLE,{x:0.7,y:yy,w:0.45,h:0.45,fill:{color:C.primary}});
  s.addText(k,{x:0.7,y:yy,w:0.45,h:0.45,margin:0,fontFace:DISP,bold:true,fontSize:16,color:C.onDark,align:"center",valign:"middle"});
  s.addText(t,{x:1.3,y:yy,w:8.0,h:0.45,margin:0,fontFace:LIGHT,fontSize:14,color:C.body,valign:"middle"}); };
alt("A","코일 선경을 0.3–0.4 mm로 완화  (권장 — 전 스펙 충족)",1.75);
alt("B","한쪽 사양 절충 — 공진 또는 변위 양보",2.4);
alt("C","경량 유연기구(Mg·중공) + 강한 자기회로(N52·FeCo) 조합 → γ ≈ 0.99",3.05);
alt("D","전류 상향 (전원·발열 부담 동반)",3.7);
s.addText("한계 — 권장 설계는 1 mm 유지 시 ~9.3 W 발열로 방열 설계 필요. 레버·브리지 증폭은 lost motion·공진 저하로 비권장.",
  {x:0.7,y:4.5,w:8.7,h:0.6,margin:0,fontFace:LIGHT,fontSize:12.5,color:C.muted,lineSpacing:18});
footer(s,13,false);

/* ============ S14 — 결론 (DARK) ============ */
s=P.addSlide(); bg(s,true); label(s,"CONCLUSION",0.7,0.6,true);
s.addText("결론",{x:0.65,y:1.05,w:8,h:0.7,margin:0,fontFace:DISP,bold:true,fontSize:34,color:C.onDark});
const con=(n,t,yy)=>{ s.addText(n,{x:0.7,y:yy,w:0.6,h:0.5,margin:0,fontFace:DISP,bold:true,fontSize:22,color:C.primary});
  s.addText(t,{x:1.4,y:yy,w:8.0,h:0.6,margin:0,fontFace:LIGHT,fontSize:15,color:C.onDark,valign:"top",lineSpacing:20}); };
con("01","Al6061이 최적이며 권장 설계(d_c=0.3)는 f₁=105 Hz, 1.05 mm@2A로 전 스펙 충족",2.1);
con("02","해석모델 → FEMM 검증 → FEA-정보 surrogate 3단 정밀화로 추력 추정 정확화 (B_g 오차 160→12 mT)",2.85);
con("03","d_c>0.5 mm 강제 시 동시 충족 불가를 물리적 한계식으로 규명하고 대안 제시",3.75);
s.addText("감사합니다.",{x:0.7,y:4.6,w:6,h:0.5,margin:0,fontFace:DISP,bold:true,fontSize:20,color:C.onDarkSoft});
footer(s,14,true);

P.writeFile({fileName:"소상연_초정밀시스템설계_발표.pptx"}).then(f=>console.log("saved",f));
