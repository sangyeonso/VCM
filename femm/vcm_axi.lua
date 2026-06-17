-- vcm_axi.lua  —  원통형 대칭 VCM 공극 자속밀도 B_g 검증 (FEMM 축대칭 자기정역학)
-- 실행: FEMM 실행 -> File > Open Lua Script (또는 콘솔에서) -> 이 파일 선택.
--       또는 명령창:  femm  ->  open("...\\vcm_axi.lua") 후 dofile.
-- 좌표: r=가로(축대칭 반경, x), z=세로(축, y). 축 r=0.
-- 자속경로: 자석(+z) -> 중앙요크 -> 공극(반경방향) -> 측면요크 -> 하단요크 -> 자석.
-- 코일은 무전류(자석만의 공극장 B_g 측정).
--
-- ===== CASE 선택: "paper"(논문 재현, ~108mT 기대) 또는 "opt"(우리 최적해) =====
CASE = "paper"     -- <-- 검증 먼저 "paper"로 ~108mT 확인 후, "opt"로 바꿔 재실행

if CASE == "paper" then            -- 2023 논문 최종설계(half 모델 치수, mm)
    r_m=3.0; h_m=8.0; h_y1=9.0; t_g=2.5; t_y3=1.0; h_y2=3.0
else                               -- 우리 최적해(Al, d_c=0.3) — 모델 B_g=450mT(가드)
    r_m=6.23; h_m=25.0; h_y1=5.1; t_g=2.77; t_y3=1.0; h_y2=3.0
end

Br = 1.17                          -- NdFeB N35 잔류자속밀도 [T]
mu_yoke = 1400                     -- 요크 비투자율(선형, 모델과 동일 가정)

-- ---- z 좌표 스택 ----
z0 = 0
z1 = h_y2                          -- 하단요크 위
z2 = z1 + h_m                      -- 자석 위(=중앙요크 아래)
z3 = z2 + h_y1                     -- 중앙요크 위(=측면요크 위)
r_side_i = r_m + t_g               -- 측면요크 내경
r_side_o = r_m + t_g + t_y3        -- 측면요크 외경

-- ===== 문서/물성 =====
newdocument(0)                     -- 0 = 자기(magnetics)
mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)

mi_getmaterial("Air")
-- NdFeB: Hc = Br/(mu0*mur), mur~1.05
Hc = Br/(4*math.pi*1e-7*1.05)
mi_addmaterial("NdFeB117", 1.05, 1.05, Hc, 0,0,0,0,0,0,0,0)
-- 요크: 선형 비투자율
mi_addmaterial("Yoke", mu_yoke, mu_yoke, 0, 0,0,0,0,0,0,0,0)

-- ===== 사각형 그리기 helper =====
function rect(r1,zz1,r2,zz2)
    mi_addnode(r1,zz1); mi_addnode(r2,zz1); mi_addnode(r2,zz2); mi_addnode(r1,zz2)
    mi_addsegment(r1,zz1,r2,zz1); mi_addsegment(r2,zz1,r2,zz2)
    mi_addsegment(r2,zz2,r1,zz2); mi_addsegment(r1,zz2,r1,zz1)
end

rect(0,       z0, r_side_o, z1)     -- 하단요크
rect(0,       z1, r_m,      z2)     -- 자석
rect(0,       z2, r_m,      z3)     -- 중앙요크
rect(r_side_i,z1, r_side_o, z3)     -- 측면요크
-- (공극 r_m..r_side_i 및 자석 옆 공간은 공기)

-- ===== 외부 공기 + 개방경계(ABC) =====
Rout = 4*r_side_o
mi_makeABC(7, Rout, 0, (z0+z3)/2, 0)   -- 점근경계(개방)

-- ===== 블록 라벨/물성 =====
function setblk(r,zz,mat,magdir,grp)
    mi_addblocklabel(r,zz); mi_selectlabel(r,zz)
    mi_setblockprop(mat, 1, 0, "<None>", magdir, grp, 0)
    mi_clearselected()
end
setblk(r_m*0.5,  (z0+z1)/2, "Yoke", 0, 1)              -- 하단요크
setblk(r_m*0.5,  (z1+z2)/2, "NdFeB117", 90, 2)         -- 자석(+z 축자화)
setblk(r_m*0.5,  (z2+z3)/2, "Yoke", 0, 1)              -- 중앙요크
setblk((r_side_i+r_side_o)/2, (z1+z3)/2, "Yoke", 0, 1) -- 측면요크
setblk((r_m+r_side_i)/2, (z2+z3)/2, "Air", 0, 0)       -- 공극(코일위치, 무전류)
setblk(Rout*0.8, (z0+z3)/2, "Air", 0, 0)               -- 외부 공기

-- ===== 해석 =====
mi_zoomnatural()
mi_saveas("vcm_axi_tmp.fem")
mi_analyze(1)
mi_loadsolution()

-- ===== 공극 자속밀도 B_g 측정 (중앙요크 높이, 중간반경에서 반경방향 Br 평균) =====
r_mid = r_m + t_g/2
N = 21
sumB = 0
for i=0,N-1 do
    zz = z2 + (z3-z2)*(i+0.5)/N
    Br_i, Bz_i = mo_getb(r_mid, zz)     -- 축대칭: (Bx,By)=(Br,Bz)
    sumB = sumB + math.abs(Br_i)
end
Bg = sumB/N

-- ===== 결과 출력 =====
msg = string.format(
 "CASE=%s\n r_m=%.2f h_m=%.2f h_y1=%.2f t_g=%.2f mm\n 공극 평균 |B_r| (B_g) = %.1f mT\n",
 CASE, r_m, h_m, h_y1, t_g, Bg*1000)
print(msg)
-- 파일로도 저장
f = openfile("vcm_axi_result.txt","a")
write(f, msg.."\n"); closefile(f)
-- showpointprops 로 화면 확인도 가능
-- mi_close()  -- 필요시 주석 해제
