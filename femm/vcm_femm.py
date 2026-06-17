"""
vcm_femm.py  -  pyFEMM으로 원통형 VCM 공극 자속밀도 B_g 자동 해석 (VS Code 실행용)
실행:  (GPSM 환경)  python vcm_femm.py
  - FEMM을 백그라운드로 띄워 여러 설계를 루프로 해석, 결과를 콘솔 + CSV로 출력.
  - GUI 조작/CASE 수동변경 불필요. DESIGNS / SWEEP만 수정해 '뺑뺑이' 돌리면 됨.
필요: FEMM 4.2 설치(C:\\femm42), pyfemm, pywin32  (GPSM venv에 설치됨)
"""
import os, csv, math, sys
try:
    sys.stdout.reconfigure(encoding="utf-8")   # 콘솔 한글 깨짐 방지
except Exception:
    pass
import femm

TMP_FEM = r"C:\femm_vcm\vcm_tmp.fem"          # ASCII 임시경로(FEMM은 한글경로 싫어함)
OUT_CSV = os.path.join(os.path.dirname(__file__), "vcm_femm_results.csv")
os.makedirs(r"C:\femm_vcm", exist_ok=True)

MU0 = 4*math.pi*1e-7
BR  = 1.17          # NdFeB N35 [T]
MU_YOKE = 1400      # 요크 선형 비투자율(해석모델과 동일)

def compute_Bg(r_m, h_m, h_y1, t_g, t_y3=1.0, h_y2=3.0):
    """축대칭 자기정역학으로 공극 평균 |Br| [mT] 반환. 치수 단위 mm."""
    z0, z1 = 0.0, h_y2
    z2 = z1 + h_m
    z3 = z2 + h_y1
    r_si = r_m + t_g            # 측면요크 내경
    r_so = r_m + t_g + t_y3     # 측면요크 외경

    femm.newdocument(0)         # 0 = 자기
    femm.mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)
    femm.mi_getmaterial("Air")
    Hc = BR/(MU0*1.05)
    femm.mi_addmaterial("NdFeB117", 1.05, 1.05, Hc, 0,0,0,0,0,0,0,0)
    femm.mi_addmaterial("Yoke", MU_YOKE, MU_YOKE, 0, 0,0,0,0,0,0,0,0)

    def rect(r1, zz1, r2, zz2):
        femm.mi_addnode(r1,zz1); femm.mi_addnode(r2,zz1)
        femm.mi_addnode(r2,zz2); femm.mi_addnode(r1,zz2)
        femm.mi_addsegment(r1,zz1,r2,zz1); femm.mi_addsegment(r2,zz1,r2,zz2)
        femm.mi_addsegment(r2,zz2,r1,zz2); femm.mi_addsegment(r1,zz2,r1,zz1)

    rect(0,    z0, r_so, z1)    # 하단요크
    rect(0,    z1, r_m,  z2)    # 자석
    rect(0,    z2, r_m,  z3)    # 중앙요크
    rect(r_si, z1, r_so, z3)    # 측면요크

    Rout = 4*r_so
    femm.mi_makeABC(7, Rout, 0, (z0+z3)/2, 0)

    def blk(r, zz, mat, magdir, grp):
        femm.mi_addblocklabel(r, zz); femm.mi_selectlabel(r, zz)
        femm.mi_setblockprop(mat, 1, 0, "<None>", magdir, grp, 0)
        femm.mi_clearselected()

    blk(r_m*0.5,        (z0+z1)/2, "Yoke",     0,  1)
    blk(r_m*0.5,        (z1+z2)/2, "NdFeB117", 90, 2)   # +z 축자화
    blk(r_m*0.5,        (z2+z3)/2, "Yoke",     0,  1)
    blk((r_si+r_so)/2,  (z1+z3)/2, "Yoke",     0,  1)
    # 공극은 바깥 공기와 이어진 한 영역 -> Air 라벨 1개로 충분(이중라벨 방지)
    blk(Rout*0.8,       (z0+z3)/2, "Air",      0,  0)

    femm.mi_zoomnatural()
    femm.mi_saveas(TMP_FEM)
    femm.mi_analyze(1)
    femm.mi_loadsolution()

    # 공극 평균 |Br| (중앙요크 높이, 중간반경)
    r_mid = r_m + t_g/2
    N, s = 21, 0.0
    for i in range(N):
        zz = z2 + (z3-z2)*(i+0.5)/N
        vals = femm.mo_getpointvalues(r_mid, zz)   # (A, B1=Br, B2=Bz, ...)
        s += abs(vals[1])
    return s/N*1000.0  # mT


# ===== 해석할 설계들 (여기만 수정해서 뺑뺑이) =====
DESIGNS = {
    # name        : (r_m,  h_m, h_y1, t_g)
    "paper(논문)"  : (3.00,  8.0,  9.0, 2.50),   # 기대 ~108mT(논문 FEA 105)
    "opt(5단계)"   : (6.23, 25.0,  5.1, 2.77),   # 모델 450mT -> FEMM?
    "final_dc0.3" : (6.04, 25.0,  7.5, 2.96),   # 보정 최종 권장 VCM
}

# ===== (선택) 파라미터 스윕: 공극 t_g 변화에 따른 B_g =====
SWEEP_ENABLE = True
SWEEP_TG = [1.5, 2.0, 2.5, 3.0]   # mm, opt 형상에서 t_g만 변경

def main():
    femm.openfemm(1)        # 1 = 창 숨김(백그라운드)
    rows = []
    print(f"{'설계':<14}{'r_m':>6}{'h_m':>6}{'h_y1':>6}{'t_g':>6}{'B_g[mT]':>10}")
    print("-"*52)
    for name, (r_m, h_m, h_y1, t_g) in DESIGNS.items():
        Bg = compute_Bg(r_m, h_m, h_y1, t_g)
        print(f"{name:<14}{r_m:>6.2f}{h_m:>6.1f}{h_y1:>6.1f}{t_g:>6.2f}{Bg:>10.1f}")
        rows.append(["design", name, r_m, h_m, h_y1, t_g, round(Bg,2)])

    if SWEEP_ENABLE:
        print("\n[스윕] opt 형상에서 t_g 변화:")
        r_m, h_m, h_y1 = 6.23, 25.0, 5.1
        for tg in SWEEP_TG:
            Bg = compute_Bg(r_m, h_m, h_y1, tg)
            print(f"  t_g={tg:.2f} mm -> B_g={Bg:.1f} mT")
            rows.append(["sweep_tg", f"tg{tg}", r_m, h_m, h_y1, tg, round(Bg,2)])

    femm.closefemm()
    with open(OUT_CSV, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["kind","name","r_m","h_m","h_y1","t_g","Bg_mT"])
        w.writerows(rows)
    print(f"\n결과 저장: {OUT_CSV}")

if __name__ == "__main__":
    main()
