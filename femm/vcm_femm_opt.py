"""
vcm_femm_opt.py - N52 VCM의 FEMM 기반 형상 최적화 + 해석모델 검증
  - 검증된 셋업(축대칭, 요크 mu=1400)에서 자석만 N35->N52(Br=1.43)로 변경
  - (1) 해석모델 '최적'(B_g=2T 주장) 형상을 FEMM으로 검증 -> 진짜 B_g
  - (2) 형상 DOE(r_m, h_y1) 스윕 -> FEMM 기준 max B_g, max K_f/m_coil
  실행: GPSM venv python vcm_femm_opt.py
"""
import os, math, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import femm

TMP_FEM = r"C:\femm_vcm\vcm_opt.fem"
os.makedirs(r"C:\femm_vcm", exist_ok=True)
MU0 = 4*math.pi*1e-7
BR  = 1.43          # N52 [T]
MU_YOKE = 1400      # 검증된 값 유지(보수적). FeCo면 약간 더 나옴
DC = 0.5e-3; RP=1.15; CC=0.3e-3; TB=0.3e-3; NSYM=2; I=2.0

def compute_Bg(r_m, h_m, h_y1, t_g, t_y3=0.8, h_y2=4.0):
    z0, z1 = 0.0, h_y2; z2 = z1+h_m; z3 = z2+h_y1
    r_si = r_m+t_g; r_so = r_m+t_g+t_y3
    femm.newdocument(0)
    femm.mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)
    femm.mi_getmaterial("Air")
    Hc = BR/(MU0*1.05)
    femm.mi_addmaterial("NdFeB", 1.05, 1.05, Hc, 0,0,0,0,0,0,0,0)
    femm.mi_addmaterial("Yoke", MU_YOKE, MU_YOKE, 0, 0,0,0,0,0,0,0,0)
    def rect(r1,zz1,r2,zz2):
        femm.mi_addnode(r1,zz1); femm.mi_addnode(r2,zz1)
        femm.mi_addnode(r2,zz2); femm.mi_addnode(r1,zz2)
        femm.mi_addsegment(r1,zz1,r2,zz1); femm.mi_addsegment(r2,zz1,r2,zz2)
        femm.mi_addsegment(r2,zz2,r1,zz2); femm.mi_addsegment(r1,zz2,r1,zz1)
    rect(0,z0,r_so,z1); rect(0,z1,r_m,z2); rect(0,z2,r_m,z3); rect(r_si,z1,r_so,z3)
    Rout=4*r_so; femm.mi_makeABC(7, Rout, 0, (z0+z3)/2, 0)
    def blk(r,zz,mat,md,grp):
        femm.mi_addblocklabel(r,zz); femm.mi_selectlabel(r,zz)
        femm.mi_setblockprop(mat,1,0,"<None>",md,grp,0); femm.mi_clearselected()
    blk(r_m*0.5,(z0+z1)/2,"Yoke",0,1); blk(r_m*0.5,(z1+z2)/2,"NdFeB",90,2)
    blk(r_m*0.5,(z2+z3)/2,"Yoke",0,1); blk((r_si+r_so)/2,(z1+z3)/2,"Yoke",0,1)
    blk(Rout*0.8,(z0+z3)/2,"Air",0,0)
    femm.mi_zoomnatural(); femm.mi_saveas(TMP_FEM); femm.mi_analyze(1); femm.mi_loadsolution()
    r_mid=r_m+t_g/2; N=21; s=0.0
    for i in range(N):
        zz=z2+(z3-z2)*(i+0.5)/N
        s+=abs(femm.mo_getpointvalues(r_mid, zz)[1])
    return s/N   # Tesla

def metrics(r_m,h_m,h_y1,t_g,t_y3=0.8):
    Bg = compute_Bg(r_m,h_m,h_y1,t_g,t_y3)          # T
    rm,hm,hy,tg = r_m*1e-3,h_m*1e-3,h_y1*1e-3,t_g*1e-3
    dcs=RP*DC; w_c=tg-2*CC; h_c=hy
    n_half=(h_c/dcs)*((2/math.sqrt(3))*(max(w_c,0)/dcs))
    L_m=2*math.pi*(rm+CC+TB+w_c/2)
    n_eff=min(NSYM*(hy+0.08*hm)/h_c*n_half, NSYM*n_half)
    Kf=n_eff*Bg*L_m; mc=8960*(NSYM*n_half)*(math.pi/4)*DC**2*L_m
    A_g=2*math.pi*(rm+tg/2)*hy; By1=Bg*A_g/(math.pi*rm**2)
    return Bg*1e3, Kf, mc*1e3, Kf/mc, By1   # mT,N/A,g,(N/A/kg),T

def main():
    femm.openfemm(1)
    print("=== FEMM 검증: 해석모델 '최적' vs 기준 ===")
    print(f"{'case':<22}{'r_m':>5}{'h_y1':>5}{'t_g':>5}{'Bg[mT]':>9}{'Kf':>7}{'mc[g]':>7}{'Kf/mc':>8}{'By1':>6}")
    print("-"*74)
    for name,(rm,hm,hy,tg,ty3) in {
        "기준(현설계)":   (6.5,25,5.8,2.5,1.0),
        "해석최적(2T주장)":(8.0,100,2.0,1.2,0.8),
    }.items():
        try:
            Bg,Kf,mc,kpm,By1=metrics(rm,hm,hy,tg,ty3)
            print(f"{name:<22}{rm:>5.1f}{hy:>5.1f}{tg:>5.1f}{Bg:>9.0f}{Kf:>7.2f}{mc:>7.2f}{kpm:>8.0f}{By1:>6.2f}")
        except Exception as e:
            print(f"{name:<22}{rm:>5.1f}{hy:>5.1f}{tg:>5.1f}   FEMM 실패: {str(e)[:30]}")

    print("\n=== 형상 DOE (N52, t_g=1.2, h_m=40, t_y3=0.8) : FEMM 기준 max ===")
    print(f"{'r_m':>5}{'h_y1':>6}{'Bg[mT]':>9}{'Kf':>7}{'mc[g]':>7}{'Kf/mc':>9}{'By1':>6}")
    print("-"*49)
    best=None
    for rm in [6,7,8]:
        for hy in [2,3,4,5]:
            try:
                Bg,Kf,mc,kpm,By1=metrics(rm,40,hy,1.2,0.8)
            except Exception as e:
                print(f"{rm:>5.1f}{hy:>6.1f}   실패: {str(e)[:25]}"); continue
            tag=""
            if best is None or kpm>best[0]: best=(kpm,rm,hy,Bg,Kf,mc,By1); tag=" *"
            print(f"{rm:>5.1f}{hy:>6.1f}{Bg:>9.0f}{Kf:>7.2f}{mc:>7.2f}{kpm:>9.0f}{By1:>6.2f}{tag}")
    femm.closefemm()
    k=best
    print(f"\nFEMM 기준 best K_f/m_coil: r_m={k[1]}, h_y1={k[2]} -> Bg={k[3]:.0f}mT, Kf={k[4]:.2f}N/A, mc={k[5]:.2f}g, Kf/mc={k[0]:.0f}, By1={k[6]:.2f}T")

if __name__=="__main__":
    main()
