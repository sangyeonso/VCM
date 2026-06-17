"""
bg_surrogate.py  -  FEMM 격자 스윕으로 B_g(r_m, h_y1, t_g) 2차 다항 surrogate 적합
실행: (GPSM)  python bg_surrogate.py
  - h_m=25mm 고정(최적해가 항상 max를 고름), t_y3=1, h_y2=3 고정.
  - 결과: ../matlab/bg_coeffs.txt (10계수, MATLAB이 읽음), bg_grid.csv(원자료)
"""
import os, csv, math, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import numpy as np
import femm

TMP = r"C:\femm_vcm\vcm_tmp.fem"; os.makedirs(r"C:\femm_vcm", exist_ok=True)
HERE = os.path.dirname(__file__)
MU0 = 4*math.pi*1e-7; BR = 1.17; MU_YOKE = 1400; H_M = 25.0

def compute_Bg(r_m, h_y1, t_g, t_y3=1.0, h_y2=3.0, h_m=H_M):
    z0,z1 = 0.0,h_y2; z2=z1+h_m; z3=z2+h_y1
    r_si=r_m+t_g; r_so=r_m+t_g+t_y3
    femm.newdocument(0)
    femm.mi_probdef(0,"millimeters","axi",1e-8,0,30)
    femm.mi_getmaterial("Air")
    femm.mi_addmaterial("NdFeB117",1.05,1.05,BR/(MU0*1.05),0,0,0,0,0,0,0,0)
    femm.mi_addmaterial("Yoke",MU_YOKE,MU_YOKE,0,0,0,0,0,0,0,0,0)
    def rect(r1,a,r2,b):
        femm.mi_addnode(r1,a);femm.mi_addnode(r2,a);femm.mi_addnode(r2,b);femm.mi_addnode(r1,b)
        femm.mi_addsegment(r1,a,r2,a);femm.mi_addsegment(r2,a,r2,b)
        femm.mi_addsegment(r2,b,r1,b);femm.mi_addsegment(r1,b,r1,a)
    rect(0,z0,r_so,z1); rect(0,z1,r_m,z2); rect(0,z2,r_m,z3); rect(r_si,z1,r_so,z3)
    # 공극 밀봉(공기-공기 내부경계 -> 항상 독립영역, B_g 불변)
    for rr,zz in [(r_m,z2),(r_si,z2),(r_m,z3),(r_si,z3)]: femm.mi_addnode(rr,zz)
    femm.mi_addsegment(r_m,z2,r_si,z2); femm.mi_addsegment(r_m,z3,r_si,z3)
    Rout=2.5*max(r_so,z3); femm.mi_makeABC(7,Rout,0,(z0+z3)/2,0)  # VCM 전체 포함
    def blk(r,z,mat,md,g):
        femm.mi_addblocklabel(r,z);femm.mi_selectlabel(r,z)
        femm.mi_setblockprop(mat,1,0,"<None>",md,g,0);femm.mi_clearselected()
    blk(r_m*0.5,(z0+z1)/2,"Yoke",0,1); blk(r_m*0.5,(z1+z2)/2,"NdFeB117",90,2)
    blk(r_m*0.5,(z2+z3)/2,"Yoke",0,1); blk((r_si+r_so)/2,(z1+z3)/2,"Yoke",0,1)
    blk((r_m+r_si)/2,(z2+z3)/2,"Air",0,0)   # 공극
    blk((r_m+r_si)/2,(z1+z2)/2,"Air",0,0)   # 자석 옆
    blk(Rout*0.8,(z0+z3)/2,"Air",0,0)       # 바깥
    femm.mi_zoomnatural(); femm.mi_saveas(TMP); femm.mi_analyze(1); femm.mi_loadsolution()
    r_mid=r_m+t_g/2; N=21; s=0.0
    for i in range(N):
        zz=z2+(z3-z2)*(i+0.5)/N
        s+=abs(femm.mo_getpointvalues(r_mid,zz)[1])
    return s/N*1000.0  # mT

def basis(rm,hy,tg):
    return [1, rm,hy,tg, rm*rm,hy*hy,tg*tg, rm*hy,rm*tg,hy*tg]

# ===== 격자 (최적화 탐색영역 커버) =====
RM  = [2.0, 3.5, 5.0, 6.5]
HY1 = [4.0, 8.0, 12.0, 16.0]
TG  = [1.0, 2.0, 3.0]

def main():
    femm.openfemm(1)
    rows=[]; X=[]; y=[]
    print(f"격자 {len(RM)}x{len(HY1)}x{len(TG)} = {len(RM)*len(HY1)*len(TG)}점 FEMM 해석...")
    for rm in RM:
        for hy in HY1:
            for tg in TG:
                Bg=compute_Bg(rm,hy,tg)
                rows.append([rm,hy,tg,round(Bg,2)]); X.append(basis(rm,hy,tg)); y.append(Bg)
    femm.closefemm()
    X=np.array(X); y=np.array(y)
    c,_,_,_=np.linalg.lstsq(X,y,rcond=None)
    pred=X@c; err=pred-y
    R2=1-np.sum(err**2)/np.sum((y-y.mean())**2)
    print(f"\n2차 다항 surrogate 적합: R^2={R2:.4f}, 최대오차={np.max(np.abs(err)):.1f} mT, RMSE={np.sqrt(np.mean(err**2)):.1f} mT")
    print("계수(c0..c9):", np.array2string(c,precision=4,separator=', '))
    # 저장
    np.savetxt(os.path.join(HERE,"..","matlab","bg_coeffs.txt"), c.reshape(1,-1), fmt="%.8e")
    with open(os.path.join(HERE,"bg_grid.csv"),"w",newline="",encoding="utf-8-sig") as f:
        w=csv.writer(f); w.writerow(["r_m","h_y1","t_g","Bg_mT"]); w.writerows(rows)
    print(f"저장: matlab/bg_coeffs.txt (10계수), femm/bg_grid.csv ({len(rows)}점)")
    # 검증점 비교
    for rm,hy,tg,name in [(3,9,2.5,"paper"),(6.23,5.1,2.77,"opt"),(6.04,7.5,2.96,"final")]:
        bs=float(np.array(basis(rm,hy,tg))@c)
        print(f"  surrogate {name}: B_g={bs:.1f} mT")

if __name__=="__main__":
    main()
