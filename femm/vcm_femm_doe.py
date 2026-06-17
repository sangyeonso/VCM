"""
vcm_femm_doe.py - N52 VCM B_g surrogate 재적합용 FEMM DOE
  격자: r_m x h_y1 x t_g (Ø20 제약 r_m+t_g+t_y3<=10 만족분만), h_m=40,t_y3=0.8,h_y2=4 고정
  출력: 2차 다항 B_g[mT](r_m,h_y1,t_g) 계수 + R^2 + CSV
"""
import os, math, sys, csv
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import femm, numpy as np

TMP_FEM=r"C:\femm_vcm\vcm_doe.fem"; os.makedirs(r"C:\femm_vcm", exist_ok=True)
OUT=os.path.join(os.path.dirname(__file__),"bg_surrogate_n52.csv")
MU0=4*math.pi*1e-7; BR=1.43; MU_YOKE=1400

def compute_Bg(r_m,h_m,h_y1,t_g,t_y3=0.8,h_y2=4.0):
    z0,z1=0.0,h_y2; z2=z1+h_m; z3=z2+h_y1; r_si=r_m+t_g; r_so=r_m+t_g+t_y3
    femm.newdocument(0); femm.mi_probdef(0,"millimeters","axi",1e-8,0,30)
    femm.mi_getmaterial("Air"); Hc=BR/(MU0*1.05)
    femm.mi_addmaterial("NdFeB",1.05,1.05,Hc,0,0,0,0,0,0,0,0)
    femm.mi_addmaterial("Yoke",MU_YOKE,MU_YOKE,0,0,0,0,0,0,0,0,0)
    def rect(r1,zz1,r2,zz2):
        femm.mi_addnode(r1,zz1); femm.mi_addnode(r2,zz1); femm.mi_addnode(r2,zz2); femm.mi_addnode(r1,zz2)
        femm.mi_addsegment(r1,zz1,r2,zz1); femm.mi_addsegment(r2,zz1,r2,zz2)
        femm.mi_addsegment(r2,zz2,r1,zz2); femm.mi_addsegment(r1,zz2,r1,zz1)
    rect(0,z0,r_so,z1); rect(0,z1,r_m,z2); rect(0,z2,r_m,z3); rect(r_si,z1,r_so,z3)
    Rout=4*r_so; femm.mi_makeABC(7,Rout,0,(z0+z3)/2,0)
    def blk(r,zz,mat,md,grp):
        femm.mi_addblocklabel(r,zz); femm.mi_selectlabel(r,zz)
        femm.mi_setblockprop(mat,1,0,"<None>",md,grp,0); femm.mi_clearselected()
    blk(r_m*0.5,(z0+z1)/2,"Yoke",0,1); blk(r_m*0.5,(z1+z2)/2,"NdFeB",90,2)
    blk(r_m*0.5,(z2+z3)/2,"Yoke",0,1); blk((r_si+r_so)/2,(z1+z3)/2,"Yoke",0,1)
    blk(Rout*0.8,(z0+z3)/2,"Air",0,0)
    femm.mi_zoomnatural(); femm.mi_saveas(TMP_FEM); femm.mi_analyze(1); femm.mi_loadsolution()
    r_mid=r_m+t_g/2; N=21; s=0.0
    for i in range(N):
        zz=z2+(z3-z2)*(i+0.5)/N; s+=abs(femm.mo_getpointvalues(r_mid,zz)[1])
    return s/N*1000.0  # mT

def main():
    femm.openfemm(1)
    pts=[]
    for t_g in [1.2,1.6,2.0,2.5]:
        rmax=9.2-t_g
        for r_m in [5,6,7,8]:
            if r_m> rmax+1e-9: continue
            for h_y1 in [2,3.5,5,7]:
                try:
                    Bg=compute_Bg(r_m,40,h_y1,t_g)
                    pts.append((r_m,h_y1,t_g,Bg))
                except Exception as e:
                    print(f"  실패 r_m={r_m} h_y1={h_y1} t_g={t_g}: {str(e)[:25]}")
    femm.closefemm()
    print(f"성공 {len(pts)}점")
    A=np.array([[1,rm,hy,tg,rm*rm,hy*hy,tg*tg,rm*hy,rm*tg,hy*tg] for (rm,hy,tg,_) in pts])
    y=np.array([b for (_,_,_,b) in pts])
    c,_,_,_=np.linalg.lstsq(A,y,rcond=None)
    pred=A@c; ss_res=np.sum((y-pred)**2); ss_tot=np.sum((y-y.mean())**2); R2=1-ss_res/ss_tot
    print("\n2차 다항 B_g[mT] 계수 (1,rm,hy,tg,rm^2,hy^2,tg^2,rm*hy,rm*tg,hy*tg):")
    print("Cs=["+",".join(f"{v:.4f}" for v in c)+"];")
    print(f"R^2={R2:.4f},  B_g 범위 {y.min():.0f}~{y.max():.0f} mT")
    with open(OUT,"w",newline="",encoding="utf-8-sig") as f:
        w=csv.writer(f); w.writerow(["r_m","h_y1","t_g","Bg_mT"]); w.writerows(pts)
    print(f"저장: {OUT}")

if __name__=="__main__": main()
