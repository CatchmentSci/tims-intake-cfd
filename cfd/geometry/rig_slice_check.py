"""Check of the rig slice reading (Results, calibration rig; the 0.585 m/s reading).

The bore-averaged velocity of the symmetry-plane slice in the glass exceeds the continuity value
(R U / 3.0625) by 0.4 / 0.6 / 1.0 / 1.7 % at 0.081 / 0.154 / 0.3 / 0.585 m/s. Three tests:
1. glass: face-weighted and trapezoidal radial quadrature, centreline over |z| < 0.4 and < 0.2 mm;
2. glass: top and bottom halves separately, and the peak height (is the flow axisymmetric?);
3. sampler bore (4 mm, x -0.08..0), where the flow equals R U exactly: same slice bulk, and the
   centreline / Poiseuille fit.
Run from cfd/: python geometry/rig_slice_check.py
"""
import sys; sys.path.insert(0,"geometry")
import numpy as np, rig_reading as rr
for vel in rr.VELS:
    U=float(vel); R,_=rr.run_R(rr.RIG,vel); ug=R*U/rr.AREA_RATIO
    z=np.load(rr.SLICES/f"rig_U{vel}_field.npz"); off=z["off"]; n=np.diff(off)
    x=np.add.reduceat(z["xz"][:,0],off[:-1])/n; zz=np.add.reduceat(z["xz"][:,1],off[:-1])/n
    a,u=z["area"],z["u"]; names=[str(q) for q in z["names"]]
    plane=np.array(["flume" not in q for q in names])[z["zone"]]
    s=x-rr.XS["glass_mouth"]; RG=rr.RG
    ratios_face=[];ratios_trap=[];nf=[];cen=[];cen2=[]
    for s0 in np.arange(0.05, rr.PASSAGE-0.02, 0.005):
        m=plane&(s>=s0)&(s<s0+0.005)&(np.abs(zz)<RG)
        if m.sum()<5: continue
        w=np.abs(zz[m])*a[m]; ratios_face.append(np.sum(u[m]*w)/np.sum(w)/ug)
        r=np.abs(zz[m]); o=np.argsort(r); r=r[o]; uu=u[m][o]
        # bin by radius (average faces at the same r), add wall u=0 and axis
        rb=np.linspace(0,RG,60); ub=np.interp(rb, r, uu, left=uu[0], right=0.0); ub[-1]=0
        ratios_trap.append(2*np.trapezoid(ub*rb,rb)/RG**2/ug)
        nf.append(len(np.unique(np.round(r,6))))
        c=m&(np.abs(zz)<0.0004); cen.append(np.mean(u[c])/ug if c.any() else np.nan)
        c2=m&(np.abs(zz)<0.0002); cen2.append(np.mean(u[c2])/ug if c2.any() else np.nan)
    print(vel, 'faces-weighted bulk %.4f'%np.mean(ratios_face), 'trapz bulk %.4f'%np.mean(ratios_trap),
          'radial points %d'%np.median(nf), 'uc(<0.4mm) %.3f uc(<0.2mm) %.3f'%(np.nanmean(cen),np.nanmean(cen2)), 'RG mm', RG*1000)

for vel in rr.VELS:
    U=float(vel); R,_=rr.run_R(rr.RIG,vel); ug=R*U/rr.AREA_RATIO
    z=np.load(rr.SLICES/f"rig_U{vel}_field.npz"); off=z["off"]; n=np.diff(off)
    x=np.add.reduceat(z["xz"][:,0],off[:-1])/n; zz=np.add.reduceat(z["xz"][:,1],off[:-1])/n
    a,u=z["area"],z["u"]; names=[str(q) for q in z["names"]]
    plane=np.array(["flume" not in q for q in names])[z["zone"]]
    s=x-rr.XS["glass_mouth"]; RG=rr.RG
    out=[]
    for s0 in (0.05,0.15,0.25,0.35,0.45):
        m=plane&(s>=s0)&(s<s0+0.01)&(np.abs(zz)<RG)
        up=m&(zz>0); dn=m&(zz<0)
        f=lambda k: np.sum(u[k]*np.abs(zz[k])*a[k])/np.sum(np.abs(zz[k])*a[k])/ug
        # peak location
        i=np.argmax(np.where(m,u,-1)); out.append(f"s{int(s0*1000)}: top {f(up):.3f} bot {f(dn):.3f} peak z {zz[i]*1000:+.2f}mm")
    print(vel, " | ".join(out))

for vel in rr.VELS:
    U=float(vel); R,_=rr.run_R(rr.RIG,vel); ub_true=R*U
    z=np.load(rr.SLICES/f"rig_U{vel}_field.npz"); off=z["off"]; n=np.diff(off)
    x=np.add.reduceat(z["xz"][:,0],off[:-1])/n; zz=np.add.reduceat(z["xz"][:,1],off[:-1])/n
    a,u=z["area"],z["u"]; names=[str(q) for q in z["names"]]
    plane=np.array(["flume" not in q for q in names])[z["zone"]]
    res=[]
    for x0 in np.arange(-0.08,0.0,0.01):
        m=plane&(x>=x0)&(x<x0+0.01)&(np.abs(zz)<0.002)
        w=np.abs(zz[m])*a[m]; res.append(np.sum(u[m]*w)/np.sum(w)/ub_true)
    print(vel, 'sampler bore slice bulk / R U: %.4f (range %.4f-%.4f)'%(np.mean(res),min(res),max(res)))
print("centreline in sampler bore, u_c/(R U), and profile shape check")
for vel in rr.VELS:
    U=float(vel); R,_=rr.run_R(rr.RIG,vel); ub=R*U
    z=np.load(rr.SLICES/f"rig_U{vel}_field.npz"); off=z["off"]; n=np.diff(off)
    x=np.add.reduceat(z["xz"][:,0],off[:-1])/n; zz=np.add.reduceat(z["xz"][:,1],off[:-1])/n
    a,u=z["area"],z["u"]; names=[str(q) for q in z["names"]]
    plane=np.array(["flume" not in q for q in names])[z["zone"]]
    m=plane&(x>-0.05)&(x<-0.01)&(np.abs(zz)<0.002)
    c=m&(np.abs(zz)<0.0002)
    # fit u = uc (1 - (r/R)^2) to check Poiseuille, and integrate the fitted vs raw
    r=np.abs(zz[m]); uu=u[m]/ub
    A=np.c_[np.ones_like(r), r**2]; co=np.linalg.lstsq(A,uu,rcond=None)[0]
    wall_r = np.sqrt(-co[0]/co[1])
    # number of faces across radius and max r
    print(vel, 'uc %.4f'%(np.mean(u[c])/ub), 'parabola fit: uc %.4f, zero at r = %.4f mm (bore 2.000), max face |z| %.4f mm'%(co[0], wall_r*1000, r.max()*1000), 'n radial', len(np.unique(np.round(r,7))))
