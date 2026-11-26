"""Screening model of TIMS intake: a 1-D pipe network driven by the external pressure field.

The sampler is a series circuit. Water enters at the nose (pressure coefficient Cp_in), flows
through the inlet tube, expands into the body, contracts into the outlet and discharges into the
wake (Cp_base). The steady energy balance, written in inlet-tube velocity heads, is

    (Cp_in - Cp_base) * U^2 = u^2 * [K_entry + f_app*L/D + alpha + K_out*(A_i/A_o)^2 + ...]

and is solved for R = u/U. Laminar developing-flow friction follows Shah (1978); the exit
kinetic-energy coefficient alpha rises from 1 to 2 as the laminar profile develops. Turbulent
friction follows Blasius. This model screens the CFD design space; it does not replace it.

Shah RK. 1978. A correlation for laminar hydrodynamic entry length solutions for circular and
noncircular ducts. J. Fluids Eng. 100: 177-179.
"""
import math
import numpy as np
from scipy.optimize import brentq

NU = 1.0e-6        # m2/s, water at 20 C
RHO_S, RHO_W, MU, G = 2650.0, 998.0, 1.0e-3, 9.81


def tube_loss(u, D, L):
    """Return (friction + developing-flow loss, exit KE coefficient) in velocity heads."""
    Re = max(u * D / NU, 1e-9)
    if Re < 2300:
        xp = L / (D * Re)
        fapp_re = (3.44 / math.sqrt(xp)
                   + (1.25 / (4 * xp) + 16 - 3.44 / math.sqrt(xp)) / (1 + 2.1e-4 / xp**2))
        alpha = 1 + (1 - math.exp(-xp / 0.03))
        return 4 * fapp_re / Re * L / D, alpha
    return 0.316 * Re**-0.25 * L / D, 1.05


def solve_R(U, d, dCp=1.2):
    """Solve for the inlet velocity ratio R = u/U for sampler geometry dict d."""
    Ai = math.pi * d["Di"]**2 / 4

    def residual(R):
        u = R * U
        fi, ai = tube_loss(u, d["Di"], d["Li"])
        k = d.get("K_entry", 0.5) + fi + ai                       # inlet tube -> body (Borda-Carnot)
        if d.get("glass") is not None:                            # Phillips calibration add-on tube
            Dg, Lg = d["glass"]
            ug = u * (d["Di"] / Dg)**2
            fg, _ = tube_loss(ug, Dg, Lg)
            k = 0.5 * (ug / u)**2 + fg * (ug / u)**2 + 0.4 + fi + ai
        Ao = math.pi * d["Do"]**2 / 4
        uo = u * Ai / Ao
        if d["outlet"] == "tube":
            fo, ao = tube_loss(uo, d["Do"], d["Lo"])
            k += (0.8 + fo + ao) * (uo / u)**2                    # re-entrant entry, friction, exit
        else:                                                     # thin-plate orifice, Cd ~ 0.61
            k += (1 / 0.61**2) * (uo / u)**2
        return R**2 * k - dCp

    return brentq(residual, 1e-6, 3.0)


SAMPLERS = {
    "Phillips2000": dict(Di=0.004, Li=0.170, Do=0.004, Lo=0.170, outlet="tube", Db=0.098, Lb=1.0),
    "Perks2014":    dict(Di=0.008, Li=0.170, Do=0.008, Lo=0.170, outlet="tube", Db=0.090, Lb=1.0),
    "McDonald2010": dict(Di=0.002, Li=0.080, Do=0.003, Lo=0.0, outlet="orifice", Db=0.0635, Lb=0.228),
}
# Assumed values: tube lengths for Perks2014 (Phillips 150 mm + 20 mm protrusion) and the
# McDonald2010 inlet (funnel length + 20 mm protrusion, ~80 mm) are not reported. Confirm them.


def hazen_cut(Q, d):
    """Ideal horizontal-flow settling cut size (Stokes, quartz) for overflow rate Q/(Lb*Db)."""
    v0 = Q / (d["Lb"] * d["Db"])
    return math.sqrt(18 * MU * v0 / ((RHO_S - RHO_W) * G)) * 1e6


if __name__ == "__main__":
    print("1. Phillips (2000) Fig. 2 calibration points (read from figure) vs model")
    pts = [(1.19, 0.23), (1.33, 0.62), (1.46, 0.81), (1.55, 1.12), (1.75, 1.40), (1.76, 1.48)]
    cal = dict(SAMPLERS["Phillips2000"], glass=(0.007, 0.600))
    print(" U cm/s | R published | R if x3.06/1.75 | model (calibration rig) | model (deployed)")
    for lx, ly in pts:
        U, u = 10**lx, 10**ly
        print(f"  {U:5.1f} |   {u/U:.3f}     |     {u/U*3.0625/1.75:.3f}       |"
              f"        {solve_R(U/100, cal):.3f}            |     {solve_R(U/100, SAMPLERS['Phillips2000']):.3f}")

    print("\n2. Sampler variants (dCp = 1.2; sensitivity to dCp 1.1-1.5 in brackets)")
    for name, d in SAMPLERS.items():
        V = math.pi * d["Db"]**2 / 4 * d["Lb"]
        print(f" {name}: area ratio {(d['Db']/d['Di'])**2:.0f}, body volume {V*1e3:.2f} L")
        for U in [0.1, 0.3, 0.5, 1.0, 1.5, 2.5]:
            R = solve_R(U, d)
            lo, hi = solve_R(U, d, 1.1), solve_R(U, d, 1.5)
            Q = R * U * math.pi * d["Di"]**2 / 4
            print(f"   U={U:3.1f}  R={R:.2f} [{lo:.2f}-{hi:.2f}]  Re_in={R*U*d['Di']/NU:6.0f}"
                  f"  Q={Q*6e4:6.3f} L/min  T_res={V/Q/60:7.1f} min  ideal cut={hazen_cut(Q, d):5.1f} um")

    print("\n3. Frictionless upper bound on R (K = entry + Borda exits only)")
    for dCp in [1.1, 1.2, 1.5]:
        print(f"   dCp={dCp}: R_max = {math.sqrt(dCp/(0.5+1+0.8+1)):.2f}")

    print("\n4. Density contrast of sediment-laden inflow vs gravity-current speed (h = 0.05 m)")
    for C in [50, 100, 500, 1000, 1700]:
        drho = C * 1e-3 * (1 - RHO_W / RHO_S)
        gp = G * drho / RHO_W
        print(f"   SSC={C:5d} mg/L: drho={drho:.3f} kg/m3, (g'h)^0.5={math.sqrt(gp*0.05)*1e3:5.1f} mm/s")
