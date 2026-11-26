"""Parametric TIMS geometry and named surface mesh for Fluent Meshing.

Builds the half-domain fluid volume (symmetry plane y = 0, gravity along -z lies in that plane;
--domain open_full builds the whole cross-section instead)
around a time-integrated mass-flux sampler, splits it into sizing regions, classifies every
surface into a named zone, triangulates it with gmsh, and writes:

  <case>.msh   Fluent ASCII boundary mesh with named face zones (input to Fluent Meshing)
  <case>_<zone>.stl  one ASCII STL per zone, only with --stl (fallback input)
  <case>.step  fluid volume (CAD fallback; carries no zone names)
  <case>_zones.json  zone list, areas, and the spec used

Coordinates are in metres. x is streamwise (flow towards +x) and z is vertical. x = 0 is the outer
face of the upstream end cap.

Usage:  python tims_geometry.py phillips2000 [--size-scale 1.0] [--domain open] [--out DIR]
"""
import argparse
import json
import math
from dataclasses import dataclass, asdict, field, replace
from pathlib import Path

import gmsh
import numpy as np


@dataclass
class SamplerSpec:
    name: str
    body_id: float          # internal diameter of the main body
    body_od: float          # external diameter of the main body
    body_len: float         # internal length of the settling chamber
    cap_t: float            # end-cap thickness
    cone_len: float         # nose-cone (funnel) length, base diameter = body_od
    tube_id: float          # inlet/outlet tube bore
    tube_od: float
    inlet_len: float        # total inlet tube length
    inlet_in: float         # inlet protrusion into the chamber
    outlet_len: float
    outlet_in: float
    tail_cone_len: float = 0.0   # 0 = blunt rear end cap
    # Outlet tube bore and o.d. where they differ from the inlet's (0 = tube_id / tube_od;
    # standard construction only, e.g. the Foets et al. 2020 2 mm outlet)
    outlet_id: float = 0.0
    outlet_od: float = 0.0
    # Optional parts (0 = absent, which gives the Phillips construction):
    collar_od: float = 0.0       # socketed end caps fitted over the body: collar o.d. ...
    collar_len: float = 0.0      # ... and length from the cap's outer face
    stem_d: float = 0.0          # funnel stem o.d.; cone runs from the cap to the stem
    stem_len: float = 0.0        # stem length beyond the cone (both funnels)
    nozzle_id: float = 0.0       # panel-mount insert through each cap: bore, ...
    nozzle_od: float = 0.0       # ... o.d. of the part inside the chamber, ...
    nozzle_len: float = 0.0      # ... and total length; the tubing bore ends where it starts
    # Elliott-type construction (elbow_c > 0): caps with half-elliptical crowns (height cap_t
    # beyond the pipe end) over a hub, a dished floor one cap_wall thinner, compression fittings
    # whose thread ends are flush with the floor, chamfered tube tips and an L-shaped outlet
    cap_wall: float = 0.0
    hub_len: float = 0.0         # cap socket over the pipe; hub o.d. is collar_od
    fit_d: float = 0.0           # fitting hex/nut, modelled round, and its length (inlet)
    fit_len: float = 0.0
    pass_d: float = 0.0          # through-bore of both fittings, between tube seat and floor
    seat: float = 0.0            # tube inserted this far past each fitting's outer end
    chamfer_half: float = 0.0    # tube-tip taper half-angle (deg) and the wall left at the tip
    land: float = 0.0
    elbow_c: float = 0.0         # elbow: thread end to port centreline, thread axis to nut face,
    elbow_arm: float = 0.0       #  body width; the tube points down (-z) to out_tip_z
    elbow_w: float = 0.0
    out_tip_z: float = 0.0
    # McDonald-type construction (hole_d > 0): hollow funnel over the tube end (no upstream cap;
    # its interior is part of the chamber, the inlet tube held in its spout), plug cap with an
    # off-axis outlet hole at height hole_z
    hole_d: float = 0.0
    hole_z: float = 0.0
    funnel_wall: float = 0.0
    funnel_half: float = 0.0     # cone half-angle (deg)
    spout_d: float = 0.0
    spout_len: float = 0.0
    # In-chamber in-line device on the inlet (dev_len > 0; Perks flow-meter variant): a tube of
    # tube_od / tube_id from the insert's inner end, dev_gap long, then the device's port (o.d.,
    # length), its body (a box dev_len x dev_w x dev_h with rounded axial edges dev_r, the flow
    # axis dev_top below its top face) and its outlet port, from which the flow leaves; bore dev_bore
    dev_gap: float = 0.0
    dev_port_od: float = 0.0
    dev_port_len: float = 0.0
    dev_bore: float = 0.0
    dev_len: float = 0.0
    dev_w: float = 0.0
    dev_h: float = 0.0
    dev_top: float = 0.0
    dev_r: float = 0.0
    # Calibration glass tube over the inlet (glass_len > 0; Phillips et al. 2000 rig, plain
    # construction only): bore, o.d. and length, and how far it is sleeved back over the inlet tube
    # from its tip. The sleeve's gap round the tube is filled (the seal), so the flow sees the glass
    # bore end in a flat step onto the inlet tube's tip
    glass_id: float = 0.0
    glass_od: float = 0.0
    glass_len: float = 0.0
    glass_overlap: float = 0.0
    sources: dict = field(default_factory=dict)


SPECS = {
    "phillips2000": SamplerSpec(
        name="phillips2000", body_id=0.098, body_od=0.110, body_len=1.000, cap_t=0.010,
        cone_len=0.100, tube_id=0.004, tube_od=0.006,
        inlet_len=0.150, inlet_in=0.020, outlet_len=0.150, outlet_in=0.020,
        sources={
            "body_id, body_len, tube_id, inlet/outlet_len, *_in": "Phillips et al. 2000, p. 2590-2591",
            "body_od": "ASSUMED: 110 mm OD UK PVC pipe; Smith & Owens 2014 Fig. 3 shows ~10 cm",
            "cap_t": "ASSUMED; Smith & Owens report 8.68 L internal volume vs 7.54 L for the tube alone",
            "cone_len": "ASSUMED: funnel base = body OD; length 100 mm, so the inlet tube protrudes "
                        "20 mm beyond the tip, consistent with 150 mm tube - 20 mm - cap - cone",
            "tube_od": "ASSUMED: standard 6 x 4 mm nylon pneumatic tubing",
            "tail": "blunt rear cap (Phillips et al. 2000 Fig. 1)",
        }),
    # Perks et al. (2014) design (the builders' parts list).
    # Each end: a CPC PLC42006 panel-mount insert through the cap, Nalgene 980 tubing pushed onto
    # its barb and running through the stem of a Nalgene 4256-0414 funnel that fits over the cap,
    # the tubing ending 20 mm beyond the stem tip. Funnel interior modelled solid (sealed).
    "perks2014": SamplerSpec(
        name="perks2014", body_id=0.0807, body_od=0.0889, body_len=1.000, cap_t=0.005,
        cone_len=0.080, tube_id=0.00794, tube_od=0.0127,
        inlet_len=0.1615, inlet_in=0.0195, outlet_len=0.1615, outlet_in=0.0195,
        tail_cone_len=0.080, collar_od=0.097, collar_len=0.040, stem_d=0.018, stem_len=0.037,
        nozzle_id=0.00635, nozzle_od=0.0165, nozzle_len=0.0485,
        sources={
            "body_id, body_od": "3in grey PVC Class C pipe: 88.9 mm OD x 4.1 mm wall (supplier listings)",
            "body_len": "Perks et al. 2014 p. 4799: length 1000 mm (taken as the internal length)",
            "tube_id, tube_od": "Nalgene 980 braided PVC 5/16 x 1/2 in (Fisher 10448612)",
            "cone_len, stem_d, stem_len": "Nalgene 4256-0414 / Fisher 10709041: 98 mm top, 18 mm x 37 mm "
                                          "stem, 117 mm high, so an 80 mm cone",
            "collar_od": "ASSUMED: 3in plain PVC solvent cap, ~4 mm wall over the 88.9 mm pipe; the "
                         "funnel rim (98 mm) is taken equal to it to avoid a 0.5 mm step",
            "collar_len, cap_t": "ASSUMED: 40 mm socket, 5 mm end",
            "nozzle_*": "CPC PLC42006 (RS 138-363): PLC series 1/4in flow, length 1.91 in, body 0.65 in "
                        "(CPC catalogue 2004 p. 21); ASSUMED centred so 19.5 mm projects into the chamber",
            "inlet_len, outlet_len": "tubing 20 mm beyond the 117 mm funnel (builder, ~2 cm) to the "
                                     "insert's inner end 19.5 mm inside the 5 mm cap",
            "tail": "second funnel as a tail cone, outlet tubing 20 mm beyond it (builder)",
        }),
    # Elliott et al. (2017) bi-directional TIMS, as built per Appendix A.3. x = 0 at the front
    # cap's outer apex.
    "elliott2017": SamplerSpec(
        name="elliott2017", body_id=0.10226, body_od=0.1143, body_len=1.000, cap_t=0.03572,
        cone_len=0.0, tube_id=0.0040, tube_od=0.009525,
        inlet_len=0.178418, inlet_in=-0.02970, outlet_len=0.028448, outlet_in=-0.02970,
        collar_od=0.12779, cap_wall=0.00602, hub_len=0.0508,
        fit_d=0.022225, fit_len=0.042418, pass_d=0.00635, seat=0.014,
        chamfer_half=22.5, land=0.0005,
        elbow_c=0.028448, elbow_arm=0.032512, elbow_w=0.020638, out_tip_z=-0.063894,
        sources={
            "body_id, body_od, body_len": "4in Sch 40 PVC pipe, 1 m (Elliott et al. 2017 A.3; ASTM D1785 "
                                          "4.500 x 0.237 in). Text's 98 mm is the Phillips figure",
            "collar_od, cap_t, hub_len, cap_wall": "Spears 447-040 sheet: M 5-1/32in, W 3-13/32in; socket 2.000in "
                                                  "(ASTM D2466); wall 0.237in. Crown W - socket, half-elliptical; "
                                                  "pipe end at the socket bottom ASSUMED",
            "tube_id, tube_od, chamfer": "3/8in nylon bored to 4 mm, 150 mm; 45 deg chamfer to 0.5 mm wall "
                                         "(A.3), read as 22.5 deg per side from photo A.2c",
            "fit_*, pass_d": "McMaster 50775K335: 1.67in long, 7/8in nut, 0.250in bore; thread end flush with "
                             "the cap floor (6 mm engaged = cap wall)",
            "seat": "ASSUMED 14 mm",
            "elbow_*": "McMaster 50775K349: 1.12in thread end to port centreline, 1.28in thread axis to nut "
                       "face (centre-to-end); body width ASSUMED 13/16in",
            "out_tip_z": "tip level with the cap bottom (Figs. 1b-c), so ~31 mm of tube exposed; "
                         "the text's 150 mm outflow tube would reach below the sampler",
            "omitted": "vents (5/16-18 screws, assumed closed), uprights, clamps, the paired sampler",
        }),
    # McDonald et al. (2010) trap.
    # x = 0 lies 10 mm upstream of the polycarbonate tube's upstream end (cav0 = cap_t = 10 mm), so
    # the plug cap fits the standard stations: tube 0.010-0.238, plug inner face cav1 = 0.228.
    "mcdonald2010": SamplerSpec(
        name="mcdonald2010", body_id=0.0635, body_od=0.06985, body_len=0.218, cap_t=0.010,
        cone_len=0.0, tube_id=0.002, tube_od=0.004,
        inlet_len=0.096162, inlet_in=0.020, outlet_len=0.010, outlet_in=0.0,
        collar_len=0.010, hole_d=0.003, hole_z=0.015875,
        funnel_wall=0.0015, funnel_half=30.0, spout_d=0.008, spout_len=0.020,
        sources={
            "body": "22.8 cm polycarbonate tube; 63.5 mm bore, matching the stated ~1000 area ratio (p. 227, "
                    "Fig. 2); wall ASSUMED 1/8in (2.5 x 2.75in tube)",
            "inlet": "2 mm ID rigid polystyrene tube extended 20 mm into the body (p. 227); o.d. 4 mm ASSUMED; "
                     "flush with the spout tip (Fig. 2)",
            "funnel": "fixed over the upstream end, hollow and open to the body (assumed); 60 deg cone, "
                      "1.5 mm wall, 8 x 20 mm spout, 10 mm collar over the tube: ASSUMED",
            "outlet": "plastic cap with a 3 mm hole near the top (p. 228); cap ASSUMED a 10 mm plug flush with "
                      "the tube end; hole centre 3/4 up the bore (Fig. 2: 0.71 and 0.75)",
            "omitted": "boom and pivot-arm mounts, 6 mm ice deflector rods",
        }),
}

# Perks et al. (2014) with a Cynergy3 UF25B ultrasonic flow meter inside the chamber on the
# inlet line: an additional variant, perks2014 itself
# unchanged. 20 mm of tubing force-fitted from the insert to the meter's inlet port; no tube after it.
SPECS["perks2014_flowmeter"] = replace(
    SPECS["perks2014"], name="perks2014_flowmeter",
    dev_gap=0.020, dev_port_od=0.01666, dev_port_len=0.014, dev_bore=0.010,
    dev_len=0.069, dev_w=0.038, dev_h=0.045, dev_top=0.0115, dev_r=0.004,
    sources={**SPECS["perks2014"].sources,
             "dev": "Cynergy3 UF25B datasheet (UF25B-2014): 97.0 overall, 69.0 x 38.0 x 45.0 body, 10 mm bore, "
                    "3/8in BSP male ports (14 mm, 16.66 major dia.); flow axis ~11.5 below the top (measured on the "
                    "drawing); builder: inside the sampler, 20 mm force-fitted tubing before it, none after",
             "dev ASSUMED": "on the sampler axis, cable up (cable not modelled); centred across its width (the "
                            "true axis is 17/21 mm across) so the half model is exact; 4 mm edge radius; "
                            "force-fit joints as butt joints"})


# Foets et al. (2020) Phillips sampler with a 2 mm outlet: an additional variant, phillips2000
# itself unchanged.
SPECS["foets2020"] = replace(
    SPECS["phillips2000"], name="foets2020", outlet_id=0.002, outlet_od=0.004,
    sources={**SPECS["phillips2000"].sources,
             "outlet_id": "Foets et al. 2020 p. 4711: Phillips design 'except that the diameter of the outlet "
                          "was 2 mm instead of 4 mm'",
             "outlet_od": "ASSUMED: standard 4 x 2 mm tubing (range 3-6 mm)",
             "outlet_len, outlet_in": "ASSUMED as Phillips (150 mm, 20 mm into the chamber): only the "
                                      "diameter is said to differ"})

# Phillips et al. (2000) calibration rig: the sampler with the glass tube over its inlet, run in the
# flume_phillips domain
SPECS["phillips2000_rig"] = replace(
    SPECS["phillips2000"], name="phillips2000_rig",
    glass_id=0.007, glass_od=0.010, glass_len=0.600, glass_overlap=0.021,
    sources={**SPECS["phillips2000"].sources,
             "glass_id, glass_len": "Phillips et al. 2000 p. 2592: glass tube 7 mm (ID) x 600 mm placed over "
                                    "the inlet tube, water-tight seal at the sampler end, open end facing the flow",
             "glass_od": "ASSUMED: 1.5 mm wall (common 7 mm-bore tubing; 1.0-1.5 mm);",
             "glass_overlap": "ASSUMED: sleeved over the 20 mm protrusion to the cone (glass end 1 mm into "
                              "the cone), the gap round the tube filled;"})
# Goharrokhi et al. (2019) sampler in their flume (flume_goharrokhi domain). Its build details are
# not reported; it is ASSUMED to be built as Phillips et al. (2000): a solid nose cone
# fitted directly to the body, base = pipe o.d., the inlet tube its only passage, running from the
# nozzle through the cone into the chamber; nominal 4-inch PVC body, 98 mm i.d. x 1000 mm; 4 mm
# outlet tube. Their Fig. 1 is not to scale and gives no dimensions.
SPECS["goharrokhi2019"] = replace(
    SPECS["phillips2000"], name="goharrokhi2019", body_od=0.1143,
    sources={**SPECS["phillips2000"].sources,
             "body_id, body_len": "Goharrokhi et al. 2019 p. 2049 (98 mm ID x 100 cm); ASSUMED nominal "
                                  "4-inch PVC",
             "body_od": "ASSUMED: nominal 4-inch North American PVC pipe is 4.5 in (114.3 mm) o.d.; "
                        "98 mm i.d. matches Schedule 80",
             "tube_id, inlet_len": "Goharrokhi et al. 2019 p. 2049 (4 mm x 150 mm); outlet 4 mm i.d. "
                                   "ASSUMED as the inlet",
             "nose": "ASSUMED as Phillips et al. (2000): solid cone, base = pipe o.d., the inlet "
                     "tube passing through it into the chamber",
             "cone_len": "ASSUMED: not given. 100 mm keeps the published 150 mm inlet tube with 20 mm "
                         "beyond the tip and 20 mm into the chamber, as the Phillips model; range to "
                         "300 mm (goharrokhi2019_longnose)"})
# Nose-length sensitivity: the cone lengthened so the tip is 320 mm ahead of the front cap (the
# ~2.9 body diameters read from their Fig. 2 photo); the inlet tube lengthens to
# keep 20 mm beyond the tip and 20 mm in the chamber. Upper end of the assumed cone-length range.
SPECS["goharrokhi2019_longnose"] = replace(
    SPECS["goharrokhi2019"], name="goharrokhi2019_longnose", cone_len=0.300, inlet_len=0.350,
    sources={**SPECS["goharrokhi2019"].sources,
             "cone_len, inlet_len": "ASSUMED sensitivity: 300 mm cone (tip 320 mm ahead of the cap, "
                                    "~2.9 D in their Fig. 2 photo), tube 350 mm"})


# Flush-tip nose: the 150 mm inlet tube, 20 mm in the chamber and through the
# 10 mm cap, leaves 120 mm for the cone, so a 120 mm cone ends flush with the tube tip, no protrusion.
# An SRWP photo of a deployed Phillips sampler shows a long conical nose with little or no visible
# tube. Sensitivity on the Phillips model (open domain) and on the Goharrokhi flume model.
SPECS["phillips2000_flushtip"] = replace(
    SPECS["phillips2000"], name="phillips2000_flushtip", cone_len=0.120,
    sources={**SPECS["phillips2000"].sources,
             "cone_len": "ASSUMED sensitivity: 120 mm, the inlet tube flush at the cone tip "
                         "(150 mm tube - 20 mm in the chamber - 10 mm cap)"})
SPECS["goharrokhi2019_flushtip"] = replace(
    SPECS["goharrokhi2019"], name="goharrokhi2019_flushtip", cone_len=0.120,
    sources={**SPECS["goharrokhi2019"].sources,
             "cone_len": "ASSUMED sensitivity: 120 mm cone from the 114.3 mm body to the 6 mm tube, "
                         "the tube flush at the tip; the author"})
# McDonald et al. (2010) with the outlet hole at the top of the end plug (their p. 228, "at the
# top"), a sensitivity test: hole centre
# 29.25 mm above the axis, so its top edge is 1 mm inside the 31.75 mm bore radius.
SPECS["mcdonald2010_holetop"] = replace(
    SPECS["mcdonald2010"], name="mcdonald2010_holetop", hole_z=0.02925,
    sources={**SPECS["mcdonald2010"].sources,
             "outlet": "SENSITIVITY (20 Nov review): 3 mm hole 'at the top' (p. 228), centre 29.25 mm "
                       "above the axis, top edge 1 mm inside the bore; production model 3/4 up (15.9 mm)"})

# Radii of the sizing regions (m); interior faces are named by which region they bound
R_JET, R_TIP, R_NEAR, R_WAKE = 0.012, 0.030, 0.110, 0.160


@dataclass
class DomainSpec:
    upstream: float = 1.0     # from the inlet tube tip to the domain inlet
    downstream: float = 3.0   # from the outlet tube tip to the domain outlet
    half_width: float = 0.55  # y extent (half model); z extent is +/- this value
    # Confined (flume) domain: side wall at y = half_width and floor at z = -half_width (currently
    # free-slip, see the classification below), rigid-lid symmetry at the water surface
    # z = +half_width. The unconfined domain has symmetry on all three.
    confined: bool = False
    name: str = ""
    # Full domain (no symmetry plane at y = 0), to test whether the half model suppresses
    # antisymmetric wake shedding; y then spans +/- half_width
    full: bool = False
    # Flume with a bed (Goharrokhi et al. 2019): floor this far below the sampler axis and
    # rigid lid this far above it (0 = +/- half_width). flume_walls: the side wall and floor are
    # no-slip walls. They are named symmetry_flumewall_* so Fluent Meshing grows no prism layers on
    # them (prisms there crashed it) and are retyped as walls in the solver, resolved for wall
    # functions by their face size (run_tims_case_v3.py FACE_SIZES "flumewall")
    floor: float = 0.0
    lid: float = 0.0
    flume_walls: bool = False


DOMAINS = {
    "open": DomainSpec(),
    # Phillips et al. (2000, p. 2592) calibration flume: Armfield 0.3 x 0.3 x 6.0 m. ASSUMED:
    # sampler centred in the cross-section, flume running full (0.3 m deep), 1 m upstream and
    # 3 m downstream of the tube tips (the 6 m flume less entry and tailgate lengths).
    "flume_phillips": DomainSpec(upstream=1.0, downstream=3.0, half_width=0.15, confined=True,
                                 name="flume"),
    "open_full": DomainSpec(full=True, name="full"),
    # Domain-extent test: every boundary
    # moved out, upstream 1.0 -> 1.5 m, downstream 3.0 -> 4.5 m, side and top/bottom 0.55 -> 1.10 m
    # (blockage 0.8 % -> 0.2 %). Half model, as production.
    "open_big": DomainSpec(upstream=1.5, downstream=4.5, half_width=1.10, name="big"),
    # Goharrokhi et al. (2019, p. 2052): 0.95 m wide flume,
    # 0.30 m deep, sampler on the centreline with the inlet tube 0.145 m above the bed. Half flume:
    # symmetry on the centreline, glass side wall at y = 0.475 m, sheet-metal bed, rigid lid. 1 m
    # upstream of the tip and 3 m downstream of the outlet tip (as the Phillips flume domain)
    "flume_goharrokhi": DomainSpec(upstream=1.0, downstream=3.0, half_width=0.475, confined=True,
                                   floor=0.145, lid=0.155, flume_walls=True, name="gflume"),
}


def x_stations(s: SamplerSpec):
    """Key streamwise coordinates."""
    x_cav0 = s.cap_t
    x_cav1 = s.cap_t + s.body_len
    x_rear = x_cav1 + s.cap_t
    x_in_inner = x_cav0 + s.inlet_in
    x_in_tip = x_in_inner - s.inlet_len
    x_out_inner = x_cav1 - s.outlet_in
    x_out_tip = x_out_inner + s.outlet_len
    tail_stem = s.stem_len if s.tail_cone_len > 0 else 0.0
    return dict(cav0=x_cav0, cav1=x_cav1, rear=x_rear, in_inner=x_in_inner, in_tip=x_in_tip,
                out_inner=x_out_inner, out_tip=x_out_tip, cone_tip=-s.cone_len - s.stem_len,
                tail_tip=x_rear + s.tail_cone_len + tail_stem,
                # where each tubing bore gives way to the narrower insert bore
                noz_in=x_in_inner - s.nozzle_len, noz_out=x_out_inner + s.nozzle_len,
                in_end=(x_in_inner + s.dev_gap + 2 * s.dev_port_len + s.dev_len) if s.dev_len else x_in_inner,
                # calibration glass tube: its sampler end and its open mouth
                glass_end=x_in_tip + s.glass_overlap,
                glass_mouth=x_in_tip + s.glass_overlap - s.glass_len,
                front=(x_in_tip + s.glass_overlap - s.glass_len) if s.glass_len else x_in_tip)



def elliott_geometry(s: SamplerSpec, xs: dict):
    """Extra stations of the Elliott-type construction (see SamplerSpec)."""
    g = dict(x_fit=xs["in_inner"] - s.fit_len)                 # inlet fitting's outer end
    g["x_seat"] = g["x_fit"] + s.seat                           # inlet tube seat
    g["xe"] = xs["out_tip"]                                     # elbow port centreline
    g["z_nut"] = -s.elbow_arm                                   # elbow nut face
    g["z_seat"] = g["z_nut"] + s.seat                           # outlet tube seat
    g["z_tip"] = s.out_tip_z
    g["chamfer"] = (s.tube_od / 2 - s.tube_id / 2 - s.land) / math.tan(math.radians(s.chamfer_half))
    return g


def elliott_solid(s: SamplerSpec, xs: dict, g: dict, occ):
    """Solid sampler of the Elliott-type construction: pipe, caps with half-elliptical crowns over
    hubs, dished floors, compression fittings, chamfered tubes and an L-shaped outlet (-z)."""
    rb_o, rb_i, rt_o, rt_i = s.body_od / 2, s.body_id / 2, s.tube_od / 2, s.tube_id / 2
    rc, rfb, c = s.collar_od / 2, s.pass_d / 2, g["chamfer"]
    cav0, cav1, xe = xs["cav0"], xs["cav1"], g["xe"]

    # Hollow body (crowns, hubs, pipe, dished floors, pipe ends inside the caps) as one revolved
    # section. Fusing dilated spheres to the hubs failed: the crowns meet the hubs tangentially and
    # the result could not be cut from the domain (box - solid returned the whole box).
    h_i, w, rear = s.cap_t - s.cap_wall, s.cap_wall, xs["rear"]
    pt = lambda x, z: occ.addPoint(x, 0, z)
    c0, c1 = pt(cav0, 0), pt(cav1, 0)
    P = [pt(0, 0), pt(cav0, rc), pt(cav0 + s.hub_len, rc), pt(cav0 + s.hub_len, rb_o),
         pt(cav1 - s.hub_len, rb_o), pt(cav1 - s.hub_len, rc), pt(cav1, rc), pt(rear, 0),
         pt(rear - w, 0), pt(cav1, rb_o), pt(cav1, rb_i), pt(cav0, rb_i), pt(cav0, rb_o), pt(w, 0)]
    curves = [occ.addEllipseArc(P[0], c0, P[1], P[1])]                  # front crown
    curves += [occ.addLine(P[i], P[i + 1]) for i in range(1, 6)]         # hubs and pipe outside
    curves += [occ.addEllipseArc(P[6], c1, P[6], P[7])]                  # rear crown
    curves += [occ.addLine(P[7], P[8]),                                  # rear end wall (axis)
               occ.addEllipseArc(P[8], c1, P[9], P[9])]                  # rear dished floor
    curves += [occ.addLine(P[i], P[i + 1]) for i in range(9, 12)]       # pipe ends and bore
    curves += [occ.addEllipseArc(P[12], c0, P[12], P[13]),               # front dished floor
               occ.addLine(P[13], P[0])]                                 # front end wall (axis)
    face = occ.addPlaneSurface([occ.addCurveLoop(curves)])
    rev = occ.revolve([(2, face)], 0, 0, 0, 1, 0, 0, 2 * math.pi)
    shell = [e for e in rev if e[0] == 3]
    occ.remove([(2, face), (0, c0), (0, c1)])            # the section itself would stay as a stray
    # inlet: fitting body (overlapping the crown), tube with a chamfered tip
    x_tip = xs["in_tip"]
    parts = shell + [
        (3, occ.addCylinder(g["x_fit"], 0, 0, 0.002 - g["x_fit"], 0, 0, s.fit_d / 2)),
        (3, occ.addCylinder(x_tip + c, 0, 0, g["x_fit"] + 0.001 - x_tip - c, 0, 0, rt_o)),
        (3, occ.addCone(x_tip, 0, 0, c, 0, 0, rt_i + s.land, rt_o)),
        # outlet: elbow body along the axis, nut and tube downwards, chamfered tip
        (3, occ.addCylinder(xs["rear"] - 0.002, 0, 0, xe + s.elbow_w / 2 - xs["rear"] + 0.002, 0, 0,
                            s.elbow_w / 2)),
        (3, occ.addCylinder(xe, 0, 0, 0, 0, g["z_nut"], s.fit_d / 2)),
        (3, occ.addCylinder(xe, 0, g["z_nut"] + 0.001, 0, 0, g["z_tip"] + c - g["z_nut"] - 0.001, rt_o)),
        (3, occ.addCone(xe, 0, g["z_tip"] + c, 0, 0, -c, rt_o, rt_i + s.land))]
    solid, _ = occ.fuse([parts[0]], parts[1:])
    # The narrower tube bores reach 0.1 mm into the wider fitting bores, so each seat step lies
    # exactly at the seat (an overlap the other way shifts the step off it)
    bores = [  # inlet: tube bore to the seat, then the fitting bore through the floor
        occ.addCylinder(x_tip - 1e-3, 0, 0, g["x_seat"] + 1e-4 - x_tip + 1e-3, 0, 0, rt_i),
        occ.addCylinder(g["x_seat"], 0, 0, xs["in_inner"] + 1e-3 - g["x_seat"], 0, 0, rfb),
        # outlet: fitting bore along the axis and down to the seat, then the tube bore
        occ.addCylinder(xs["out_inner"] - 1e-3, 0, 0, xe + rfb - xs["out_inner"] + 1e-3, 0, 0, rfb),
        occ.addCylinder(xe, 0, rfb, 0, 0, g["z_seat"] - rfb, rfb),
        # (the outlet tube runs in -z, so "into the fitting bore" is above the seat)
        occ.addCylinder(xe, 0, g["z_seat"] + 1e-4, 0, 0, g["z_tip"] - 1e-3 - g["z_seat"] - 1e-4, rt_i)]
    solid, _ = occ.cut(solid, [(3, b) for b in bores])
    return solid


def elliott_wall_zone(s: SamplerSpec, xs: dict, g: dict, box, tol=1e-5):
    """Zone of a sampler wall face in the Elliott-type construction, from its bounding box."""
    bx0, by0, bz0, bx1, by1, bz1 = box
    rmax = max(abs(by0), abs(by1), abs(bz0), abs(bz1))
    rt_o, rt_i, rfb, rb_o, rb_i = s.tube_od / 2, s.tube_id / 2, s.pass_d / 2, s.body_od / 2, s.body_id / 2
    xe = g["xe"]
    hx, hy = max(abs(bx0 - xe), abs(bx1 - xe)), max(abs(by0), abs(by1))
    # outlet: fitting bore down to the seat (incl. the seat step), vertical tube and its bore
    if hx <= rfb + tol and hy <= rfb + tol and bz0 >= g["z_seat"] - tol and bz1 <= rfb + tol:
        return "wall_downstream_bore"
    if hx <= rt_o + tol and hy <= rt_o + tol and bz0 >= g["z_tip"] - tol and bz1 <= g["z_seat"] + tol:
        return "wall_downstream_bore" if (hx <= rt_i + tol and hy <= rt_i + tol) else "wall_downstream_tube"
    if rmax <= rfb + tol and bx0 >= xs["out_inner"] - 1e-3 and bx1 <= xe + rfb + tol:
        return "wall_downstream_bore"
    if bx0 >= xs["rear"] - 2e-3 - tol:
        return "wall_fitting"                                    # elbow body and nut
    # inlet
    if rmax <= rt_i + tol and bx0 >= xs["in_tip"] - tol and bx1 <= g["x_seat"] + tol:
        return "wall_upstream_bore"
    if rmax <= rfb + tol and bx0 >= g["x_seat"] - tol and bx1 <= xs["in_inner"] + 1e-3:
        return "wall_upstream_bore"
    if rmax <= rt_o + tol and bx0 >= xs["in_tip"] - tol and bx1 <= g["x_fit"] + tol:
        return "wall_upstream_tube"
    if rmax <= s.fit_d / 2 + tol and bx0 >= g["x_fit"] - tol and bx1 <= 2e-3 + tol:
        return "wall_fitting"
    # chamber: pipe bore, the dished floors and the pipe ends inside the caps
    if rmax <= rb_i + tol and bx0 >= xs["cav0"] - tol and bx1 <= xs["cav1"] + tol:
        return "wall_cavity"
    if rmax <= rb_o + tol and ((bx1 <= xs["cav0"] + tol and bx0 >= xs["in_inner"] - 1e-3) or
                               (bx0 >= xs["cav1"] - tol and bx1 <= xs["out_inner"] + 1e-3)):
        return "wall_cavity"
    return "wall_body"                                           # pipe outside, hubs, crowns



def mcdonald_geometry(s: SamplerSpec, xs: dict):
    """Extra stations of the McDonald-type construction (see SamplerSpec)."""
    t30 = math.tan(math.radians(s.funnel_half))
    rf, ro = s.body_od / 2 + s.funnel_wall, s.body_od / 2
    g = dict(x_spout=xs["cav0"] - (rf - s.spout_d / 2) / t30,          # cone meets the spout
             x_ic=xs["cav0"] - (ro - s.tube_od / 2) / t30,              # inner cone meets the tube
             rf=rf, ov=0.0005)                                          # overlap of fused parts
    g["x_tip"] = g["x_spout"] - s.spout_len
    if abs(g["x_tip"] - xs["in_tip"]) > 1e-6:
        raise ValueError(f"inlet_len gives the tube tip at {xs['in_tip']:.6f}, the spout tip is at {g['x_tip']:.6f}")
    g["xe"], g["z_tip"] = xs["rear"], s.hole_z                         # outlet "tip": the hole exit
    return g


def mcdonald_solid(s: SamplerSpec, xs: dict, g: dict, occ):
    """Solid sampler of the McDonald-type construction, built from revolved sections: hollow funnel
    with collar and spout, inlet tube through the spout, body tube, plug cap with an outlet hole."""
    ro, ri, rf, sr = s.body_od / 2, s.body_id / 2, g["rf"], s.spout_d / 2
    to, ti, ov = s.tube_od / 2, s.tube_id / 2, g["ov"]
    x0, x1, rear, xt = xs["cav0"], xs["cav1"], xs["rear"], xs["in_tip"]

    def revolved(pts):
        tags = [occ.addPoint(x, 0, z) for x, z in pts]
        lines = [occ.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
        face = occ.addPlaneSurface([occ.addCurveLoop(lines)])
        vol = [e for e in occ.revolve([(2, face)], 0, 0, 0, 1, 0, 0, 2 * math.pi) if e[0] == 3]
        occ.remove([(2, face)])
        return vol

    parts = revolved([(x0 + s.collar_len, ro - ov), (x0 + s.collar_len, rf), (x0, rf), (g["x_spout"], sr),
                      (xt, sr), (xt, to - ov), (g["x_ic"], to - ov), (x0, ro - ov)])          # funnel
    parts += revolved([(xt, ti), (xs["in_inner"], ti), (xs["in_inner"], to), (xt, to)])     # inlet tube
    parts += revolved([(x0, ri), (rear, ri), (rear, ro), (x0, ro)])                          # body tube
    parts += revolved([(x1, 0), (rear, 0), (rear, ri + ov), (x1, ri + ov)])                  # plug cap
    solid, _ = occ.fuse([parts[0]], parts[1:])
    hole = occ.addCylinder(x1 - 1e-3, 0, s.hole_z, rear - x1 + 2e-3, 0, 0, s.hole_d / 2)
    solid, _ = occ.cut(solid, [(3, hole)])
    return solid


def mcdonald_wall_zone(s: SamplerSpec, xs: dict, g: dict, box, tol=1e-5):
    """Zone of a sampler wall face in the McDonald-type construction, from its bounding box."""
    bx0, by0, bz0, bx1, by1, bz1 = box
    rmax = max(abs(by0), abs(by1), abs(bz0), abs(bz1))
    ro, ri, hr = s.body_od / 2, s.body_id / 2, s.hole_d / 2
    flat_x = abs(bx1 - bx0) < tol
    if (bx0 >= xs["cav1"] - tol and bx1 <= xs["rear"] + tol and max(abs(by0), abs(by1)) <= hr + tol
            and abs(bz0 - s.hole_z) <= hr + tol and abs(bz1 - s.hole_z) <= hr + tol):
        return "wall_downstream_bore"                            # the outlet hole
    if rmax <= s.tube_id / 2 + tol and bx0 >= xs["in_tip"] - tol and bx1 <= xs["in_inner"] + tol:
        return "wall_upstream_bore"
    if flat_x and abs(bx0 - xs["in_tip"]) < tol and rmax <= s.spout_d / 2 + tol:
        return "wall_upstream_tube"                              # tube and spout tip, one annulus
    if rmax <= s.tube_od / 2 + tol and bx0 >= xs["in_tip"] - tol and bx1 <= xs["in_inner"] + tol:
        return "wall_upstream_tube"
    if flat_x and abs(bx0 - xs["rear"]) < tol:
        return "wall_rear"
    if (rmax <= ri + tol and bx0 >= xs["cav0"] - tol and bx1 <= xs["cav1"] + tol) or \
            (flat_x and abs(bx0 - xs["cav0"]) < tol and rmax <= ro + tol):
        return "wall_cavity"                                     # bore, plug face, tube end in the funnel
    if bx1 <= xs["cav0"] + s.collar_len + tol:
        return "wall_cone"                                       # funnel: collar, cone, spout, inside
    return "wall_body"


def build(s: SamplerSpec, dom: DomainSpec, out_dir: Path, size_scale: float = 1.0, gui=False,
          write_stl=False):
    xs = x_stations(s)
    if xs["in_tip"] > xs["cone_tip"] + 1e-9:
        raise ValueError("inlet tube does not reach beyond the nose cone tip")
    rb_o, rb_i, rt_o, rt_i = s.body_od / 2, s.body_id / 2, s.tube_od / 2, s.tube_id / 2
    ro_o, ro_i = (s.outlet_od or s.tube_od) / 2, (s.outlet_id or s.tube_id) / 2   # outlet tube
    if (s.outlet_id or s.outlet_od) and (s.elbow_c or s.hole_d or s.nozzle_len or s.tail_cone_len):
        raise ValueError("a separate outlet bore is built for the plain coupled construction only")
    rg_i, rg_o = s.glass_id / 2, s.glass_od / 2                                   # glass tube
    if s.glass_len and (s.elbow_c or s.hole_d or s.nozzle_len or s.dev_len
                        or s.stem_len or rg_i <= rt_o or xs["glass_end"] <= xs["cone_tip"]):
        raise ValueError("the glass tube is built over a plain stemless inlet tube only, wider than "
                         "the tube and reaching into the cone")
    # the fetch is measured from whatever meets the flow first: the glass mouth or the tube tip
    x0, x1 = xs["front"] - dom.upstream, xs["out_tip"] + dom.downstream
    H = dom.half_width
    z_lo, z_hi = (-dom.floor, dom.lid) if dom.floor else (-H, H)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    case_name = s.name
    if dom.name:
        case_name += f"_{dom.name}"
    gmsh.model.add(case_name)
    occ = gmsh.model.occ

    # --- solid sampler ------------------------------------------------------------------
    if s.elbow_c:                                             # Elliott-type construction
        # The crowns and floors are dilated spheres (B-spline surfaces), whose OCC bounding boxes
        # can be twice their true size; the zone classification needs tight boxes
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
        g = elliott_geometry(s, xs)
        solid = elliott_solid(s, xs, g, occ)
    elif s.hole_d:                                            # McDonald-type construction
        gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
        g = mcdonald_geometry(s, xs)
        solid = mcdonald_solid(s, xs, g, occ)
    else:
        body = occ.addCylinder(0, 0, 0, xs["rear"], 0, 0, rb_o)
        if s.collar_od:          # socketed caps over the body, fused before the cavity is cut out
            collars = [(3, occ.addCylinder(0, 0, 0, s.collar_len, 0, 0, s.collar_od / 2)),
                       (3, occ.addCylinder(xs["rear"] - s.collar_len, 0, 0, s.collar_len, 0, 0,
                                           s.collar_od / 2))]
            (_, body), *_ = occ.fuse([(3, body)], collars)[0]
        cavity = occ.addCylinder(xs["cav0"], 0, 0, s.body_len, 0, 0, rb_i)
        shell, _ = occ.cut([(3, body)], [(3, cavity)])            # hollow body first, so the
        # Nose (and tail) cone from the cap's outer diameter to the stem, or to the tube if no stem
        r_base = s.collar_od / 2 if s.collar_od else rb_o
        r_tip = s.stem_d / 2 if s.stem_d else rt_o
        cone = occ.addCone(0, 0, 0, -s.cone_len, 0, 0, r_base, r_tip)  # tube protrusions survive
        parts = shell + [(3, cone)]
        if s.stem_len:
            parts.append((3, occ.addCylinder(-s.cone_len, 0, 0, -s.stem_len, 0, 0, r_tip)))
        if s.tail_cone_len > 0:
            parts.append((3, occ.addCone(xs["rear"], 0, 0, s.tail_cone_len, 0, 0, r_base, r_tip)))
            if s.stem_len:
                parts.append((3, occ.addCylinder(xs["rear"] + s.tail_cone_len, 0, 0, s.stem_len, 0, 0, r_tip)))
        tube_in = occ.addCylinder(xs["in_tip"], 0, 0, s.inlet_len, 0, 0, rt_o)
        tube_out = occ.addCylinder(xs["out_inner"], 0, 0, s.outlet_len, 0, 0, ro_o)
        parts += [(3, tube_in), (3, tube_out)]
        if s.glass_len:              # solid glass to the cone; its bore is cut with the others
            parts.append((3, occ.addCylinder(xs["glass_mouth"], 0, 0, s.glass_len, 0, 0, rg_o)))
        if s.nozzle_len:                                          # insert bodies inside the chamber
            r_n = s.nozzle_od / 2
            parts.append((3, occ.addCylinder(xs["noz_in"], 0, 0, s.nozzle_len, 0, 0, r_n)))
            parts.append((3, occ.addCylinder(xs["out_inner"], 0, 0, s.nozzle_len, 0, 0, r_n)))
        if s.dev_len:                                             # in-line device in the chamber
            p0 = xs["in_inner"] + s.dev_gap
            m0, m1 = p0 + s.dev_port_len, p0 + s.dev_port_len + s.dev_len
            parts.append((3, occ.addCylinder(xs["in_inner"] - 1e-4, 0, 0, s.dev_gap + 2e-4, 0, 0, rt_o)))
            parts.append((3, occ.addCylinder(p0, 0, 0, s.dev_port_len + 1e-4, 0, 0, s.dev_port_od / 2)))
            parts.append((3, occ.addCylinder(m1 - 1e-4, 0, 0, s.dev_port_len + 1e-4, 0, 0, s.dev_port_od / 2)))
            box = occ.addBox(m0, -s.dev_w / 2, s.dev_top - s.dev_h, s.dev_len, s.dev_w, s.dev_h)
            occ.synchronize()
            edges = []                                            # round the four edges along x
            for _, f in gmsh.model.getBoundary([(3, box)], oriented=False):
                for _, c in gmsh.model.getBoundary([(2, f)], oriented=False):
                    b = gmsh.model.getBoundingBox(1, abs(c))
                    if b[3] - b[0] > 0.9 * s.dev_len and abs(c) not in edges:
                        edges.append(abs(c))
            parts.append(occ.fillet([box], edges, [s.dev_r])[0])
        solid, _ = occ.fuse([parts[0]], parts[1:])
        if s.nozzle_len:                                        # tubing bore, then the insert bore
            r_nb = s.nozzle_id / 2
            if s.dev_len:     # device: tube bore from the insert, then the device bore out to the chamber
                p0 = xs["in_inner"] + s.dev_gap
                extra = [occ.addCylinder(xs["in_inner"] - 1e-4, 0, 0, s.dev_gap + 2e-4, 0, 0, rt_i),
                         occ.addCylinder(p0, 0, 0, xs["in_end"] - p0 + 1e-3, 0, 0, s.dev_bore / 2)]
            else:
                extra = []
            bores = [occ.addCylinder(xs["in_tip"] - 1e-3, 0, 0, xs["noz_in"] - xs["in_tip"] + 1e-3, 0, 0, rt_i),
                     occ.addCylinder(xs["noz_in"] - 1e-4, 0, 0, s.nozzle_len + 1e-3, 0, 0, r_nb),
                     occ.addCylinder(xs["out_inner"] - 1e-3, 0, 0, s.nozzle_len + 1.1e-3, 0, 0, r_nb),
                     occ.addCylinder(xs["noz_out"], 0, 0, xs["out_tip"] - xs["noz_out"] + 1e-3, 0, 0, rt_i)] + extra
        else:
            bores = [occ.addCylinder(xs["in_tip"] - 1e-3, 0, 0, s.inlet_len + 2e-3, 0, 0, rt_i),
                     occ.addCylinder(xs["out_inner"] - 1e-3, 0, 0, s.outlet_len + 2e-3, 0, 0, ro_i)]
            if s.glass_len:          # glass bore from the mouth to the tube tip: the step lies there
                bores.append(occ.addCylinder(xs["glass_mouth"] - 1e-3, 0, 0,
                                             xs["in_tip"] - xs["glass_mouth"] + 1e-3, 0, 0, rg_i))
        solid, _ = occ.cut(solid, [(3, b) for b in bores])

    # --- half-domain fluid --------------------------------------------------------------
    y0 = -H if dom.full else 0.0
    box = occ.addBox(x0, y0, z_lo, x1 - x0, H - y0, z_hi - z_lo)
    fluid, _ = occ.cut([(3, box)], solid)

    # --- sizing regions and reporting faces (fragmented into the fluid) -----------------
    def half_cyl(xa, length, r):
        c = occ.addCylinder(xa, 0, 0, length, 0, 0, r)
        if dom.full:
            return [(3, c)]
        cut_box = occ.addBox(xa - 1, -2 * r - 1, -2 * r - 1, length + 2, 2 * r + 1, 4 * r + 2)
        out, _ = occ.cut([(3, c)], [(3, cut_box)])
        return out

    def half_sphere(xc, r, zc=0.0):
        sph = occ.addSphere(xc, 0, zc, r)
        cut_box = occ.addBox(xc - 2 * r, -2 * r, zc - 2 * r, 4 * r, 2 * r, 4 * r)
        if dom.full:
            # Two hemispheres: a whole sphere centred on the tube axis failed to imprint at the
            # downstream tip (no tip-sizing face there), which would coarsen the mouth mesh
            other = occ.copy([(3, sph)])
            top, _ = occ.cut([(3, sph)], [(3, cut_box)])
            cut_top = occ.addBox(xc - 2 * r, 0, zc - 2 * r, 4 * r, 2 * r, 4 * r)
            bottom, _ = occ.cut(other, [(3, cut_top)])
            return top + bottom
        out, _ = occ.cut([(3, sph)], [(3, cut_box)])
        return out

    def half_disc(xc, r):
        d = occ.addDisk(xc, 0, 0, r, r)
        occ.rotate([(2, d)], xc, 0, 0, 0, 1, 0, math.pi / 2)    # normal along x
        if dom.full:
            return [(2, d)]
        cut_box = occ.addBox(xc - 1, -2 * r, -2 * r, 2, 2 * r, 4 * r)
        out, _ = occ.cut([(2, d)], [(3, cut_box)])
        return out

    # No sizing-region or reporting face may lie on a tube end face. A disc coplanar with an
    # end face meets the 90-degree lip, the topology of the rear-edge fault: harmless at s1.00,
    # but at s0.64 the lip prisms collapsed (3.4e-14 m3) and held 2-3 m/s at the inlet tip's
    # inner lip and the inlet tube's inner end. So the jet region ends inside each tube's
    # protrusion, where its discs cross the tube walls at right angles (or, for a protrusion
    # under 4 mm, 5 mm out in the chamber), and the mouth discs sit one radius inside the bores.
    def jet_end(x_tube_end, protrusion, into_chamber):
        if protrusion >= 0.004:
            return x_tube_end - into_chamber * min(0.010, 0.5 * protrusion)
        return x_tube_end + into_chamber * 0.005
    x_jet0 = jet_end(xs["in_end"], xs["in_end"] - xs["cav0"], +1)
    x_jet1 = jet_end(xs["out_inner"], xs["cav1"] - xs["out_inner"], -1)
    xs["mouth_up"], xs["mouth_down"] = xs["in_tip"] + rt_i, xs["out_tip"] - ro_i
    tools = []
    tools += half_cyl(x_jet0, x_jet1 - x_jet0, R_JET)                                 # jet core
    tools += half_sphere(xs["in_tip"], R_TIP)                                         # inlet tip
    tools += half_sphere(xs["out_tip"], R_TIP, g["z_tip"] if (s.elbow_c or s.hole_d) else 0.0)        # outlet tip
    near_end = xs["out_tip"] + 0.15
    r_near = R_NEAR
    if r_near >= rb_o + 0.010:
        tools += half_cyl(xs["front"] - 0.15, near_end - xs["front"] + 0.15, r_near)  # near body
    if s.glass_len:
        tools += half_sphere(xs["glass_mouth"], R_TIP)                                # glass mouth
    # The wake region must stay clear of the domain walls: in the 0.15 m half-width flume a
    # 0.16 m radius cut the cross-section corners in thin wedges (~20 deg) and the steady
    # solution converged ~10x worse than in the open domain
    r_wake = min(R_WAKE, 0.8 * min(H, -z_lo, z_hi))
    # The wake region starts where the near region ends, not at the body's rear face. Started at
    # the rear face, its first disc was coplanar with that face and met the side wall at the
    # 90-degree rear edge, so the prism layers from the side wall and the rear face could not
    # fan round the edge. They collapsed into ~190 degenerate cells (faces down to 1e-14 m2)
    # holding 300-770 m/s in every mesh, and the flume runs diverged from them. Starting
    # at near_end also avoids the thin annulus a 0.12 m flume wake region left against the
    # 0.11 m near region.
    x_wake = near_end
    tools += half_cyl(x_wake, min(1.5, x1 - x_wake - 0.2), r_wake)                    # wake
    discs = []
    discs += half_disc(xs["mouth_up"], rt_i)                                          # interior_mouth_upstream
    if s.elbow_c:                                     # mouth one radius up the L-tube (z-normal)
        xs["mouth_down_z"] = g["z_tip"] + rt_i
        dz = occ.addDisk(g["xe"], 0, xs["mouth_down_z"], rt_i, rt_i)
        discs += ([(2, dz)] if dom.full else occ.cut([(2, dz)], [(3, occ.addBox(
            g["xe"] - 1, -2 * rt_i, xs["mouth_down_z"] - 1, 2, 2 * rt_i, 2))])[0])
    elif s.hole_d:                                    # mouth one hole radius inside the hole (x-normal)
        hr = s.hole_d / 2
        xs["mouth_down"] = xs["rear"] - hr
        dh = occ.addDisk(xs["mouth_down"], 0, s.hole_z, hr, hr)
        occ.rotate([(2, dh)], xs["mouth_down"], 0, s.hole_z, 0, 1, 0, math.pi / 2)
        discs += ([(2, dh)] if dom.full else occ.cut([(2, dh)], [(3, occ.addBox(
            xs["mouth_down"] - 1, -2 * hr, s.hole_z - 1, 2, 2 * hr, 2))])[0])
    else:
        discs += half_disc(xs["mouth_down"], ro_i)                                    # interior_mouth_downstream
    discs += half_disc(0.5 * (xs["cav0"] + xs["cav1"]), rb_i)                        # interior_body_mid
    if s.glass_len:                                   # one glass radius inside the glass mouth
        xs["mouth_glass"] = xs["glass_mouth"] + rg_i
        discs += half_disc(xs["mouth_glass"], rg_i)                                   # interior_mouth_glass
    occ.synchronize()

    # Fragment the fluid by the tools, then discard tool pieces lying inside the solid sampler
    _, out_map = occ.fragment(fluid, tools + discs)
    keep = {t for i in range(len(fluid)) for d, t in out_map[i] if d == 3}
    drop = [(3, t) for d, t in occ.getEntities(3) if t not in keep]
    occ.remove(drop, recursive=True)
    occ.synchronize()

    vols = [t for _, t in gmsh.model.getEntities(3)]
    surfs = [t for _, t in gmsh.model.getEntities(2)]

    # --- classify surfaces into zones ----------------------------------------------------
    tol = 1e-6
    zones = {}

    def add(zone, tag):
        zones.setdefault(zone, []).append(tag)

    x_mid = 0.5 * (xs["cav0"] + xs["cav1"])
    skipped = []                                         # bounding boxes of discarded faces
    for sf in surfs:
        up, _ = gmsh.model.getAdjacencies(2, sf)
        nvol = len([v for v in up if v in vols])
        bx0, by0, bz0, bx1, by1, bz1 = gmsh.model.getBoundingBox(2, sf)
        rmax = max(abs(by1), abs(by0), abs(bz0), abs(bz1))
        if nvol == 0:
            skipped.append([round(v, 5) for v in (bx0, by0, bz0, bx1, by1, bz1)])
            continue                                     # stray (e.g. trimmed tool face)
        if nvol == 2:                                    # interior face
            if s.elbow_c and abs(bz1 - bz0) < tol and abs(bz0 - xs["mouth_down_z"]) < 1e-5 \
                    and max(abs(bx0 - g["xe"]), abs(bx1 - g["xe"])) <= rt_i + 1e-5:
                add("interior_mouth_downstream", sf)
            elif abs(bx1 - bx0) < tol and rmax <= rt_i + 1e-5 and abs(bx0 - xs["mouth_up"]) < 1e-5:
                add("interior_mouth_upstream", sf)
            elif s.hole_d and abs(bx1 - bx0) < tol and abs(bx0 - xs["mouth_down"]) < 1e-5 \
                    and abs(bz0 - s.hole_z) <= s.hole_d / 2 + 1e-5 and abs(bz1 - s.hole_z) <= s.hole_d / 2 + 1e-5:
                add("interior_mouth_downstream", sf)
            elif s.hole_d and bx0 >= xs["cav1"] - 1e-5 and bx1 <= xs["rear"] + 1e-5 \
                    and abs(bz0 - s.hole_z) <= s.hole_d / 2 + 1e-5 and abs(bz1 - s.hole_z) <= s.hole_d / 2 + 1e-5 \
                    and max(abs(by0), abs(by1)) <= s.hole_d / 2 + 1e-5:
                add("interior_bore", sf)
            elif abs(bx1 - bx0) < tol and rmax <= rt_i + 1e-5 and abs(bx0 - xs["mouth_down"]) < 1e-5:
                add("interior_mouth_downstream", sf)
            elif abs(bx1 - bx0) < tol and abs(bx0 - x_mid) < 1e-5 and rmax <= rb_i + 1e-5:
                add("interior_body_mid", sf)
            elif s.glass_len and abs(bx1 - bx0) < tol and abs(bx0 - xs["mouth_glass"]) < 1e-5 \
                    and rmax <= rg_i + 1e-5:
                add("interior_mouth_glass", sf)
            elif s.glass_len and rmax <= rg_i + 1e-5 and bx0 >= xs["glass_mouth"] - 1e-5 \
                    and bx1 <= xs["in_tip"] + 1e-5:
                add("interior_bore", sf)          # glass-bore cross-sections cut by sizing regions
            elif rmax <= max(rt_i, s.pass_d / 2) + 1e-5 or (s.elbow_c and max(abs(bx0 - g["xe"]), abs(bx1 - g["xe"])) <= s.pass_d / 2 + 1e-5
                                                       and bz1 <= s.pass_d / 2 + 1e-5) \
                    or (s.dev_len and rmax <= s.dev_bore / 2 + 1e-5 and bx0 >= xs["in_inner"] - 1e-5
                        and bx1 <= xs["in_end"] + 1e-5):
                add("interior_bore", sf)          # bore cross-sections cut by sizing regions
            elif (s.elbow_c or s.hole_d) and bx0 >= g["xe"] - R_TIP - 1e-5 and bx1 <= g["xe"] + R_TIP + 1e-5 \
                    and bz0 >= g["z_tip"] - R_TIP - 1e-5 and bz1 <= g["z_tip"] + R_TIP + 1e-5 and rmax > R_JET + 1e-5:
                add("interior_tip", sf)                   # off-axis outlet tip sphere
            elif rmax <= R_JET + 1e-5:
                add("interior_jet", sf)
            elif rmax <= R_TIP + 1e-5:
                add("interior_tip", sf)
            elif rmax <= R_NEAR + 1e-5:
                add("interior_near", sf)
            else:
                add("interior_wake", sf)
            continue
        if abs(by1) < tol and abs(by0) < tol:
            add("symmetry", sf)
        elif abs(bx1 - x0) < tol:
            add("inlet", sf)
        elif abs(bx0 - x1) < tol:
            add("outlet", sf)
        # Flume side wall and floor are free-slip (symmetry) in this first blockage test: prism
        # layers on them crashed Fluent Meshing 2025 R1 (SIGSEGV), and blockage is essentially
        # inviscid. Wall boundary layers (displacement ~4 mm after 1 m) are a stated omission.
        elif dom.confined and abs(by0 - H) < tol:
            add("symmetry_flumewall_side" if dom.flume_walls else "symmetry_flume_side", sf)
        elif dom.confined and abs(bz1 - z_lo) < tol:
            add("symmetry_flumewall_floor" if dom.flume_walls else "symmetry_flume_floor", sf)
        elif (abs(by0 - H) < tol or abs(by1 + H) < tol or abs(bz0 - z_lo) < tol
              or abs(bz1 - z_hi) < tol):
            add("symmetry_far", sf)
        elif s.elbow_c:
            add(elliott_wall_zone(s, xs, g, (bx0, by0, bz0, bx1, by1, bz1)), sf)
        elif s.hole_d:
            add(mcdonald_wall_zone(s, xs, g, (bx0, by0, bz0, bx1, by1, bz1)), sf)
        else:                                            # sampler walls
            in_tube_x = bx0 >= xs["in_tip"] - 1e-5 and bx1 <= xs["in_end"] + 1e-5
            out_tube_x = bx0 >= xs["out_inner"] - 1e-5 and bx1 <= xs["out_tip"] + 1e-5
            in_glass_x = s.glass_len and bx0 >= xs["glass_mouth"] - 1e-5 and bx1 <= xs["glass_end"] + 1e-5
            flat_x = abs(bx1 - bx0) < tol
            if in_glass_x and flat_x and abs(bx0 - xs["in_tip"]) < 1e-5 and rmax <= rg_i + 1e-5:
                add("wall_glass_step", sf)                    # glass bore onto the tube tip
            elif in_glass_x and flat_x and abs(bx0 - xs["glass_mouth"]) < 1e-5:
                add("wall_glass_mouth", sf)                   # end annulus facing the flow
            elif in_glass_x and rmax <= rg_i + 1e-5 and bx1 <= xs["in_tip"] + 1e-5:
                add("wall_glass_bore", sf)
            elif in_glass_x and rmax <= rg_o + 1e-5 and rmax > rg_i + 1e-5:
                add("wall_glass_outer", sf)
            elif s.dev_len and rmax <= s.dev_bore / 2 + 1e-5 and bx0 >= xs["in_inner"] - 1e-5 \
                    and bx1 <= xs["in_end"] + 1e-5:
                add("wall_upstream_bore", sf)                 # device bore and its steps
            elif rmax <= rt_i + 1e-5 and in_tube_x:
                add("wall_upstream_bore", sf)
            elif rmax <= ro_i + 1e-5 and out_tube_x:
                add("wall_downstream_bore", sf)
            elif s.dev_len and bx0 >= xs["in_inner"] + s.dev_gap - 1e-5 and bx1 <= xs["in_end"] + 1e-5:
                add("wall_fitting", sf)                       # device body and ports (outside)
            elif rmax <= rt_o + 1e-5 and in_tube_x:
                add("wall_upstream_tube", sf)
            elif rmax <= ro_o + 1e-5 and out_tube_x:
                add("wall_downstream_tube", sf)
            elif bx1 <= 1e-5:
                add("wall_cone", sf)
            elif rmax <= rb_i + 1e-5 and bx0 >= xs["cav0"] - 1e-5 and bx1 <= xs["cav1"] + 1e-5:
                add("wall_cavity", sf)
            elif bx0 >= xs["rear"] - 1e-5 and s.tail_cone_len > 0:
                add("wall_tail_cone", sf)
            elif abs(bx0 - xs["rear"]) < 1e-5 and abs(bx1 - xs["rear"]) < 1e-5:
                add("wall_rear", sf)
            else:
                add("wall_body", sf)

    for name, tags in zones.items():
        gmsh.model.addPhysicalGroup(2, tags, name=name)

    # --- mesh size fields -----------------------------------------------------------------
    k = size_scale
    f = gmsh.model.mesh.field
    tube_faces = sum((zones.get(z, []) for z in
                      ("wall_upstream_bore", "wall_downstream_bore", "wall_upstream_tube",
                       "wall_downstream_tube", "interior_mouth_upstream", "interior_mouth_downstream",
                       "outlet_tube_up", "inlet_tube_down", "inlet_bore", "outlet_bore",
                       "wall_glass_bore", "wall_glass_step", "interior_mouth_glass")), [])
    # Sizes grow linearly with distance from each surface family (size = s0 + g * d)
    d_tube = f.add("Distance"); f.setNumbers(d_tube, "SurfacesList", tube_faces)
    f.setNumber(d_tube, "Sampling", 40)
    t_tube = f.add("MathEval")
    f.setString(t_tube, "F", f"{0.0003 * k} + 0.12*F{d_tube}")
    d_body = f.add("Distance")
    f.setNumbers(d_body, "SurfacesList", sum((zones.get(z, []) for z in
                 ("wall_cone", "wall_body", "wall_rear", "wall_cavity", "wall_tail_cone", "wall_fitting")), []))
    t_body = f.add("MathEval")
    f.setString(t_body, "F", f"{0.003 * k} + 0.15*F{d_body}")
    jet = f.add("Box")
    for key, val in dict(VIn=0.0012 * k, VOut=1.0, XMin=xs["in_end"], XMax=xs["out_inner"],
                         YMin=-0.0121 if dom.full else 0, YMax=0.0121, ZMin=-0.0121, ZMax=0.0121).items():
        f.setNumber(jet, key, val)
    fields = [t_tube, t_body, jet]
    glass_out = sum((zones.get(z, []) for z in ("wall_glass_outer", "wall_glass_mouth")), [])
    if glass_out:                        # the 600 mm outer glass surface: 1 mm faces, not bore-sized
        d_glass = f.add("Distance"); f.setNumbers(d_glass, "SurfacesList", glass_out)
        t_glass = f.add("MathEval")
        f.setString(t_glass, "F", f"{0.001 * k} + 0.12*F{d_glass}")
        fields.append(t_glass)
    fmin = f.add("Min"); f.setNumbers(fmin, "FieldsList", fields)
    f.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.04 * k)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / case_name
    gmsh.write(str(stem.with_suffix(".step")))
    report = write_outputs(stem, zones, write_stl)
    report["spec"] = asdict(s)
    report["domain"] = asdict(dom)
    report["stations_m"] = xs
    report["size_scale"] = size_scale
    report["skipped_faces"] = skipped
    report["volumes"] = [[round(v, 4) for v in gmsh.model.getBoundingBox(3, t)] for t in vols]
    stem.with_name(case_name + "_zones.json").write_text(json.dumps(report, indent=2))
    if gui:
        gmsh.fltk.run()
    gmsh.finalize()
    return report


FLUENT_BC = {"inlet": ("a", "velocity-inlet"), "outlet": ("5", "pressure-outlet"),
             "symmetry": ("7", "symmetry"), "symmetry_far": ("7", "symmetry"),
             "symmetry_flume_side": ("7", "symmetry"), "symmetry_flume_floor": ("7", "symmetry"),
             "symmetry_flumewall_side": ("7", "symmetry"), "symmetry_flumewall_floor": ("7", "symmetry")}


def write_outputs(stem: Path, zones: dict, write_stl=False):
    """Write the Fluent ASCII boundary mesh and a multi-solid STL from the gmsh triangulation."""
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)
    index = {int(t): i + 1 for i, t in enumerate(node_tags)}     # Fluent nodes are 1-based
    zone_tris, report = {}, {"zones": {}}
    for name, tags in zones.items():
        tris = []
        for sf in tags:
            types, _, enodes = gmsh.model.mesh.getElements(2, sf)
            for et, en in zip(types, enodes):
                if et != 2:
                    raise RuntimeError(f"non-triangle surface element in {name}")
                tris.append(np.array(en, dtype=np.int64).reshape(-1, 3))
        tri = np.vstack(tris)
        zone_tris[name] = tri
        p = xyz[np.vectorize(index.get)(tri) - 1]
        area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum()
        report["zones"][name] = {"faces": int(len(tri)), "area_m2": float(area)}

    n_faces = sum(len(t) for t in zone_tris.values())
    lines = ['(0 "TIMS boundary mesh written by tims_geometry.py")', "(2 3)",
             f"(10 (0 1 {len(xyz):x} 0 3))", f"(10 (1 1 {len(xyz):x} 1 3)("]
    lines += [f"{x:.9e} {y:.9e} {z:.9e}" for x, y, z in xyz]
    lines += ["))", f"(13 (0 1 {n_faces:x} 0))"]
    start, zid, names = 1, 3, []
    for name, tri in zone_tris.items():
        internal = (name in ("interior_mouth_upstream", "interior_mouth_downstream", "interior_body_mid")
                    or name.startswith("interior"))
        # Fluent Meshing skips "interior" faces when reading a boundary mesh, so internal
        # surfaces are written as walls; the watertight workflow retypes faces shared by two
        # fluid regions as interior when it builds the regions.
        bc_hex, bc_name = FLUENT_BC.get(name, ("3", "wall"))
        end = start + len(tri) - 1
        lines.append(f"(13 ({zid:x} {start:x} {end:x} {bc_hex} 3)(")
        lines += [f"{index[a]:x} {index[b]:x} {index[c]:x} 0 0" for a, b, c in tri]
        lines.append("))")
        names.append(f"(45 ({zid} {bc_name} {name})())")
        report["zones"][name]["fluent_id"] = zid
        report["zones"][name]["fluent_type"] = "interior" if internal else bc_name
        start, zid = end + 1, zid + 1
    stem.with_suffix(".msh").write_text("\n".join(lines + names) + "\n")

    for name, tri in (zone_tris.items() if write_stl else []):
        with open(stem.with_name(f"{stem.name}_{name}.stl"), "w") as fh:
            fh.write(f"solid {name}\n")
            for a, b, c in tri:
                p0, p1, p2 = xyz[index[a] - 1], xyz[index[b] - 1], xyz[index[c] - 1]
                nrm = np.cross(p1 - p0, p2 - p0)
                nrm = nrm / (np.linalg.norm(nrm) or 1.0)
                fh.write(f" facet normal {nrm[0]:e} {nrm[1]:e} {nrm[2]:e}\n  outer loop\n")
                for q in (p0, p1, p2):
                    fh.write(f"   vertex {q[0]:.9e} {q[1]:.9e} {q[2]:.9e}\n")
                fh.write("  endloop\n endfacet\n")
            fh.write(f"endsolid {name}\n")
    report["nodes"], report["faces"] = int(len(xyz)), int(n_faces)
    return report


def build_duct(dom: DomainSpec, out_dir: Path, size_scale: float = 1.0, length: float = 2.0):
    """Empty flume section with the domain's cross-section, for the inflow precursor.

    The developed flow of the bare flume is computed in this duct by recycling its outlet profile
    to its inlet until it stops changing (run_tims_case_v3.py precursor); the result is mapped onto
    the sampler domain's inlet. Its boundaries are named as in the sampler domain, so the flume
    walls get the same face size (and so the same wall-function resolution) in both.
    """
    if not dom.confined:
        raise ValueError("a duct precursor is built for a confined (flume) domain only")
    H = dom.half_width
    z_lo, z_hi = (-dom.floor, dom.lid) if dom.floor else (-H, H)
    name = f"duct_{dom.name}"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    occ = gmsh.model.occ
    occ.addBox(0.0, 0.0, z_lo, length, H, z_hi - z_lo)
    occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.02 * size_scale)
    gmsh.model.mesh.generate(2)
    zones, tol = {}, 1e-6
    for _, sf in gmsh.model.getEntities(2):
        bx0, by0, bz0, bx1, by1, bz1 = gmsh.model.getBoundingBox(2, sf)
        if abs(bx1) < tol:
            zone = "inlet"
        elif abs(bx0 - length) < tol:
            zone = "outlet"
        elif abs(by1) < tol:
            zone = "symmetry"
        elif abs(by0 - H) < tol:
            zone = "symmetry_flumewall_side" if dom.flume_walls else "symmetry_flume_side"
        elif abs(bz1 - z_lo) < tol:
            zone = "symmetry_flumewall_floor" if dom.flume_walls else "symmetry_flume_floor"
        else:
            zone = "symmetry_far"
        zones.setdefault(zone, []).append(sf)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / name
    report = write_outputs(stem, zones)
    report.update(domain=asdict(dom), size_scale=size_scale, duct=True, length=length,
                  skipped_faces=[])
    stem.with_name(name + "_zones.json").write_text(json.dumps(report, indent=2))
    gmsh.finalize()
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=sorted(SPECS))
    ap.add_argument("--size-scale", type=float, default=1.0)
    ap.add_argument("--out", default=str(Path(__file__).parent / "out"))
    ap.add_argument("--gui", action="store_true")
    ap.add_argument("--stl", action="store_true")
    ap.add_argument("--domain", choices=sorted(DOMAINS), default="open")
    ap.add_argument("--duct", action="store_true",
                    help="build only the empty flume section of --domain (inflow precursor)")
    a = ap.parse_args()
    if a.duct:
        rep = build_duct(DOMAINS[a.domain], Path(a.out), a.size_scale)
    else:
        rep = build(SPECS[a.case], DOMAINS[a.domain], Path(a.out), a.size_scale, a.gui, a.stl)
    print(json.dumps({k: v for k, v in rep["zones"].items()}, indent=1))
    print("nodes", rep["nodes"], "faces", rep["faces"])
