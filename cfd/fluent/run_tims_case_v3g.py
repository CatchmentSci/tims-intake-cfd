"""Mesh and solve a TIMS hydrodynamics case in Ansys Fluent through PyFluent.

v3g: run_tims_case_v3.py plus the tube-area fix of the decoupled stages,
the repeatable --laminar-bore with the split after --init-data, and an overwrite guard (--overwrite).
--laminar-zone takes whole the zones lying wholly inside a laminar-bore cylinder (the
re-screens left them SST-gamma), and --bore-volume stops the run at setup unless each bore's
laminar zones hold 99-101 % of its cylinder's cell volume (geometry/laminar_zones.py).
The split after --init-data appended its new interface
surfaces to every surface report already defined (mouth pressures 228/158 Pa instead of 433/40 at
the same instant). After the split every report's surface list and the force-report zones are
reset to their setup values, and the run stops unless every report reads the same after the split
as before it (relative 1e-9; the field does not change in between).

Stages (run in order on the HPC; each stage writes its own outputs and log):

  mesh    Fluent Meshing watertight workflow from the named boundary mesh made by
          tims_geometry.py: local face sizing, surface remesh, prism layers, poly-hexcore.
  setup   Reads the volume mesh, applies models and boundary conditions, runs 20 iterations at
          one velocity and writes the case. Use it to confirm the settings paths on your
          Fluent version before committing to the sweep.
  sweep   Steady solutions over a list of ambient velocities. Each velocity starts from the
          previous solution and runs until the inlet-tube flow rate is stationary.
  transient  URANS run at one velocity: steady pre-iterations, a short-step ramp, then time
             stepping with R averaged over the second half of the run (the production runs).
  precursor  Developed flume inflow for the Goharrokhi flume (recycled duct profile).

Outputs per velocity (in --out): case/data files and one row in results.csv holding
R = Q_tube / (A_tube U), the pressure coefficients at the tube mouths, the sampler drag
coefficient, the mass imbalance between the two mouths and the convergence history.

Settings paths in the Fluent settings API have been renamed between releases. Every write
goes through set_first(), which tries the known spellings in order and logs the one accepted,
so a renamed path fails with a clear message instead of silently doing nothing.

Example:
  python run_tims_case_v3.py mesh --geometry ../geometry/out/v5/scale_1.00/phillips2000.msh --cores 32
  python run_tims_case_v3.py transient --mesh out/v5/phillips2000_s1.00.msh.h5 --velocities 0.3
      --zones ../geometry/out/v5/scale_1.00/phillips2000_zones.json --cores 32   (one command)
"""
import argparse
import csv
import json
import logging
import math
import os
import random
import re
import sys
import time
from pathlib import Path

import ansys.fluent.core as pyfluent
import grpc

# gRPC keepalive: runs died with "keepalive watchdog timeout" while the solution was
# healthy, several on different nodes in the same minute: the Fluent host stopped
# answering for longer than gRPC's default 20 s ping-ack timeout, most likely in a shared-filesystem
# stall. PyFluent sets no keepalive options, so its channels are created here with a 10 min ack
# timeout to ride out such stalls. Only the timeout is changed; the ping interval stays at gRPC's
# default (a shorter interval could trip the server's too-many-pings limit).
# This alone did not help: the watchdog that fires is the Fluent server's (it pings the client
# every 10 s), which no client option changes; see dual_time() for the fix. Kept, harmless.
KEEPALIVE_OPTS = [("grpc.keepalive_timeout_ms", 600000)]
_insecure_channel = grpc.insecure_channel


def _channel_with_keepalive(target, options=None, compression=None):
    opts = list(options or [])
    names = {k for k, _ in opts}
    opts += [o for o in KEEPALIVE_OPTS if o[0] not in names]
    return _insecure_channel(target, options=opts, compression=compression)


grpc.insecure_channel = _channel_with_keepalive

# Loopback server address: PyFluent infers REMOTING_SERVER_ADDRESS by a 2 s gRPC health
# check on "localhost"; on a loaded node that check can time out, Fluent then serves on the node's
# network address, and PyFluent refuses it ("Connecting to remote Fluent instances is not allowed").
# The client always runs on the Fluent host node,
# so the address is pinned to loopback, which also keeps the session off the node's network.
if not pyfluent.config.launch_fluent_ip:
    pyfluent.config.launch_fluent_ip = "127.0.0.1"

LOG = logging.getLogger("tims")

# Water at 20 C (Kestin et al. 1978; IAPWS): Phillips et al. (2000) calibrated in a laboratory flume
RHO, MU = 998.2, 1.002e-3

WALLS = ["wall_upstream_bore", "wall_upstream_tube", "wall_downstream_bore",
         "wall_downstream_tube", "wall_cone", "wall_body", "wall_rear", "wall_cavity",
         "wall_tail_cone"]
EXTERNAL_WALLS = ["wall_upstream_tube", "wall_downstream_tube", "wall_cone", "wall_body",
                  "wall_rear", "wall_tail_cone"]

# Sweep velocities (m/s): Phillips et al. (2000) calibration range 0.154-0.585 m/s, Smith & Owens
# (2014) flume 0.28 m/s, field maxima ~1 m/s (Phillips) and >2.5 m/s (McDonald et al. 2010)
DEFAULT_U = [0.10, 0.154, 0.20, 0.28, 0.30, 0.40, 0.50, 0.585, 0.75, 1.00, 1.50, 2.00, 2.50]


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------
def ride(fn, what, max_wait=3600):
    """Call fn(); if the connection to Fluent drops, repeat it every 30 s until Fluent answers.
    Only for idempotent settings reads and writes (a drop during the time-step assignment after
    the steady phase otherwise ends the run with "failed to connect to all addresses"). Matched
    on the message whatever the exception type, since the settings API wraps gRPC errors."""
    t0 = None
    while True:
        try:
            out = fn()
            if t0 is not None:
                LOG.warning("%s answered again after %.0f s", what, time.time() - t0)
            return out
        except Exception as exc:                       # noqa: BLE001
            if not any(k in str(exc) for k in SESSION_DROP):
                raise
            if t0 is None:
                t0 = time.time()
                LOG.warning("Fluent dropped the connection during %s (%s); retrying", what, exc)
            elif time.time() - t0 > max_wait:
                raise RuntimeError(f"{what}: no answer from Fluent in {max_wait} s") from exc
            time.sleep(30)


def resolve(root, path):
    """Walk a dotted settings path such as "setup.materials.fluid['water-liquid'].density"."""
    def walk():
        obj = root
        for attr, key in re.findall(r"(\w+)(?:\['([^']+)'\])?", path):
            obj = getattr(obj, attr)
            if key:
                obj = obj[key]
        return obj
    return ride(walk, f"settings lookup {path}")


def set_first(root, paths, value, what):
    """Assign value at the first settings path that exists; raise listing all attempts."""
    errors = []
    for p in paths:
        try:
            parent, _, leaf = p.rpartition(".")
            target = resolve(root, parent) if parent else root
            ride(lambda: setattr(target, leaf, value), f"setting {what}")
            LOG.info("set %-40s %s = %r", what, p, value)
            return p
        except Exception as exc:          # noqa: BLE001 - try the next spelling
            errors.append(f"{p}: {type(exc).__name__}: {exc}")
    raise RuntimeError(f"could not set {what}; tried:\n  " + "\n  ".join(errors))


def find_child(root, path, contains):
    """Dotted path of the first setting below path whose name contains the given text."""
    stack = [(path, resolve(root, path))]
    while stack:
        p, obj = stack.pop(0)
        for name in getattr(obj, "child_names", []):
            if contains in name:
                return f"{p}.{name}"
            try:
                stack.append((f"{p}.{name}", getattr(obj, name)))
            except Exception:                          # noqa: BLE001 - inactive branch
                pass
    raise RuntimeError(f"no setting containing {contains!r} below {path}")


def choose_below(root, path, wanted, what):
    """choose() on the first active option below path whose allowed values match wanted."""
    stack = [path]
    while stack:
        p = stack.pop(0)
        try:
            obj = resolve(root, p)
            if hasattr(obj, "is_active") and not ride(obj.is_active, f"{what} lookup"):
                continue
            allowed = ride(obj.allowed_values, f"{what} lookup") if hasattr(obj, "allowed_values") else []
        except Exception:                              # noqa: BLE001 - inactive branch
            continue
        if any(w.lower() in v.lower() for w in wanted for v in allowed):
            return choose(root, p, wanted, what)
        stack += [f"{p}.{n}" for n in getattr(obj, "child_names", [])]
    raise RuntimeError(f"{what}: no option offering {wanted} below {path}")


def choose(root, path, wanted, what):
    """Set a string option to the first allowed value containing one of the wanted substrings."""
    obj = resolve(root, path)
    allowed = ride(obj.allowed_values, f"{what} lookup")
    for w in wanted:
        for v in allowed:
            if w.lower() in v.lower():
                ride(lambda: obj.set_state(v), f"setting {what}")
                LOG.info("set %-40s %s = %r (allowed: %s)", what, path, v, allowed)
                return v
    raise RuntimeError(f"{what}: none of {wanted} in allowed values {allowed}")


def launch(mode, cores, precision="double"):
    extra = os.environ.get("TIMS_FLUENT_ARGS", "")      # e.g. "-cnf=hosts.txt -mpi=openmpi"
    # start_timeout: the default 60 s timed out on a busy node;
    # dropped by the fallback spellings if this PyFluent does not know it.
    attempts = [dict(mode=mode, precision=precision, processor_count=cores, ui_mode="no_gui",
                     additional_arguments=extra, cleanup_on_exit=True, start_timeout=600),
                dict(mode=mode, precision=precision, processor_count=cores, ui_mode="no_gui",
                     additional_arguments=extra, cleanup_on_exit=True),
                dict(mode=mode, precision=precision, processor_count=cores, show_gui=False,
                     additional_arguments=extra, cleanup_on_exit=True)]
    last = None
    for tries in range(3):                             # retry a failed launch, staggered
        for kw in attempts:
            try:
                session = pyfluent.launch_fluent(**kw)
                LOG.info("launched Fluent %s (%s cores) with %s", mode, cores, sorted(kw))
                return session
            except TypeError as exc:                   # keyword not known to this PyFluent
                last = exc
            except Exception as exc:                   # noqa: BLE001  launch refused or timed out
                last = exc
                wait = 30 + random.uniform(0, 60)
                LOG.warning("Fluent launch failed (%s: %s); retry in %.0f s", type(exc).__name__,
                            exc, wait)
                time.sleep(wait)
                break
    raise last


def guard_out(out, stem, args):
    """Overwrite guard (v3g): refuse to start a run whose
    output folder already holds a final, checkpoint, steady state or summary for this stem,
    unless --overwrite is given. Prints the job that wrote the existing files; makes no copy."""
    clash = [f"{stem}_{s}" for s in ("final.cas.h5", "final.dat.h5", "ckpt.cas.h5", "steady.cas.h5", "summary.json")
             if (out / f"{stem}_{s}").exists()]
    if not clash or getattr(args, "overwrite", False):
        if clash:
            LOG.warning("--overwrite: replacing %s in %s", clash, out)
        return
    prev = "unknown (metadata has no job id)"
    meta = out / f"{stem}_run_metadata.json"
    if meta.exists():
        try:
            m = json.loads(meta.read_text())
            prev = f"{m.get('slurm_job') or 'unrecorded'} ({m.get('time', 'no time')})"
        except Exception:                              # noqa: BLE001
            pass
    msg = (f"refusing to overwrite {clash} in {out}: written by job {prev}. "
           f"Use a new TAG, or --overwrite if this is intended.")
    LOG.error(msg)
    raise SystemExit(msg)


def metadata(session, args, extra=None):
    meta = {"args": vars(args), "pyfluent": pyfluent.__version__, "time": time.ctime(),
            "rho": RHO, "mu": MU, "slurm_job": os.environ.get("SLURM_JOB_ID", "")}
    try:
        meta["fluent"] = str(session.get_fluent_version())
    except Exception:                                  # noqa: BLE001
        meta["fluent"] = "unknown"
    meta.update(extra or {})
    return meta


# ---------------------------------------------------------------------------------------------
# stage 1: meshing
# ---------------------------------------------------------------------------------------------
# Face sizes (m) by geometry zone name, scaled with the geometry size factor for the GCI study.
# Bore: >= 16 cells across the 4 mm diameter; jet core: ~10 cells across a 12 mm radius.
FACE_SIZES = [
    ("bore", ("wall_upstream_bore", "wall_downstream_bore", "interior_mouth_upstream",
              "interior_mouth_downstream", "interior_bore",
              # calibration glass tube (Phillips rig): 28 faces across the 7 mm bore
              "wall_glass_bore", "wall_glass_step", "interior_mouth_glass"), 0.00025),
    ("glass", ("wall_glass_outer", "wall_glass_mouth"), 0.001),
    ("tube", ("wall_upstream_tube", "wall_downstream_tube"), 0.0004),
    ("tip", ("interior_tip",), 0.0015),
    ("jet", ("interior_jet",), 0.0012),
    ("cone", ("wall_cone", "wall_tail_cone"), 0.002),
    ("body", ("wall_body", "wall_rear", "wall_cavity", "interior_body_mid"), 0.003),
    ("near", ("interior_near",), 0.006),
    ("wake", ("interior_wake",), 0.010),
    ("fitting", ("wall_fitting",), 0.0015),          # Elliott fittings (22 mm nut, 21 mm elbow)
    # Goharrokhi flume bed and side wall: no prism layers (they crashed Fluent Meshing on
    # flume walls); 8 mm faces put the first cell centre ~3 mm off the wall, y+ ~20 at
    # 0.15 m/s to ~75 at 0.6 m/s, inside the range of SST's y+-insensitive wall treatment
    ("flumewall", ("symmetry_flumewall_side", "symmetry_flumewall_floor"), 0.008),
]
# Flume walls are meshed under symmetry names (no prisms grown) and retyped as walls in the solver
FLUME_WALLS = ("symmetry_flumewall_side", "symmetry_flumewall_floor")
# Bore-only refinement test: the Phillips GCI trio showed the bore's wall friction, not
# the outer flow, carries the mesh error. TIMS_BORE_SCALE < 1 shrinks the bore face size only.
BORE_SCALE = float(os.environ.get("TIMS_BORE_SCALE", "1"))
FACE_SIZES = [(n, z, s * BORE_SCALE if n == "bore" else s) for n, z, s in FACE_SIZES]
# Reporting planes are written as walls (Fluent Meshing drops "interior" faces from a boundary
# mesh). Describe Geometry retypes zones named *interior*; if a plane still reaches the solver as
# a wall/shadow baffle (which would seal the sampler), setup_solver() converts it to interior.
INTERNAL_ZONES = ("interior_mouth_upstream", "interior_mouth_downstream", "interior_body_mid",
                  "interior_mouth_glass")


def face_zones(meshing):
    """{zone name: zone type} for all face zones. Region splitting renames pieces 'name:id'."""
    mu = meshing.meshing_utilities
    return {str(meshing.scheme_eval.scheme_eval(f"(zone-id->name {z})")): mu.get_zone_type(zone_id=z)
            for z in mu.get_face_zones(filter="*")}


def solid_regions(meshing):
    """Regions bounded by walls alone: the sampler's own solid parts.

    In the half model every solid part is open to the y = 0 plane, which carries no faces, so
    region detection never closes it. In the full domain the solid parts are closed and the
    fluid-only workflow meshed them as fluid: sealed zones behind two-sided walls with no
    pressure reference, which diverged within ~60 iterations. Every fluid region
    touches an interior, inlet, outlet or symmetry face; a solid one touches none.
    """
    mu = meshing.meshing_utilities
    open_ = [z for z in face_zones(meshing)
             if z.split(":")[0].startswith(("interior", "inlet", "outlet", "symmetry"))]
    every = set(mu.get_regions_of_face_zones(face_zone_name_pattern="*"))
    fluid = set(mu.get_regions_of_face_zones(face_zone_name_list=open_))
    return sorted(every - fluid), sorted(fluid)


def build_volume_mesh(meshing, geo, k):
    """Watertight workflow on a named boundary mesh (validated on Fluent 2025 R1, PyFluent 0.37).

    A .msh import arrives as one object whose labels are the auto-detected *regions*, not the
    face zones, so local sizing and boundary types are scoped by zone name instead.
    """
    wf = meshing.workflow
    wf.InitializeWorkflow(WorkflowType="Watertight Geometry")
    imp = wf.TaskObject["Import Geometry"]
    imp.Arguments = dict(FileFormat="Mesh", MeshFileName=str(geo), LengthUnit="m", MeshUnit="m")
    imp.Execute()

    zt = face_zones(meshing)
    ls = wf.TaskObject["Add Local Sizing"]
    for name, bases, size in FACE_SIZES:
        zl = [z for z in zt if z.split(":")[0] in bases]
        if not zl:
            continue
        ls.Arguments = dict(AddChild="yes", BOIControlName=f"fs_{name}", BOIZoneorLabel="zone",
                            BOIFaceZoneList=zl, BOIGrowthRate=1.15, BOISize=size * k,
                            BOIExecution="Face Size")
        ls.AddChildAndUpdate()
        LOG.info("face size %-5s %.5f m on %d zones", name, size * k, len(zl))

    srf = wf.TaskObject["Generate the Surface Mesh"]
    srf.Arguments = dict(CFDSurfaceMeshControls=dict(MinSize=0.00025 * k * min(BORE_SCALE, 1), MaxSize=0.04 * k,
                                                     GrowthRate=1.15, CellsPerGap=3,
                                                     CurvatureNormalAngle=12))
    srf.Execute()

    solid, fluid = solid_regions(meshing)
    LOG.info("regions: %d fluid, solid (not meshed): %s", len(fluid), solid)

    # Names containing inlet/outlet/symmetry/interior are typed automatically here
    dg = wf.TaskObject["Describe Geometry"]
    dg.UpdateChildTasks(SetupTypeChanged=False)
    if solid:
        dg.Arguments = dict(SetupType="The geometry consists of both fluid and solid regions "
                                      "and/or voids", CappingRequired="No")
    else:
        dg.Arguments = dict(SetupType="The geometry consists of only fluid regions with no voids")
    dg.UpdateChildTasks(SetupTypeChanged=True)
    dg.Execute()

    ur = wf.TaskObject["Update Regions"]
    if solid:
        state = dict(ur.Arguments() or {})
        LOG.info("Update Regions defaults: %s", state)
        names = list(state.get("RegionNameList") or fluid + solid)
        types = list(state.get("RegionTypeList") or ["fluid"] * len(names))
        missing = [r for r in solid if r not in names]
        if missing:
            raise RuntimeError(f"solid regions {missing} not among Update Regions {names}")
        ur.Arguments = dict(OldRegionNameList=names, OldRegionTypeList=types, RegionNameList=names,
                            RegionTypeList=["dead" if n in solid else "fluid" for n in names])
    ur.Execute()
    if solid:
        LOG.info("Update Regions: %s", ur.Arguments())

    # Prism layers on all walls: first height 25 um gives y+ ~ 1 on the body at 1 m/s
    # (u_tau ~ 0.04 m/s); 10 layers at 1.2 span ~0.65 mm, i.e. < 1/3 of the bore radius.
    # The flume precursor duct has no walls (its flume walls carry symmetry names): no layers.
    if any(z.split(":")[0].startswith("wall") for z in face_zones(meshing)):
        bl = wf.TaskObject["Add Boundary Layers"]
        bl.Arguments = dict(BLControlName="walls", OffsetMethodType="uniform",
                            FirstHeight=2.5e-5 * k, NumberOfLayers=10, Rate=1.2)
        bl.AddChildAndUpdate()
    else:
        LOG.info("no wall zones: no boundary layers")

    vm = wf.TaskObject["Generate the Volume Mesh"]
    vm.Arguments = dict(VolumeFill="poly-hexcore",
                        VolumeFillControls=dict(HexMaxCellLength=0.04 * k, PeelLayers=1))
    vm.Execute()
    meshing.tui.mesh.check_mesh()


def stage_mesh(args):
    geo = Path(args.geometry).resolve()
    zones = json.loads(geo.with_name(geo.stem + "_zones.json").read_text())
    k = zones["size_scale"]
    out = Path(args.out).resolve() / f"{geo.stem}_s{k:.2f}.msh.h5"
    out.parent.mkdir(parents=True, exist_ok=True)
    meshing = launch("meshing", args.cores)
    build_volume_mesh(meshing, geo, k)
    # write_mesh asks before overwriting, and in batch the old file silently stayed (Elliott
    # remesh), so remove it first and check the new one is there
    t0 = time.time()
    out.unlink(missing_ok=True)
    meshing.tui.file.write_mesh(str(out))
    if not out.exists() or out.stat().st_mtime < t0:
        raise RuntimeError(f"mesh not written: {out}")
    Path(str(out) + ".json").write_text(json.dumps(
        metadata(meshing, args, {"face_sizes": FACE_SIZES, "size_scale": k, "bore_scale": BORE_SCALE}), indent=2))
    LOG.info("wrote %s", out)
    meshing.exit()


# ---------------------------------------------------------------------------------------------
# stage 2/3: solver
# ---------------------------------------------------------------------------------------------
def tube_area(geo_meta, tube_id):
    """(tube cross-section inside the model, mirror factor to the whole sampler)."""
    mirror = 1 if geo_meta.get("domain", {}).get("full") else 2
    return math.pi * tube_id**2 / 4 / mirror, mirror


def setup_solver(solver, args):
    s = solver.settings if hasattr(solver, "settings") else solver
    mesh = str(Path(args.mesh).resolve())
    read_mesh(solver, s, mesh)
    solver.tui.mesh.check()

    # Reporting planes must be open interior faces, never wall baffles (validated on 2025 R1)
    bc = s.setup.boundary_conditions
    try:
        fw = [z for z in bc.symmetry.keys() if z.split(":")[0] in FLUME_WALLS]
    except AttributeError:
        fw = []
    if fw:
        bc.set_zone_type(zone_list=fw, new_type="wall")
        LOG.info("flume walls retyped as no-slip smooth walls: %s", fw)
    # Fluent Meshing keeps a solid region that it cannot delete because fluid encloses it (the
    # sampler body in the full domain); it then comes as a fluid zone, retyped solid here
    czc = s.setup.cell_zone_conditions
    solid = [z for z in czc.fluid.keys() if z.split(":")[0] in args.solid_zones]
    if solid:
        try:
            czc.set_zone_type(zone_list=solid, new_type="solid")
        except (AttributeError, TypeError):
            for z in solid:
                solver.tui.define.boundary_conditions.zone_type(z, "solid")
        LOG.info("retyped cell zones %s as solid", solid)
    # Any other two-sided wall means a solid part was meshed as fluid (see solid_regions)
    shadows = [w for w in bc.wall.keys()
               if w.endswith("-shadow") and w.split(":")[0].split("-shadow")[0] not in INTERNAL_ZONES]
    if shadows and solid:
        LOG.info("walls coupled to the solid zones: %s", shadows)
    elif shadows:
        raise RuntimeError(f"mesh has fluid on both sides of sampler walls {shadows}: solid "
                           "regions were meshed as fluid; remesh")
    baffles = [w for w in bc.wall.keys()
               if w.split(":")[0] in INTERNAL_ZONES and not w.endswith("-shadow")]
    if baffles:
        LOG.warning("converting reporting-plane baffles to interior: %s", baffles)
        bc.set_zone_type(zone_list=baffles, new_type="interior")
    interiors = list(bc.interior.keys())

    def planes(base):
        found = [z for z in interiors if z.split(":")[0] == base]
        if not found:
            raise RuntimeError(f"no interior zone named {base}; interiors are {interiors}")
        return found

    # Laminar-bore control (several cylinders and a converged start):
    # the cells inside each cylinder along a bore are split into their own fluid zones, made
    # laminar with zero turbulent viscosity once the turbulence model is on. With --init-data the
    # split comes after the data are read (the data belong to the unsplit zones), so the run
    # starts from the converged SST-gamma field, not from a hybrid initialisation, which diverged
    # in a test.
    lam_zones, lam_bores = [], []
    if args.laminar_bore and not args.init_data:
        lam_zones, lam_bores = laminar_bores(s, args)

    # Material: constant-property water
    s.setup.materials.database.copy_by_name(type="fluid", name="water-liquid")
    set_first(s, ["setup.materials.fluid['water-liquid'].density.value"], RHO, "density")
    set_first(s, ["setup.materials.fluid['water-liquid'].viscosity.value"], MU, "viscosity")
    for name in s.setup.cell_zone_conditions.fluid.keys():
        set_first(s, [f"setup.cell_zone_conditions.fluid['{name}'].general.material",
                      f"setup.cell_zone_conditions.fluid['{name}'].material"],
                  "water-liquid", f"material {name}")

    # A restart from a solution with another turbulence model starts in that model, reads the
    # data, then switches (Fluent converts omega to epsilon and back when the model changes)
    set_turbulence(s, args.init_turbulence or args.turbulence, args)
    make_laminar(s, solver, lam_zones)

    # Solution methods: coupled, pseudo-transient, second-order (defaults in recent releases)
    try:
        choose(s, "solution.methods.p_v_coupling.flow_scheme", ["coupled"], "p-v coupling")
    except Exception as exc:                           # noqa: BLE001
        LOG.warning("p-v coupling left at default: %s", exc)

    # Pseudo-time stepping. The default global step is set by the domain length scale and the
    # ambient velocity, so the chamber (through-flow ~1e-4 m/s, flushing time ~2 h) advances
    # by only a small fraction of its flushing time per iteration. The 2025 R1 pressure-based
    # coupled solver offers only 'off' and 'global-time-step' (no local step), so the lever is
    # the scale factor on the automatic global step. Requested settings must hold: a silent
    # fallback would make the comparison meaningless.
    if args.pseudo != "default":
        choose(s, "solution.methods.pseudo_time_method.formulation.coupled_solver",
               [f"{args.pseudo}-time-step"], "pseudo-time method")
    pts = "solution.run_calculation.pseudo_time_settings"
    if args.time_scale_factor != 1.0:
        set_first(s, [find_child(s, pts, "scale_factor")], args.time_scale_factor,
                  "pseudo-time scale factor")
    try:
        LOG.info("pseudo-time settings: %s", resolve(s, pts)())
    except Exception as exc:                           # noqa: BLE001 - record only
        LOG.warning("could not read pseudo-time settings: %s", exc)

    # Iteration control is by the report-based criterion below, not residual thresholds
    try:
        choose(s, "solution.monitor.residual.options.criterion_type", ["none"], "residual stop")
    except Exception as exc:                           # noqa: BLE001
        LOG.warning("residual criteria left at default (may stop early): %s", exc)

    # Report definitions
    rd = s.solution.report_definitions
    surf = {
        "mdot_mouth_up": ("surface-massflowrate", None, planes("interior_mouth_upstream")),
        "mdot_mouth_down": ("surface-massflowrate", None, planes("interior_mouth_downstream")),
        "p_mouth_up": ("surface-areaavg", "pressure", planes("interior_mouth_upstream")),
        "p_mouth_down": ("surface-areaavg", "pressure", planes("interior_mouth_downstream")),
        "pt_mouth_up": ("surface-areaavg", "total-pressure", planes("interior_mouth_upstream")),
        "p_inlet": ("surface-areaavg", "pressure", ["inlet"]),
        "mdot_body_mid": ("surface-massflowrate", None, planes("interior_body_mid")),
    }
    for name, (rtype, field, surfaces) in surf.items():
        rd.surface[name] = {}
        rd.surface[name].report_type = rtype
        if field:
            rd.surface[name].field = field
        rd.surface[name].surface_names = surfaces
    walls = [w for w in bc.wall.keys() if w.split(":")[0] in EXTERNAL_WALLS]
    names = list(surf)
    if walls:
        rd.force["drag_x"] = {}
        rd.force["drag_x"].force_vector = [1, 0, 0]
        set_first(s, ["solution.report_definitions.force['drag_x'].zones",
                      "solution.report_definitions.force['drag_x'].thread_names"], walls,
                  "drag zones")
        names.append("drag_x")
        # Full domain only (no y = 0 symmetry zone): the cross-stream force, whose oscillation
        # shows antisymmetric shedding that the half model cannot represent.
        try:
            sym = [z.split(":")[0] for z in bc.symmetry.keys()]
        except AttributeError:
            sym = ["symmetry"]
        if "symmetry" not in sym:
            rd.force["side_y"] = {}
            rd.force["side_y"].force_vector = [0, 1, 0]
            set_first(s, ["solution.report_definitions.force['side_y'].zones",
                          "solution.report_definitions.force['side_y'].thread_names"], walls,
                      "side-force zones")
            names.append("side_y")
    force = {k: walls for k in ("drag_x", "side_y") if k in names}
    if args.init_data:
        # Restart from an earlier solution on the same mesh (e.g. an unconverged pilot)
        s.file.read(file_type="data", file_name=str(Path(args.init_data).resolve()))
        LOG.info("initialised from %s", args.init_data)
        if args.init_turbulence and args.init_turbulence != args.turbulence:
            set_turbulence(s, args.turbulence, args)
            LOG.info("switched turbulence model %s -> %s", args.init_turbulence, args.turbulence)
        if args.laminar_bore:                          # split once the data are in
            before = compute(s, names)
            lam_zones, lam_bores = laminar_bores(s, args)
            reset_reports(s, surf, force)
            check_reports_unchanged(before, compute(s, names))
            make_laminar(s, solver, lam_zones)
    else:
        s.solution.initialization.hybrid_initialize()
    if lam_bores:                                      # stop here, before any iteration
        check_laminar_fraction(s, lam_bores)
    return s, names


def laminar_spec(args):
    """Per --laminar-bore cylinder: the zones to take whole, the zones to split and the reference
    volume, from --laminar-zone and --bore-volume (geometry/laminar_zones.py). ValueError if they
    do not match the cylinders."""
    bores = args.laminar_bore or []
    if not bores:
        raise ValueError("--laminar-zone and --bore-volume need --laminar-bore")
    spec = [dict(index=i, whole=[], split=[], volume=None) for i in range(len(bores))]
    for i, zone, mode in args.laminar_zone or []:
        if not i.isdigit() or int(i) >= len(bores):
            raise ValueError(f"--laminar-zone {i} {zone} {mode}: no laminar bore {i}")
        if mode not in ("whole", "split"):
            raise ValueError(f"--laminar-zone {i} {zone} {mode}: MODE must be whole or split")
        spec[int(i)][mode].append(zone)
    for i, v in args.bore_volume or []:
        if not i.isdigit() or int(i) >= len(bores):
            raise ValueError(f"--bore-volume {i} {v}: no laminar bore {i}")
        spec[int(i)]["volume"] = float(v)
    every = [z for b in spec for z in b["whole"] + b["split"]]
    whole = [z for b in spec for z in b["whole"]]
    for b in spec:
        if b["volume"] is None or not b["volume"] > 0:
            raise ValueError(f"laminar bore {b['index']}: no positive --bore-volume")
        if not b["whole"] and not b["split"]:
            raise ValueError(f"laminar bore {b['index']}: no --laminar-zone")
    if len(set(whole)) != len(whole) or set(whole) & (set(every) - set(whole)):
        raise ValueError(f"a zone taken whole is also listed elsewhere: {whole}")
    return spec


def laminar_bores(s, args):
    """Laminar zones for each --laminar-bore cylinder (x0, x1, r on the x axis).
    Zones listed 'whole' (every cell inside the cylinder, geometry/laminar_zones.py) are taken as
    they are. Every other zone is split by the cylinder's cell register, as before: a zone
    with marked cells gives a new zone '<zone>:<id>'. Previously only the new zones were made
    laminar, so zones wholly inside a cylinder, which Fluent has nothing to split, stayed SST-gamma
. The run stops if a listed zone is
    missing, a zone listed 'split' gives no new zone, or a bore gets no laminar zone.
    Returns the laminar zones and, per bore, its zones and reference volume."""
    czc = s.setup.cell_zone_conditions
    cr = s.solution.cell_registers
    spec = laminar_spec(args)
    taken = [z for b in spec for z in b["whole"]]
    zones0 = list(czc.fluid.keys())
    missing = [z for b in spec for z in b["whole"] + b["split"] if z not in zones0]
    if missing:
        raise RuntimeError(f"laminar zones not in the case: {missing}; fluid zones are {zones0}")
    new, per_bore = [], []
    for (x0, x1, rad), b in zip(args.laminar_bore, spec):
        i = b["index"]
        name = f"bore{i}"
        cr.create(name=name)
        cr[name].type.option.set_state("cylinder")
        cr[name].type.cylinder.set_state({"inside": True, "radius": rad,
                                          "axis_begin": [x0, 0.0, 0.0], "axis_end": [x1, 0.0, 0.0]})
        before = list(czc.fluid.keys())
        for z in before:                               # zones with no marked cells are skipped
            if z not in new and z not in taken:
                s.mesh.modify_zones.sep_cell_zone_mark(cell_zone_name=z, register=name, move_faces=True)
        added = [z for z in czc.fluid.keys() if z not in before]
        parents = {z.rsplit(":", 1)[0] for z in added}
        unsplit = [z for z in b["split"] if z not in parents]
        if unsplit:
            raise RuntimeError(f"laminar bore {i}: zones listed to split gave no new zone: {unsplit} "
                               f"(new zones {added})")
        extra = sorted(parents - set(b["split"]))
        if extra:
            LOG.warning("laminar bore %d: zones split that laminar_zones.py did not list: %s", i, extra)
        if not added and not b["whole"]:
            raise RuntimeError(f"laminar bore {i}: no laminar zone")
        LOG.info("laminar bore %d: x %.4f..%.4f r %.4f, taken whole %s, split %s -> new zones %s",
                 i, x0, x1, rad, b["whole"], sorted(parents), added)
        new += added
        per_bore.append(dict(index=i, zones=b["whole"] + added, volume=b["volume"]))
    for z in new:
        set_first(s, [f"setup.cell_zone_conditions.fluid['{z}'].general.material",
                      f"setup.cell_zone_conditions.fluid['{z}'].material"], "water-liquid", f"material {z}")
    return taken + new, per_bore


def reset_reports(s, surf, force):
    """Put every surface report's surface list and every force report's zones back to their
    setup values. sep_cell_zone_mark (move_faces) appends the new interface surfaces
    'interior--<zone>:<id>' to every surface report already defined,
    so after a split the mouth pressures averaged over 9 or more surfaces instead of 1."""
    rd = s.solution.report_definitions
    for name, (_rtype, _field, surfaces) in surf.items():
        got = list(rd.surface[name].surface_names())
        if got != list(surfaces):
            LOG.info("report %s: surfaces after the split %s; reset to %s", name, got, list(surfaces))
            rd.surface[name].surface_names = list(surfaces)
        now = list(rd.surface[name].surface_names())
        if now != list(surfaces):
            raise RuntimeError(f"report {name}: surfaces {now} after reset, wanted {list(surfaces)}")
    for name, zones in force.items():
        set_first(s, [f"solution.report_definitions.force['{name}'].zones",
                      f"solution.report_definitions.force['{name}'].thread_names"], zones,
                  f"{name} zones (reset after the split)")


def check_reports_unchanged(before, after, rtol=1e-9):
    """Every report must read the same after the split and reset as before the split: no
    iteration lies between, so any change means the report now covers other faces."""
    bad = {k: (before[k], after.get(k)) for k in before
           if after.get(k) is None or abs(after[k] - before[k]) > rtol * max(abs(before[k]), abs(after[k])) + 1e-15}
    for k in before:
        LOG.info("report %-16s before the split %.10g, after %.10g", k, before[k], after.get(k, float("nan")))
    if bad:
        raise RuntimeError(f"reports changed by the zone split (before, after): {bad}; run stopped before iterating")
    LOG.info("all %d reports unchanged by the zone split (relative tolerance %g)", len(before), rtol)


def check_laminar_fraction(s, per_bore):
    """S2's volume test at setup: for each bore, the volume of its laminar
    zones over the cylinder's cell volume from geometry/laminar_zones.py. Stops the run, before any
    iteration, unless every fraction is within 0.99-1.01. The volume is a volume-integral report
    of density over the bore's zones divided by the constant density RHO or, if Fluent refuses
    that, a 'volume-volume' report."""
    rd = s.solution.report_definitions
    names, div = {}, {}
    for b in per_bore:
        name = f"lamvol_bore{b['index']}"
        rd.volume[name] = {}
        errors = []
        for rtype, field, d in (("volume-integral", "density", RHO), ("volume-volume", None, 1.0)):
            try:
                rd.volume[name].report_type = rtype
                if field:
                    rd.volume[name].field = field
                div[name] = d
                LOG.info("laminar volume report %s: %s%s", name, rtype, f" of {field}" if field else "")
                break
            except Exception as exc:                   # noqa: BLE001 - try the next form
                errors.append(f"{rtype}: {exc}")
        else:
            raise RuntimeError(f"no volume report for {name}: " + "; ".join(errors))
        set_first(s, [f"solution.report_definitions.volume['{name}'].cell_zones",
                      f"solution.report_definitions.volume['{name}'].zone_names"], b["zones"],
                  f"zones {name}")
        names[name] = b
    vals = compute(s, list(names))
    bad = []
    for name, b in names.items():
        v = vals[name] / div[name]
        f = v / b["volume"]
        LOG.info("laminar volume fraction bore %d: %.4f (laminar zones %.6e m3 of cylinder cells "
                 "%.6e m3; zones %s)", b["index"], f, v, b["volume"], b["zones"])
        if not 0.99 <= f <= 1.01 or not math.isfinite(f):
            bad.append(f"bore {b['index']}: {f:.4f}")
    if bad:
        raise RuntimeError("laminar volume fraction outside 0.99-1.01, run stopped before "
                           "iterating: " + "; ".join(bad))


def make_laminar(s, solver, zones):
    """Laminar flag and zero turbulent viscosity in each laminar-bore zone."""
    for z in zones:
        set_first(s, [f"setup.cell_zone_conditions.fluid['{z}'].general.laminar",
                      f"setup.cell_zone_conditions.fluid['{z}'].laminar"], True, f"laminar {z}")
        # 'Set turbulent viscosity to zero within laminar zone': in 2025 R1 its settings path
        # (embedded_les.laminar_mut_zero) stays inactive, but the TUI keyword is accepted (probe
        # test). The screen's bore nu_t/nu (field_slice.py) checks that it took effect.
        solver.execute_tui(f"/define/boundary-conditions/set/fluid {z} () laminar-mut-zero? yes q")
        LOG.info("set mu_t = 0 in %s (TUI laminar-mut-zero?)", z)


def set_turbulence(s, model, args):
    """Viscous model and inlet turbulence.

    Baseline sst-gamma: SST k-omega with the intermittency transition model (Menter et al.
    2015), because the 4 mm tube flow is laminar to transitional (Re ~ 30-4000) while the
    external flow is turbulent. Sensitivity models: sst (fully turbulent boundary layers),
    rke (realizable k-epsilon, enhanced wall treatment), sbes (stress-blended eddy simulation
    on SST: LES in the separated wake, RANS in the attached boundary layers), laminar.
    """
    vis = "setup.models.viscous"
    if model == "laminar":
        choose(s, f"{vis}.model", ["laminar"], "viscous model")
        return
    if model == "rke":
        choose(s, f"{vis}.model", ["k-epsilon"], "viscous model")
        choose(s, f"{vis}.k_epsilon_model", ["realizable"], "k-epsilon variant")
        choose_below(s, f"{vis}.near_wall_treatment", ["enhanced"], "near-wall treatment")
    else:
        choose(s, f"{vis}.model", ["k-omega"], "viscous model")
        choose(s, f"{vis}.k_omega_model", ["sst"], "k-omega variant")
        tm = f"{vis}.transition_module"
        if model == "sst-gamma":
            choose(s, tm, ["gamma-transport", "gamma", "intermittency"], "transition model")
        else:
            try:
                choose(s, tm, ["none"], "transition model")
            except Exception as exc:                   # noqa: BLE001 - absent when never enabled
                LOG.info("transition model not reset (%s)", exc)
        if model == "sbes":
            # Fluent 2025 R1: SBES is the hybrid RANS-LES option of SST k-omega, offered for
            # transient runs only (sub-grid model WALE by default)
            choose(s, "setup.general.solver.time", ["unsteady-2nd-order"], "time formulation")
            choose(s, f"{vis}.hybrid_rans_les", ["stress-blended"], "SBES")
            # Bounded central differencing for momentum, as recommended for SBES
            try:
                choose(s, "solution.methods.spatial_discretization.discretization_scheme['mom']",
                       ["bounded-central"], "momentum scheme")
            except Exception as exc:                   # noqa: BLE001
                LOG.warning("momentum scheme left at default: %s", exc)
    try:
        LOG.info("viscous settings: %s", resolve(s, vis)())
    except Exception:                                  # noqa: BLE001 - record only
        pass

    # Inlet turbulence: intensity and viscosity ratio. Ambient inlet: --ti/--tvr (defaults 5 %,
    # 10). Tube-flow inlets (inlet_tube_down, inlet_bore) carry tube flow that is laminar at the
    # baseline Re_tube ~ 280: 2 % and 1.
    bc = s.setup.boundary_conditions
    for vi in bc.velocity_inlet.keys():
        ti, tvr = (args.ti, args.tvr) if vi.split(":")[0] == "inlet" else (0.02, 1.0)
        inlet = f"setup.boundary_conditions.velocity_inlet['{vi}']"
        choose(s, f"{inlet}.turbulence.turbulence_specification",
               ["intensity and viscosity ratio"], f"{vi} turbulence spec")
        set_first(s, [f"{inlet}.turbulence.turbulent_intensity"], ti, f"{vi} TI")
        set_first(s, [f"{inlet}.turbulence.turbulent_viscosity_ratio",
                      f"{inlet}.turbulence.turbulent_viscosity_ratio_real"], tvr, f"{vi} TVR")


def compute(s, names, max_wait=3600):
    """Report definitions by name. A report is read-only, so after a dropped connection the
    call is simply repeated every 30 s until Fluent answers (a drop during this call, after a
    completed step, otherwise ends the run)."""
    t0 = None
    while True:
        try:
            out = {}
            for item in s.solution.report_definitions.compute(report_defs=names):
                for k, v in item.items():
                    out[k] = v[0] if isinstance(v, (list, tuple)) else v
            if t0 is not None:
                LOG.warning("report compute answered again after %.0f s", time.time() - t0)
            return out
        except Exception as exc:                       # noqa: BLE001
            if not session_dropped(exc):
                raise
            if t0 is None:
                t0 = time.time()
                LOG.warning("Fluent dropped the connection during report compute (%s); retrying", exc)
            elif time.time() - t0 > max_wait:
                raise RuntimeError(f"report compute: no answer from Fluent in {max_wait} s") from exc
            time.sleep(30)


def set_velocity(s, U):
    inlet = "setup.boundary_conditions.velocity_inlet['inlet']"
    set_first(s, [f"{inlet}.momentum.velocity_magnitude.value", f"{inlet}.momentum.velocity.value",
                  f"{inlet}.vmag.value", f"{inlet}.vmag"], U, "inlet velocity")


def _set_vec(state, keys, vec):
    """Write a 3-vector into a settings state dict at the first key found (floats or
    {'option': 'value', 'value': x} entries); returns the key used, or None."""
    for k in keys:
        if k in state and isinstance(state[k], (list, tuple)) and len(state[k]) == 3:
            state[k] = [dict(e, value=v) if isinstance(e, dict) else v for e, v in zip(state[k], vec)]
            return k
    return None


def set_inlet_direction(s, U, deg):
    """Inflow at speed U turned deg about the vertical (z) axis, in the x-y plane; deg = 0 restores
    the normal inflow. For the perturbed run with no symmetry plane:
    a short antisymmetric kick tests whether the URANS wake stays on its symmetric branch.
    The settings spelling for components differs between Fluent versions, so the momentum state is
    read, the component entries are located by name and written back, and the result is read again
    and logged. Anything not found raises, so a test job fails at the first kick."""
    vi = "setup.boundary_conditions.velocity_inlet['inlet']"
    if abs(deg) < 1e-12:
        choose_below(s, f"{vi}.momentum", ["Normal to Boundary"], "inlet velocity specification")
        set_velocity(s, U)
        LOG.info("inlet direction restored to normal, U = %.4f m/s: %s", U,
                 ride(resolve(s, f"{vi}.momentum"), "inlet momentum state"))
        return
    a = math.radians(deg)
    vec = [U * math.cos(a), U * math.sin(a), 0.0]
    mom = resolve(s, f"{vi}.momentum")
    method = choose_below(s, f"{vi}.momentum", ["Components", "Magnitude and Direction"],
                          "inlet velocity specification")
    state = ride(mom, "inlet momentum state")
    LOG.info("inlet momentum state before the kick (%s): %s", method, state)
    if "component" in method.lower():
        key = _set_vec(state, [k for k in state if "component" in k], vec)
        if key is None:                                # spelled as separate x/y/z entries
            names = {ax: next((k for k in state if k.startswith(f"{ax}_velocity") or k == f"{ax}_vel"
                               or k.endswith(f"_{ax}")), None) for ax in "xyz"}
            if None in names.values():
                raise RuntimeError(f"inlet kick: no velocity components in {sorted(state)}")
            for ax, v in zip("xyz", vec):
                e = state[names[ax]]
                state[names[ax]] = dict(e, value=v) if isinstance(e, dict) else v
            key = ",".join(names.values())
    else:
        key = _set_vec(state, [k for k in state if "direction" in k], [math.cos(a), math.sin(a), 0.0])
        if key is None:
            raise RuntimeError(f"inlet kick: no flow direction in {sorted(state)}")
        set_velocity(s, U)
    ride(lambda: mom.set_state(state), "setting inlet kick")
    after = ride(mom, "inlet momentum state")
    LOG.info("inlet kick %.1f deg set via %s: %s", deg, key, after)


def solve_to_stationary(s, names, block, max_iter, tol, imbalance_tol=5e-3):
    """Iterate in blocks until the inlet-tube flow rate changes by < tol over three blocks AND
    inflow and outflow of the chamber agree within imbalance_tol.

    The chamber holds ~600 tube volumes, so its interior adjusts far more slowly than the
    external flow; a stationary inlet flow alone is not evidence of a converged solution.
    """
    hist, it, change, imb = [], 0, float("nan"), float("nan")
    while it < max_iter:
        iterate(None, s, block)
        it += block
        vals = compute(s, names)
        up, down = abs(vals["mdot_mouth_up"]), abs(vals["mdot_mouth_down"])
        hist.append([up, down])
        imb = abs(up - down) / max(up, 1e-30)
        if len(hist) >= 4:
            change = max(abs(h[0] - up) for h in hist[-4:-1]) / max(up, 1e-30)
            LOG.info("iter %5d  mdot_in=%.6e  mdot_out=%.6e  mdot_mid=%.6e  change=%.2e  "
                     "imbalance=%.2e", it, up, down, abs(vals["mdot_body_mid"]), change, imb)
            if change < tol and imb < imbalance_tol:
                return vals, it, change, hist
    LOG.warning("not converged after %d iterations (change=%.2e, imbalance=%.2e)", it, change, imb)
    return vals, it, float("nan"), hist


def stage_solve(args, sweep):
    meta = json.loads(Path(args.mesh + ".json").read_text()) if Path(args.mesh + ".json").exists() else {}
    geo_meta = {}
    if args.zones:
        geo_meta = json.loads(Path(args.zones).read_text())
    spec = geo_meta.get("spec", {})
    tube_id = spec.get("tube_id", args.tube_id)
    a_half, mirror = tube_area(geo_meta, tube_id)      # tube cross-section inside the model
    d_body = spec.get("body_od", 0.110)
    tag = Path(args.mesh).name.replace(".msh.h5", "") + f"_{args.turbulence}"
    out = Path(args.out).resolve() / tag
    out.mkdir(parents=True, exist_ok=True)

    solver, s, names = launch_and_setup(args)
    Path(out / "run_metadata.json").write_text(json.dumps(
        metadata(solver, args, {"mesh_meta": meta, "geometry": geo_meta}), indent=2))

    velocities = args.velocities if sweep else [args.velocities[0]]
    rows = []
    for i, U in enumerate(velocities):
        set_velocity(s, U)
        if not sweep:
            iterate(solver, s, 20)
            s.file.write(file_type="case", file_name=str(out / "setup_check.cas.h5"))
            LOG.info("setup stage complete: %s", compute(s, names))
            break
        first = i == 0
        vals, iters, change, hist = solve_to_stationary(
            s, names, block=args.block, max_iter=args.max_iter * (2 if first else 1), tol=args.tol)
        q = mirror * abs(vals["mdot_mouth_up"]) / RHO
        dyn = 0.5 * RHO * U**2
        row = dict(
            U=U, iterations=iters, mdot_change=change, Q_m3s=q, R=q / (mirror * a_half * U),
            Re_tube=(q / (mirror * a_half)) * tube_id * RHO / MU,
            Cp_mouth_up=(vals["p_mouth_up"] - vals["p_inlet"]) / dyn,
            Cpt_mouth_up=(vals["pt_mouth_up"] - vals["p_inlet"]) / dyn,
            Cp_mouth_down=(vals["p_mouth_down"] - vals["p_inlet"]) / dyn,
            Cd_body=mirror * vals["drag_x"] / (dyn * math.pi * d_body**2 / 4),
            mass_imbalance=(abs(vals["mdot_mouth_up"]) - abs(vals["mdot_mouth_down"]))
            / max(abs(vals["mdot_mouth_up"]), 1e-30),
            mdot_body_mid_ratio=abs(vals["mdot_body_mid"]) / max(abs(vals["mdot_mouth_up"]), 1e-30),
        )
        rows.append(row)
        LOG.info("U=%.3f m/s  R=%.4f  Re_tube=%.0f  Cp_up=%.3f  Cp_down=%.3f", U, row["R"],
                 row["Re_tube"], row["Cp_mouth_up"], row["Cp_mouth_down"])
        s.file.write(file_type="case-data", file_name=str(out / f"U{U:.3f}.cas.h5"))
        with open(out / "results.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        (out / f"U{U:.3f}_mdot_history.json").write_text(json.dumps(hist))
    solver.exit()


def add_zone_vmax(s):
    """One volume-max velocity report per cell zone, to locate where a start-up blows up."""
    rd, names = s.solution.report_definitions, []
    for i, zone in enumerate(s.setup.cell_zone_conditions.fluid.keys()):
        name = f"vmax_{zone.replace(':', '_').replace('-', '_')}"
        try:
            rd.volume[name] = {}
            rd.volume[name].report_type = "volume-max"
            rd.volume[name].field = "velocity-magnitude"
            set_first(s, [f"solution.report_definitions.volume['{name}'].cell_zones",
                          f"solution.report_definitions.volume['{name}'].zone_names"], [zone],
                      f"zone {zone}")
            names.append(name)
        except Exception as exc:                       # noqa: BLE001 - diagnostics only
            LOG.warning("no velocity report for zone %s: %s", zone, exc)
    return names


def log_zone_vmax(s, vnames, where):
    if not vnames:
        return
    try:
        v = compute(s, vnames)
        top = sorted(v.items(), key=lambda kv: -abs(kv[1]) if math.isfinite(kv[1]) else -1e300)
        LOG.info("max |U| by zone at %s: %s", where,
                 ", ".join(f"{k[5:]}={val:.3g}" for k, val in top[:6]))
    except Exception as exc:                           # noqa: BLE001 - diagnostics only
        LOG.warning("zone velocity report failed: %s", exc)


def check_finite(vals, where, R=None, s=None, dump=None):
    """Stop a diverged run (NaN or overflow) instead of time-stepping it to the end.

    Diverged runs can go to NaN in their first minutes and then hold their licensed cores for
    hours; exit code 4 is not a licence failure, so the job script
    neither waits nor resubmits. The diverged field is written to dump for inspection.
    """
    bad = [k for k, v in vals.items() if not math.isfinite(v) or abs(v) > 1e6]
    if R is not None and (not math.isfinite(R) or R > 10):
        bad.append("R")
    if bad:
        LOG.error("diverged at %s (%s): %s", where, ", ".join(bad), vals)
        if s is not None and dump:
            try:
                s.file.write(file_type="case-data", file_name=str(dump))
            except Exception as exc:                   # noqa: BLE001
                LOG.warning("could not write %s: %s", dump, exc)
        # Not sys.exit: a run exiting that way can hang for hours at "Shutting
        # down Server..." with its cores held. Slurm removes the Fluent processes at job end.
        logging.shutdown()
        sys.stdout.flush()
        os._exit(4)


TC = "solution.run_calculation.transient_controls"


# Riding out a dropped connection. The Fluent gRPC server pings this client every 10 s
# and closes the connection if a ping goes unanswered for ~20 s ("keepalive watchdog timeout";
# the string is in Fluent's libGrpcServer.so, not in the client's gRPC). That happens when the
# client process is held up on a busy node. Fluent
# itself carries on with the command it was running, and the client's channel reconnects on the
# next call, so instead of exiting the client waits for Fluent to come back, checks how far the
# flow time got and runs any steps the command did not complete. Tested by freezing the client
# with SIGSTOP in the middle of a call (ka_probe2.py).
SESSION_DROP = ("keepalive watchdog timeout", "Socket closed", "Connection reset", "UNAVAILABLE",
                "failed to connect", "Broken pipe")


def session_dropped(exc):
    return isinstance(exc, RuntimeError) and any(k in str(exc) for k in SESSION_DROP)


def rp_var(solver, name):
    return float(scheme(solver, f"(rpgetvar '{name})"))


def scheme(solver, expr):
    sc = solver.scheme
    for f in ("scheme_eval", "eval", "string_eval"):
        if hasattr(sc, f):
            return getattr(sc, f)(expr)
    raise AttributeError("no scheme evaluator on the session")


def wait_for_session(solver, exc, what, max_wait=6 * 3600):
    """Block until the Fluent session answers again; return the flow time (None if unset)."""
    LOG.warning("Fluent dropped the connection during %s (%s); Fluent keeps running the command, "
                "waiting for it", what, exc)
    t0, last = time.time(), None
    while time.time() - t0 < max_wait:
        try:
            scheme(solver, "(+ 1 1)")
        except Exception as err:                       # noqa: BLE001
            last = err
            time.sleep(30)
            continue
        try:
            ft = rp_var(solver, "flow-time")
        except Exception:                              # noqa: BLE001
            ft = None
        LOG.warning("session back after %.0f s, flow time %s", time.time() - t0, ft)
        return ft
    raise RuntimeError(f"Fluent session did not come back in {max_wait} s ({last})") from exc


def wait_idle(cmd, what, max_wait=3600):
    """Wait until a settings command is active again, i.e. Fluent has finished its command."""
    t0 = time.time()
    while time.time() - t0 < max_wait:
        try:
            if cmd.is_active():
                LOG.warning("Fluent idle again (%s active) after %.0f s", what, time.time() - t0)
                return
        except Exception:                              # noqa: BLE001
            pass
        time.sleep(15)
    raise RuntimeError(f"Fluent still busy after {max_wait} s ({what} inactive)")


def write_case_data(solver, s, path):
    """Write case and data; if the connection drops during the write, wait for Fluent to finish
    it (the command keeps running on the server) rather than lose the job."""
    try:
        s.file.write(file_type="case-data", file_name=str(path))
    except Exception as exc:                           # noqa: BLE001
        if not session_dropped(exc):
            raise
        wait_for_session(solver, exc, f"write {Path(path).name}")
        wait_idle(s.file.write, "write")
        LOG.warning("write of %s finished after the drop", Path(path).name)


def launch_and_setup(args, tries=3):
    """Launch Fluent and set the case up. A dropped connection before any iteration loses
    nothing, so Fluent is killed and the launch and setup repeated ("Socket closed" in a settings .keys() during report setup)."""
    for k in range(tries):
        solver = launch("solver", args.cores)
        try:
            s, names = setup_solver(solver, args)
            return solver, s, names
        except Exception as exc:                       # noqa: BLE001
            if not session_dropped(exc) or k == tries - 1:
                raise
            LOG.warning("Fluent dropped the connection during setup (%s); killing it and starting "
                        "again (attempt %d of %d)", exc, k + 2, tries)
            try:
                solver.force_exit()
            except Exception as err:                   # noqa: BLE001
                LOG.warning("force_exit failed: %s", err)
            time.sleep(60)


def read_mesh(solver, s, mesh):
    """Read the mesh; if the connection drops during the read (seen with
    an 870 MB mesh), wait for Fluent to finish it instead of retrying on a dead connection."""
    try:
        s.file.read_mesh(file_name=mesh)
        return
    except Exception as exc:                           # noqa: BLE001
        if not session_dropped(exc):
            s.file.read(file_type="mesh", file_name=mesh)
            return
        wait_for_session(solver, exc, "read_mesh")
        prev, t0 = None, time.time()
        while time.time() - t0 < 3600:                 # the zone list settles when the read ends
            try:
                names = sorted(s.setup.cell_zone_conditions.fluid.get_object_names())
            except Exception:                          # noqa: BLE001
                names = None
            if names and names == prev:
                wait_idle(s.file.read_mesh, "read_mesh")
                LOG.warning("mesh read completed after the drop: %d fluid zones", len(names))
                return
            prev = names
            time.sleep(60)
        raise RuntimeError("mesh not loaded after the dropped read") from exc


def iterate(solver, s, n):
    """n steady iterations. If the connection drops during the block (at iteration ~250
    of --steady-init 800), Fluent finishes the block on its
    own: wait for the session and for the iterate command to be free again, then carry on.
    solver may be None (callers that only hold the settings root); the wait then polls the
    command alone. Shared with run_f4_case.py."""
    try:
        s.solution.run_calculation.iterate(iter_count=n)
    except Exception as exc:                           # noqa: BLE001
        if not session_dropped(exc):
            raise
        if solver is not None:
            wait_for_session(solver, exc, f"iterate({n})")
        else:
            LOG.warning("Fluent dropped the connection during iterate(%d) (%s); waiting for it", n, exc)
        wait_idle(s.solution.run_calculation.iterate, "iterate", max_wait=6 * 3600)
        LOG.warning("iterate(%d) finished after the drop", n)


_LAST_FT = [None]          # flow time at the end of the last dual_time() call (None if unknown)


def _flow_time(solver, max_wait=600):
    """Flow time, riding out a connection drop; None if it cannot be read."""
    try:
        return ride(lambda: rp_var(solver, "flow-time"), "flow-time read", max_wait)
    except Exception as exc:                           # noqa: BLE001
        LOG.warning("flow time unreadable (%s)", exc)
        return None


def dual_time(solver, s, n, inner, _retry=0):
    """n time steps of dual_time_iterate. If the connection drops during the command, Fluent keeps
    running it: wait for the session and the command, count the steps done from the flow time and
    run the rest. The flow time before the call is read with retries; if it is still unknown (the
    pre-call read can fail with the same drop that then hits the iterate), the
    flow time at the end of the previous call stands in for it, provided the current flow time is
    consistent with it. With no usable reference the command is taken as complete, since Fluent
    finishes a command after a drop in every case seen."""
    t0 = _flow_time(solver)
    try:
        dt = ride(lambda: rp_var(solver, "physical-time-step"), "time-step read", 600)
    except Exception:                                  # noqa: BLE001
        dt = None
    try:
        s.solution.run_calculation.dual_time_iterate(time_step_count=n, max_iter_per_step=inner)
    except (TypeError, AttributeError) as exc:         # older/newer settings API spelling
        LOG.warning("settings dual_time_iterate failed (%s); using the TUI", exc)
        solver.tui.solve.dual_time_iterate(n, inner)
    except RuntimeError as exc:
        if not session_dropped(exc) or _retry >= 5:
            raise
        ft = wait_for_session(solver, exc, f"dual_time_iterate({n})")
        if not dt:
            try:
                dt = ride(lambda: rp_var(solver, "physical-time-step"), "time-step read", 600)
            except Exception:                          # noqa: BLE001
                dt = None
        if ft is None:
            ft = _flow_time(solver)
        if t0 is None and dt and ft is not None and _LAST_FT[0] is not None                 and _LAST_FT[0] - 0.5 * dt <= ft <= _LAST_FT[0] + (n + 0.5) * dt:
            t0 = _LAST_FT[0]
            LOG.warning("flow time before the call unknown; using %.6g s, the flow time at the end "
                        "of the previous call (now %.6g s)", t0, ft)
        if t0 is None or not dt or ft is None:
            wait_idle(s.solution.run_calculation.dual_time_iterate, "dual_time_iterate")
            _LAST_FT[0] = _flow_time(solver)
            LOG.warning("flow time before the call unknown and no usable reference (dt %s, now %s); "
                        "Fluent finishes a command after a drop, so the %d steps are taken as done",
                        dt, _LAST_FT[0], n)
            return
        # Fluent answers between time steps while the command is still running (in a test it
        # was back after 4 s at 2 of 8 steps), so wait until the flow time reaches the target or
        # stops advancing for 10 min before deciding how many steps are left.
        target, last, still = t0 + n * dt, ft, time.time()
        while ft < target - 0.5 * dt and time.time() - still < 600:
            time.sleep(30)
            try:
                ft = rp_var(solver, "flow-time")
            except Exception:                          # noqa: BLE001
                continue
            if ft > last + 1e-9 * max(1.0, abs(last)):
                last, still = ft, time.time()
        # Fluent sets the flow time at the start of a step, so wait until the command object is
        # active again, i.e. the last step has finished (a fixed 60 s wait was not enough).
        wait_idle(s.solution.run_calculation.dual_time_iterate, "dual_time_iterate")
        ft = ride(lambda: rp_var(solver, "flow-time"), "flow-time read")
        _LAST_FT[0] = ft
        done = int(round((ft - t0) / dt))
        LOG.warning("the command completed %d of %d time steps", done, n)
        if done < n:
            dual_time(solver, s, n - done, inner, _retry + 1)
        return
    _LAST_FT[0] = t0 + n * dt if (t0 is not None and dt) else None


def start_unsteady(solver, s, args, names, vnames, dump, dt):
    """Steady pre-iterations, switch to unsteady, then a short-step ramp up to dt.

    Returns the simulated time used by the ramp, so that every time-dependent run starts the
    same way.
    """
    if not args.init_data and args.steady_init > 0:
        # Hybrid initialisation is a potential-flow field that drives a frictionless, far too
        # large flow through the tubes (R ~ 2.7 at the first step in the flume case, which then
        # diverged). Steady iterations first let wall friction set the tube flow.
        U = args.velocities[0]
        if args.steady_ramp_from and args.steady_ramp_from < U:
            # Velocity continuation: the Elliott full domain diverged within 6-100 steady
            # iterations started at 1.25 m/s (twice, also with the pseudo-time step x0.25), while
            # the same mesh converged at 0.3 m/s. Start there and step the inlet up geometrically,
            # re-initialising only once at the lowest speed; then the usual steady_init at U.
            n_st = args.steady_ramp_stages
            for k in range(n_st):
                u = args.steady_ramp_from * (U / args.steady_ramp_from) ** (k / n_st)
                set_velocity(s, u)
                if k == 0:
                    s.solution.initialization.hybrid_initialize()
                done = 0
                while done < args.steady_ramp_iters:
                    n = min(100, args.steady_ramp_iters - done)
                    iterate(solver, s, n)
                    done += n
                    where = f"steady velocity ramp {u:.3f} m/s iteration {done}"
                    log_zone_vmax(s, vnames, where)
                    check_finite(compute(s, names), where, s=s, dump=dump)
                LOG.info("steady velocity ramp at %.3f m/s done: %s", u, compute(s, names))
            set_velocity(s, U)
        done = 0
        while done < args.steady_init:                 # in blocks, to stop early on divergence
            n = min(100, args.steady_init - done)
            iterate(solver, s, n)
            done += n
            log_zone_vmax(s, vnames, f"steady iteration {done}")
            check_finite(compute(s, names), f"steady iteration {done}", s=s, dump=dump)
        LOG.info("steady pre-iterations (%d): %s", args.steady_init, compute(s, names))
        # Post-steady checkpoint: a run lost in the ramp can be
        # resumed from here with --init-data U<U>_steady.dat.h5 --ramp-after-init, instead of
        # repeating the steady phase
        if dump is not None and dump.name.endswith("_diverged.cas.h5"):
            ck = dump.with_name(dump.name.replace("_diverged.cas.h5", "_steady.cas.h5"))
            write_case_data(solver, s, ck)
            LOG.info("wrote post-steady checkpoint %s", ck)
    choose(s, "setup.general.solver.time", ["unsteady-2nd-order"], "time formulation")
    set_first(s, [f"{TC}.time_step_size"], dt, "time step")
    set_first(s, [f"{TC}.max_iter_per_time_step"], args.inner, "inner iterations")
    t = 0.0
    switched = bool(args.init_turbulence) and args.init_turbulence != args.turbulence
    if (not args.init_data or switched or args.ramp_after_init) and args.ramp_steps > 0:
        # Gentle start from an initialised field, or after a change of turbulence model (20
        # steady iterations straight after switching SST-gamma -> realizable k-epsilon raised
        # the downstream mouth flow 30-fold): short steps first
        dt_ramp = dt / args.ramp_factor
        set_first(s, [f"{TC}.time_step_size"], dt_ramp, "ramp time step")
        done = 0
        while done < args.ramp_steps:                  # in blocks of 10, to locate a blow-up
            n = min(10, args.ramp_steps - done)
            dual_time(solver, s, n, 2 * args.inner)
            done += n
            log_zone_vmax(s, vnames, f"ramp step {done}")
            check_finite(compute(s, names), f"ramp step {done}", s=s, dump=dump)
        t = args.ramp_steps * dt_ramp
        LOG.info("ramp (%d steps of %.4g s) done: %s", args.ramp_steps, dt_ramp, compute(s, names))
        # Then double the step towards dt, 10 steps at each size: the flume run at 0.585 m/s
        # survived the ramp at dt/20 but diverged within 0.3 s of a single jump to dt
        dt_k = 2 * dt_ramp
        while dt_k < dt * (1 - 1e-9):
            set_first(s, [f"{TC}.time_step_size"], dt_k, "ramp time step")
            dual_time(solver, s, 10, 2 * args.inner)
            t += 10 * dt_k
            log_zone_vmax(s, vnames, f"ramp dt = {dt_k:.4g} s")
            check_finite(compute(s, names), f"ramp dt = {dt_k:.4g} s", s=s, dump=dump)
            dt_k *= 2
        set_first(s, [f"{TC}.time_step_size"], dt, "time step")
    return t


def stage_transient(args):
    """Time-dependent (URANS) run of the coupled model (external flow, bores and chamber in one
    domain) at one ambient velocity: every time step is iterated to the inner-iteration
    limit, the wake behind the blunt rear is resolved in time, and the reported flow rates are
    time series whose second half is averaged.
    """
    U = args.velocities[0]
    tag = Path(args.mesh).name.replace(".msh.h5", "") + f"_{args.turbulence}"
    out = Path(args.out).resolve() / tag
    out.mkdir(parents=True, exist_ok=True)
    guard_out(out, f"U{U:.3f}", args)
    geo_meta = json.loads(Path(args.zones).read_text()) if args.zones else {}
    tube_id = geo_meta.get("spec", {}).get("tube_id", args.tube_id)
    a_half, _ = tube_area(geo_meta, tube_id)           # tube cross-section inside the model

    solver, s, names = launch_and_setup(args)
    set_velocity(s, U)
    if args.inlet_profile:
        # Developed flume inflow from the precursor (stage_precursor) at this bulk velocity
        with open(args.inlet_profile) as fh:
            rows = list(csv.DictReader(fh))
        cols = {k: [float(r[k]) for r in rows] for k in rows[0]}
        apply_inlet_profile(solver, s, cols, x_of_inlet(s), "inflow", out / f"U{U:.3f}_inflow.prof")
        s.solution.initialization.hybrid_initialize()
    vnames = add_zone_vmax(s)
    dump = out / f"U{U:.3f}_diverged.cas.h5"
    Path(out / f"U{U:.3f}_run_metadata.json").write_text(json.dumps(
        metadata(solver, args, {"geometry": geo_meta}), indent=2))

    chunk = max(1, round(args.report_every / args.dt))
    rows, t_ckpt = [], 0.0
    t = start_unsteady(solver, s, args, names, vnames, dump, args.dt)
    kick_on = bool(args.kick)
    if kick_on:                                       # perturbed run (wake-symmetry test)
        set_inlet_direction(s, U, args.kick[1])
    while t < args.sim_time - 1e-9:
        dual_time(solver, s, chunk, args.inner)
        t += chunk * args.dt
        if kick_on and t >= args.kick[0] - 1e-9:
            set_inlet_direction(s, U, 0.0)
            kick_on = False
        vals = compute(s, names)
        up, down = abs(vals["mdot_mouth_up"]), abs(vals["mdot_mouth_down"])
        rows.append(dict(t=t, R=up / RHO / (a_half * U), **vals))
        LOG.info("t=%7.3f s  mdot_in=%.6e  mdot_out=%.6e  mdot_mid=%.6e  R=%.4f", t, up, down,
                 abs(vals["mdot_body_mid"]), rows[-1]["R"])
        check_finite(vals, f"t = {t:.3f} s", rows[-1]["R"], s=s, dump=dump)
        with open(out / f"U{U:.3f}_timeseries.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        if t - t_ckpt >= args.checkpoint - 1e-9:
            write_case_data(solver, s, out / f"U{U:.3f}_ckpt.cas.h5")
            t_ckpt = t
    late = [r for r in rows if r["t"] > args.sim_time / 2]
    mean = {k: sum(r[k] for r in late) / len(late) for k in ("R", "mdot_mouth_up",
                                                             "mdot_mouth_down", "mdot_body_mid")}
    LOG.info("time-averaged over t > %.1f s: %s", args.sim_time / 2, mean)
    (out / f"U{U:.3f}_summary.json").write_text(json.dumps(
        {"U": U, "dt": args.dt, "inner": args.inner, "sim_time": args.sim_time,
         "average_from": args.sim_time / 2, "mean": mean}, indent=2))
    write_case_data(solver, s, out / f"U{U:.3f}_final.cas.h5")
    solver.exit()


# ---------------------------------------------------------------------------------------------
# flume inflow: precursor duct and inlet profiles (Goharrokhi flume)
# ---------------------------------------------------------------------------------------------
# Fluent field names written to and read from profiles, and the column names used here
PROFILE_FIELDS = [("x-velocity", "u"), ("turb-kinetic-energy", "k"),
                  ("specific-diss-rate", "omega"), ("intermittency", "gamma")]


def read_fluent_profile(path):
    """{column: [floats]} from a Fluent point profile ((name point n) (x ...) (y ...) ...)."""
    text = Path(path).read_text()
    body = text[text.index(")") + 1:]                  # after the (name point n) header
    cols = {}
    for name, vals in re.findall(r"\(\s*([\w\-]+)\s+([^()]*)\)", body):
        cols[name] = [float(v) for v in vals.split()]
    return cols


def write_fluent_profile(path, name, cols):
    n = len(next(iter(cols.values())))
    lines = [f"(({name} point {n})"]
    for k, v in cols.items():
        lines.append(f"({k}")
        lines += [f"{a:.9e}" for a in v]
        lines.append(")")
    lines.append(")")
    Path(path).write_text("\n".join(lines) + "\n")


def surface_profile(solver, surface, path):
    """Write a Fluent profile of the PROFILE_FIELDS on a boundary and read it back as columns."""
    path = Path(path).resolve()
    path.unlink(missing_ok=True)
    fields = " ".join(f for f, _ in PROFILE_FIELDS)
    solver.execute_tui(f'/file/write-profile "{path}" {surface} () {fields} ()')
    if not path.exists():
        raise RuntimeError(f"write-profile wrote nothing to {path}")
    return read_fluent_profile(path)


def apply_inlet_profile(solver, s, cols, x_in, pname, path):
    """Map profile columns (y, z, u, k, omega, gamma) onto the domain inlet at x = x_in.

    Velocity magnitude normal to the inlet, k, omega and intermittency all come from the profile
    (the developed RANS duct flow has no cross-stream velocity: an isotropic eddy viscosity drives
    no secondary currents).
    """
    out = {"x": [x_in] * len(cols["y"]), "y": cols["y"], "z": cols["z"]}
    out.update({c: cols[c] for _, c in PROFILE_FIELDS})
    write_fluent_profile(path, pname, out)
    solver.execute_tui(f'/file/read-profile "{Path(path).resolve()}"')
    vi = "setup.boundary_conditions.velocity_inlet['inlet']"
    choose(s, f"{vi}.turbulence.turbulence_specification", ["k and omega"], "inlet turbulence spec")
    def active_child(base, prefer, contains):
        # exact names first; skip inactive ones (e.g. laminar_kinetic_energy, which a plain
        # "kinetic" search found first and which cannot be set)
        obj = resolve(s, base)
        names = list(getattr(obj, "child_names", []))
        for n in prefer + [c for c in names if contains in c and c not in prefer]:
            if n not in names:
                continue
            try:
                child = getattr(obj, n)
                if hasattr(child, "is_active") and not child.is_active():
                    continue
            except Exception:                          # noqa: BLE001 - inactive branch
                continue
            return f"{base}.{n}"
        raise RuntimeError(f"no active setting like {prefer} below {base}: {names}")

    tb = f"{vi}.turbulence"
    targets = [(f"{vi}.momentum.velocity_magnitude", "u"),
               (active_child(tb, ["turbulent_kinetic_energy"], "turbulent_kinetic"), "k"),
               (active_child(tb, ["specific_dissipation_rate"], "specific_dissipation"), "omega")]
    try:
        targets.append((active_child(tb, ["intermittency"], "intermittency"), "gamma"))
    except RuntimeError:                               # no transition model
        pass
    for path_, field in targets:
        obj = resolve(s, path_)
        obj.set_state({"option": "profile", "profile_name": pname, "field_name": field})
        state = obj()
        if state.get("option") != "profile" or state.get("profile_name") != pname:
            raise RuntimeError(f"inlet profile not applied at {path_}: {state}")
        LOG.info("inlet %s <- profile %s field %s", path_.rsplit(".", 1)[-1], pname, field)


def x_of_inlet(s):
    rd = s.solution.report_definitions
    rd.surface["x_inlet"] = {}
    rd.surface["x_inlet"].report_type = "surface-areaavg"
    rd.surface["x_inlet"].field = "x-coordinate"
    rd.surface["x_inlet"].surface_names = ["inlet"]
    return compute(s, ["x_inlet"])["x_inlet"]


def stage_precursor(args):
    """Developed flow of the empty flume, by recycling a duct's outlet profile to its inlet.

    A streamwise-periodic precursor without a periodic mesh: each pass maps the outlet profile
    (u, k, omega, intermittency) onto the inlet, with u scaled so the inlet carries the bulk
    velocity, and iterates the duct; it stops when the outlet profile changes by less than
    --precursor-tol (max |du| / U) between passes. Writes U<U>_profile.csv per velocity.
    """
    geo_meta = json.loads(Path(args.zones).read_text())
    a_in = geo_meta["zones"]["inlet"]["area_m2"]
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    solver, s, _ = launch_and_setup(args)
    rd = s.solution.report_definitions
    rd.surface["mdot_in"] = {}
    rd.surface["mdot_in"].report_type = "surface-massflowrate"
    rd.surface["mdot_in"].surface_names = ["inlet"]
    x_in = x_of_inlet(s)
    LOG.info("precursor duct: inlet area %.4f m2 at x = %.4f m", a_in, x_in)
    vi = "setup.boundary_conditions.velocity_inlet['inlet']"
    for U in args.velocities:
        # constant inflow (TI/TVR as set_turbulence) for the first pass
        resolve(s, f"{vi}.momentum.velocity_magnitude").set_state({"option": "value", "value": U})
        choose(s, f"{vi}.turbulence.turbulence_specification", ["intensity and viscosity ratio"],
               "inlet turbulence spec")
        s.solution.initialization.hybrid_initialize()
        iterate(solver, s, args.precursor_iters)
        prev, scale, change, n = None, 1.0, float("inf"), 0
        for n in range(1, args.precursor_cycles + 1):
            cols = surface_profile(solver, "outlet", out / f"U{U:.3f}_outlet.prof")
            mdot = abs(compute(s, ["mdot_in"])["mdot_in"])
            scale = RHO * U * a_in / mdot                # bulk velocity held at U
            cols["u"] = [v * scale for v in cols.pop("x-velocity")]
            for f, c in PROFILE_FIELDS[1:]:
                if f in cols:
                    cols[c] = cols.pop(f)
            cols.setdefault("gamma", [1.0] * len(cols["u"]))
            change = (max(abs(a - b) for a, b in zip(cols["u"], prev["u"])) / U
                      if prev and len(prev["u"]) == len(cols["u"]) else float("inf"))
            LOG.info("U %.3f pass %d: max |du|/U %.2e, bulk scale %.5f, u max %.4f", U, n, change,
                     scale, max(cols["u"]))
            check_finite({"u_max": max(cols["u"]), "k_max": max(cols["k"])}, f"precursor pass {n}")
            prev = cols
            if change < args.precursor_tol:
                break
            apply_inlet_profile(solver, s, cols, x_in, f"rec{n}", out / f"U{U:.3f}_recycled.prof")
            iterate(solver, s, args.precursor_iters)
        else:
            LOG.warning("U %.3f: precursor not converged after %d passes", U, args.precursor_cycles)
        keys = ["y", "z", "u", "k", "omega", "gamma"]
        with open(out / f"U{U:.3f}_profile.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(keys)
            w.writerows(zip(*[prev[k] for k in keys]))
        # centreline velocity at the sampler axis height (y = 0, z = 0) over the bulk velocity
        i = min(range(len(prev["y"])), key=lambda j: prev["y"][j] ** 2 + prev["z"][j] ** 2)
        LOG.info("U %.3f developed: u(axis) / U = %.4f at y %.4f z %.4f; passes %d", U,
                 prev["u"][i] / U, prev["y"][i], prev["z"][i], n)
        Path(out / f"U{U:.3f}_precursor.json").write_text(json.dumps(dict(
            U=U, passes=n, change=change, bulk_scale=scale, u_axis_over_U=prev["u"][i] / U,
            n_points=len(prev["u"])), indent=2))
    solver.exit()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("stage", choices=["mesh", "setup", "sweep", "transient", "precursor"])
    ap.add_argument("--geometry", help="boundary mesh from tims_geometry.py (mesh stage)")
    ap.add_argument("--mesh", help="volume mesh .msh.h5 from the mesh stage (setup/sweep)")
    ap.add_argument("--zones", help="<case>_zones.json from tims_geometry.py (for tube/body size)")
    ap.add_argument("--tube-id", type=float, default=0.004)
    ap.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_NTASKS", 8)))
    ap.add_argument("--out", default="out")
    models = ["sst-gamma", "sst", "rke", "sbes", "laminar"]
    ap.add_argument("--turbulence", choices=models, default="sst-gamma")
    ap.add_argument("--init-turbulence", choices=models,
                    help="turbulence model of the --init-data solution, when it differs")
    ap.add_argument("--ti", type=float, default=0.05, help="inlet turbulence intensity")
    ap.add_argument("--tvr", type=float, default=10.0, help="inlet turbulent viscosity ratio")
    ap.add_argument("--kick", type=float, nargs=2, metavar=("SECONDS", "DEG"), default=None,
                    help="transient: turn the inflow DEG about the vertical for the first SECONDS of the time stepping, then restore it (wake-symmetry test, full domain)")
    ap.add_argument("--overwrite", action="store_true",
                    help="allow a transient run to replace an existing final/ckpt/steady/summary in "
                         "its output folder (refused by default since v3g, 23 Nov)")
    ap.add_argument("--laminar-bore", type=float, nargs=3, metavar=("X0", "X1", "R"), action="append",
                    help="make the cells within radius R of the x axis, X0 < x < X1, a laminar "
                         "zone with zero turbulent viscosity; repeat for each bore")
    ap.add_argument("--laminar-zone", nargs=3, metavar=("I", "ZONE", "MODE"), action="append",
                    help="cell zone ZONE belongs to laminar bore I (0-based, in --laminar-bore "
                         "order); MODE 'whole' (every cell inside the cylinder: taken whole) or "
                         "'split' (straddles it: its inside cells are split off). From "
                         "geometry/laminar_zones.py; required with --laminar-bore (24 Nov)")
    ap.add_argument("--bore-volume", nargs=2, metavar=("I", "V"), action="append",
                    help="volume (m3) of the cells inside laminar-bore cylinder I, from "
                         "geometry/laminar_zones.py; the run stops before iterating if the "
                         "laminar zones of bore I hold less than 99 %% (or more than 101 %%) of it")
    ap.add_argument("--velocities", type=float, nargs="+", default=DEFAULT_U)
    ap.add_argument("--block", type=int, default=100)
    ap.add_argument("--max-iter", type=int, default=3000)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--pseudo", choices=["default", "global", "local"], default="default",
                    help="coupled-solver pseudo-time method (default: leave Fluent's setting)")
    ap.add_argument("--time-scale-factor", type=float, default=1.0,
                    help="multiplier on Fluent's automatic pseudo-time step")
    ap.add_argument("--init-data", help="data file (.dat.h5) on the same mesh to start from")
    ap.add_argument("--dt", type=float, default=0.02, help="time step (s), transient stage")
    ap.add_argument("--inner", type=int, default=10, help="iterations per time step")
    ap.add_argument("--sim-time", type=float, default=20.0, help="simulated time (s)")
    ap.add_argument("--report-every", type=float, default=0.2, help="report interval (s)")
    ap.add_argument("--checkpoint", type=float, default=2.0, help="case-data write interval (s)")
    ap.add_argument("--steady-init", type=int, default=1500,
                    help="steady iterations before time stepping when not restarting (transient)")
    ap.add_argument("--steady-ramp-from", type=float, default=None,
                    help="start the steady pre-iterations at this inlet velocity (m/s) and step it "
                         "up geometrically to the run velocity (transient, not restarting)")
    ap.add_argument("--steady-ramp-stages", type=int, default=3,
                    help="velocity stages below the run velocity in the steady ramp")
    ap.add_argument("--steady-ramp-iters", type=int, default=200,
                    help="steady iterations per stage of the velocity ramp")
    ap.add_argument("--ramp-after-init", action="store_true",
                    help="run the short-step ramp even when starting from --init-data (resume from a "
                         "post-steady checkpoint, U<U>_steady.dat.h5)")
    ap.add_argument("--ramp-steps", type=int, default=50,
                    help="initial time steps at dt/ramp-factor when not restarting (transient)")
    ap.add_argument("--ramp-factor", type=float, default=5.0,
                    help="ramp time step = dt / ramp-factor (transient)")
    ap.add_argument("--inlet-profile",
                    help="precursor U<U>_profile.csv to impose on the inlet (transient)")
    ap.add_argument("--precursor-iters", type=int, default=150,
                    help="steady iterations per recycling pass (precursor)")
    ap.add_argument("--precursor-cycles", type=int, default=60,
                    help="maximum recycling passes (precursor)")
    ap.add_argument("--precursor-tol", type=float, default=2e-3,
                    help="stop when the outlet u changes by less than this x U between passes")
    ap.add_argument("--solid-zones", nargs="*", default=[],
                    help="cell zones (region names) to retype as solid when read in as fluid, "
                         "e.g. wall_body in the full domain")
    args = ap.parse_args()
    if args.laminar_bore or args.laminar_zone or args.bore_volume:
        try:
            laminar_spec(args)                         # argument errors before Fluent starts
        except ValueError as exc:
            ap.error(str(exc))
    Path(args.out).mkdir(parents=True, exist_ok=True)
    # One log per stage and mesh, so concurrent array tasks do not share a file
    tag = (Path(args.geometry).parent.name if args.stage == "mesh"
           else Path(args.mesh).name.replace(".msh.h5", "") + f"_{args.turbulence}")
    if args.stage == "transient":
        tag += f"_U{args.velocities[0]:.3f}"           # one job per velocity, run concurrently
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(Path(args.out) / f"{args.stage}_{tag}.log")])
    if args.stage == "mesh":
        stage_mesh(args)
    elif args.stage == "precursor":
        stage_precursor(args)
    elif args.stage == "transient":
        stage_transient(args)
    else:
        stage_solve(args, sweep=args.stage == "sweep")
    # Python's normal exit can wait for hours on the PyFluent connection after Fluent has shut
    # down (a mesh job sat 25 min at "Services Shutdown...Done.")
    logging.shutdown()
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
