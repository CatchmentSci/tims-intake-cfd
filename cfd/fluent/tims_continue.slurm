#!/bin/bash
# Continue a transient run that hit its time limit, from its last checkpoint (the 1.5 m/s
# fine-mesh runs reach only t ~ 3-4 of 7 s inside 48 h).
#
#   sbatch --dependency=afterany:<job> --export=ALL,CASE=perks2014,SCALE=0.64,VEL=1.500,\
#       GEOVER=v5,TAG=urans_v5,DT=0.004,SIM=7 tims_continue.slurm
#
# Run at job start, so it sees the checkpoint the first job actually left:
#  - if the first part finished (summary written), there is nothing to do;
#  - otherwise it finds the time of the last checkpoint from the first part's log, copies that
#    field to out/${TAG}_cont/<tag>/U<vel>_start.dat.h5 and runs tims_job.slurm with
#    --init-data <copy> (no steady iterations, no ramp: the run carries straight on at DT) and
#    --sim-time SIM - t_ckpt, writing to out/${TAG}_cont so nothing of the first part is touched.
# The two time series are joined afterwards (part 2 shifted by t_ckpt) and averaged over the
# second half of the full run, as for every other run.
#
#SBATCH --job-name=tims-cont
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
: "${CASE:?}" "${SCALE:?}" "${VEL:?}" "${TAG:?}" "${DT:?}" "${SIM:?}" "${FIRST:?set FIRST=<first job id>}"

RUN=${CASE}_s${SCALE}_${TURB:-sst-gamma}
SRC=out/${TAG}/${RUN}
DST=out/${TAG}_cont/${RUN}
U=U${VEL}

if [[ -f "${SRC}/${U}_summary.json" ]]; then
  echo "first part finished (${SRC}/${U}_summary.json): nothing to continue"
  exit 0
fi
[[ -f "${SRC}/${U}_ckpt.dat.h5" ]] || { echo "no checkpoint in ${SRC}"; exit 1; }

# The driver writes a checkpoint after the report at which t - t_ckpt >= 2 s. Replay that rule
# over the reports in the first part's log, then keep the last checkpoint time whose report
# came before the checkpoint file was last written (a write cut off by the time limit would
# leave the previous checkpoint's time on the file).
LOG1=logs/attempt_${FIRST}.txt
CKPT_EPOCH=$(stat -c %Y "${SRC}/${U}_ckpt.dat.h5")
T_CKPT=$(python3 - "${LOG1}" "${CKPT_EPOCH}" <<'EOF'
import re, sys, datetime
log, ckpt = sys.argv[1], float(sys.argv[2])
last, t_ck = 0.0, 0.0
for line in open(log, errors="replace"):
    m = re.match(r"(\S+ \S+),\d+ INFO t=\s*([0-9.]+) s", line)
    if not m:
        continue
    t = float(m.group(2))
    when = datetime.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S").timestamp()
    if t - last >= 2.0 - 1e-9:
        last = t
        if when <= ckpt:
            t_ck = t
print(f"{t_ck:.4f}")
EOF
)
[[ "${T_CKPT}" != "0.0000" ]] || { echo "no checkpoint time found in ${LOG1}"; exit 1; }
REMAIN=$(python3 -c "print(f'{${SIM} - ${T_CKPT}:.4f}')")
echo "continuing ${RUN} ${U} from checkpoint at t = ${T_CKPT} s for ${REMAIN} s (job ${FIRST})"
[[ -n "${DRYRUN:-}" ]] && exit 0                 # DRYRUN=1 bash tims_continue.slurm: check only

mkdir -p "${DST}"
cp -p "${SRC}/${U}_ckpt.dat.h5" "${DST}/${U}_start.dat.h5"
printf '{"first_job": "%s", "t_ckpt": %s, "sim_time_total": %s, "remaining": %s, "start_data": "%s"}\n' \
  "${FIRST}" "${T_CKPT}" "${SIM}" "${REMAIN}" "${SRC}/${U}_ckpt.dat.h5" > "${DST}/${U}_continuation.json"

export STAGE=transient DRIVER=run_tims_case_v3.py TAG="${TAG}_cont" INIT="${DST}/${U}_start.dat.h5"
# EXTRA_ARGS carries the first part's other driver options (e.g. "--inner 20",
# "--solid-zones wall_body"), which the continuation must repeat.
export ARGS="--dt ${DT} --sim-time ${REMAIN} --ramp-factor 20 --ramp-steps 100 ${EXTRA_ARGS:-}"
exec bash "${SLURM_SUBMIT_DIR}/tims_job.slurm"
