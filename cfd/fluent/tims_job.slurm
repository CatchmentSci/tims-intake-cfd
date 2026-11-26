#!/bin/bash
# Slurm job for one TIMS case stage. Site settings (modules, Python environment, Fluent paths)
# come from site.sh (copy site.example.sh); pass --account / --partition to sbatch as your site needs.
#
#   sbatch --export=ALL,STAGE=mesh,SCALE=1.00 tims_job.slurm
#   sbatch --export=ALL,STAGE=setup,SCALE=1.00 tims_job.slurm
#   sbatch --export=ALL,STAGE=sweep,SCALE=1.00 tims_job.slurm
#   sbatch --array=0-2 --export=ALL,STAGE=mesh tims_job.slurm   # all three GCI meshes
#
#SBATCH --job-name=tims
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

# --- software ----------------------------------------------------------------------------
source "${SLURM_SUBMIT_DIR:-.}/site.sh"     # site settings: modules, Python environment, Fluent paths (site.example.sh)

# Multi-node runs only: give Fluent a host file (on one node it would override -t).
# On one node, keep Intel MPI on shared memory: some nodes (e.g. hmem*) have no usable
# UCX/InfiniBand transport and the Fluent launch otherwise hangs.
TIMS_FLUENT_ARGS=""
if [[ "${SLURM_JOB_NUM_NODES:-1}" -eq 1 ]]; then
  export I_MPI_FABRICS=shm
  # Large shared nodes (here 256 cores) can hold several jobs. Two Fluent runs pinned to the same cores crawled at
  # ~6 min per iteration, so leave placement to the job's cgroup instead.
  export FLUENT_AFFINITY=0
  export I_MPI_PIN=off
  export UCX_TLS=self,sm
else
  scontrol show hostnames "${SLURM_JOB_NODELIST}" > "logs/hosts_${SLURM_JOB_ID}.txt"
  TIMS_FLUENT_ARGS="-cnf=logs/hosts_${SLURM_JOB_ID}.txt"
fi
export TIMS_FLUENT_ARGS

# --- case --------------------------------------------------------------------------------
SCALES=(1.65 1.00 0.64)                 # coarse, medium, fine (the grid-convergence set)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then SCALE=${SCALES[$SLURM_ARRAY_TASK_ID]}; fi
SCALE=${SCALE:-1.00}
CASE=${CASE:-phillips2000}
TURB=${TURB:-sst-gamma}
GEOVER=${GEOVER:-}                     # geometry version subfolder, e.g. v5 ("" = current)
GEO=../geometry/out/${GEOVER:+${GEOVER}/}scale_${SCALE}/${CASE}
MESH=out/${GEOVER:+${GEOVER}/}${CASE}_s${SCALE}.msh.h5
echo "stage=${STAGE:-} scale=${SCALE} case=${CASE} turb=${TURB} host=$(hostname) cores=${SLURM_NTASKS}"

run_stage() {
  # ARGS passes any further run_tims_case.py options, e.g. ARGS="--dt 0.02 --sim-time 20";
  # DRIVER runs another copy of the driver (to test a change while jobs use the current one)
  case "${STAGE:?set STAGE=mesh|setup|sweep|transient|external|internal}" in
    mesh)  python "${DRIVER:-run_tims_case_v3.py}" mesh  --geometry "${GEO}.msh" --cores "${SLURM_NTASKS}" ;;
    *)     python "${DRIVER:-run_tims_case_v3.py}" "${STAGE}" --mesh "${MESH}" --zones "${GEO}_zones.json" \
                  --turbulence "${TURB}" --cores "${SLURM_NTASKS}" ${VEL:+--velocities ${VEL}} \
                  ${TAG:+--out out/${TAG}} ${PSEUDO:+--pseudo ${PSEUDO}} \
                  ${TSF:+--time-scale-factor ${TSF}} ${MAXIT:+--max-iter ${MAXIT}} \
                  ${INIT:+--init-data ${INIT}} ${ARGS:-} ;;
  esac
}

# --- licences ----------------------------------------------------------------------------
# Parallel Fluent draws one anshpc licence per core beyond the 4 included with the solver, from
# a shared site-wide pool (often nearly all in use), and aborts with
# "Unable to spawn node: license not available" when the pool is short. The pool cannot be
# polled from the compute nodes (lmutil reaches lmgrd but not the vendor daemon), so the stage
# itself is the test: on a licence failure, wait 5 min and launch again, for up to LIC_WAIT_MIN
# minutes. After that the job resubmits itself (for clusters without Slurm licence tracking and
# JobRequeue=0), up to MAX_RETRIES times: by default ~7 days of trying.
ATTEMPT_LOG="logs/attempt_${SLURM_JOB_ID}.txt"

# --- shutdown watchdog -------------------------------------------------------------------
# Fluent 2025R1 sometimes segfaults while shutting down after the final write and then hangs
# (for up to 11 h, each time holding 32 cores and its dependants).
# Once the driver has logged its time average and Fluent its shutdown, and the log has not
# grown for HANG_MIN minutes, stop waiting: the job then ends normally and Slurm kills what
# is left of Fluent.
HUNG_FLAG="logs/hung_${SLURM_JOB_ID}.flag"
watch_shutdown() {
  local pid=$1 idle=0 size=-1 s
  while kill -0 "${pid}" 2>/dev/null; do
    sleep 60
    grep -q "time-averaged over" "${ATTEMPT_LOG}" || continue
    grep -qE "Shutdown Complete|Stack backtrace" "${ATTEMPT_LOG}" || continue
    s=$(stat -c %s "${ATTEMPT_LOG}")
    if (( s == size )); then idle=$(( idle + 1 )); else idle=0; size=${s}; fi
    if (( idle >= ${HANG_MIN:-10} )); then
      echo "$(date '+%F %T') Fluent still running ${idle} min after shutdown: ending the job"
      touch "${HUNG_FLAG}"
      kill -KILL "${pid}"
      return
    fi
  done
}

for (( waited=0; ; waited+=5 )); do
  set +e
  run_stage > >(tee "${ATTEMPT_LOG}") 2>&1 &
  stage_pid=$!
  watch_shutdown "${stage_pid}" &
  watch_pid=$!
  wait "${stage_pid}"
  rc=$?
  kill "${watch_pid}" 2>/dev/null
  sleep 5                                # let tee finish the log before it is searched
  [[ -f "${HUNG_FLAG}" ]] && rc=0
  set -e
  (( rc == 0 )) && exit 0
  grep -qi "license not available" "${ATTEMPT_LOG}" || exit "${rc}"
  echo "$(date '+%F %T') licence not available (waited ${waited} min, attempt ${RETRY:-0})"
  (( waited + 5 > ${LIC_WAIT_MIN:-240} )) && break
  sleep 300
done
if (( ${RETRY:-0} < ${MAX_RETRIES:-42} )); then
  # The case variables are resolved here so an array task resubmits as a plain job
  env -u SLURM_ARRAY_TASK_ID -u SLURM_MEM_PER_NODE \
    sbatch ${TIMS_EXCLUDE:+--exclude="${TIMS_EXCLUDE}"} --mail-type=END,FAIL --time="$(squeue -h -j "${SLURM_JOB_ID}" -o %l)" --job-name="${SLURM_JOB_NAME}" --ntasks="${SLURM_NTASKS}" --mem="${SLURM_MEM_PER_NODE}M" \
    --export=ALL,RETRY=$(( ${RETRY:-0} + 1 )),SCALE="${SCALE}",CASE="${CASE}",TURB="${TURB}",GEOVER="${GEOVER}" \
    "${SLURM_SUBMIT_DIR}/tims_job.slurm"
fi
exit 3
