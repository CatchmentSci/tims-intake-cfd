#!/bin/bash
# Mesh screen: (re)mesh if asked, 100 steady iterations + 1 step at 0.3 m/s, then locate fast
# and tiny cells in every zone. One job per case and scale.
#   sbatch --ntasks=16 --export=ALL,CASE=phillips2000_full,SCALE=1.00,REMESH=1,EXTRA="--solid-zones wall_body" screen.slurm
#SBATCH --job-name=tims-screen
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=96G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "${SLURM_SUBMIT_DIR:-.}/site.sh"     # site settings: modules, Python environment, Fluent paths (site.example.sh)
export I_MPI_FABRICS=shm FLUENT_AFFINITY=0 I_MPI_PIN=off UCX_TLS=self,sm
export TIMS_FLUENT_ARGS=""
SCALE=${SCALE:-1.00}
GEOVER=${GEOVER:-}                     # geometry version subfolder, e.g. v5 ("" = current)
GEO=../geometry/out/${GEOVER:+${GEOVER}/}scale_${SCALE}/${CASE}
MESHOUT=out${GEOVER:+/${GEOVER}}
MESH=${MESHOUT}/${CASE}_s${SCALE}.msh.h5
TEST=out/screen/${GEOVER:+${GEOVER}_}${CASE}_s${SCALE}
N=${SLURM_NTASKS}
echo "screen case=${CASE} scale=${SCALE} remesh=${REMESH:-0} cores=${N} host=$(hostname)"
rm -rf "${TEST}"
if [[ "${REMESH:-0}" == 1 ]]; then
  python run_tims_case_v3.py mesh --geometry "${GEO}.msh" --cores "${N}" --out "${MESHOUT}"
fi
COMMON=(--mesh "${MESH}" --zones "${GEO}_zones.json" --cores "${N}" --out "${TEST}"
        --steady-init 100 --ramp-steps 0 --inner 2 --dt 0.02 ${EXTRA:-})
case "${CASE}" in
  *_external) python run_tims_case_v3.py external --unsteady "${COMMON[@]}" --velocities 0.3 \
                --r-values 0.2 --sim-time 0.02 ;;
  *_internal) python run_tims_case_v3.py internal --unsteady "${COMMON[@]}" \
                --tube-velocities 0.06 --point-steps 1 ;;
  *)          python run_tims_case_v3.py transient "${COMMON[@]}" --velocities 0.3 \
                --sim-time 0.02 --report-every 0.02 ;;
esac
FINAL=$(ls -t "${TEST}"/*/*_final.cas.h5 | head -1)
python locate_fast_cells.py "${FINAL}" --fast 0.9
