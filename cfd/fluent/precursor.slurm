#!/bin/bash
# Flume inflow precursor (Goharrokhi flume): mesh the empty flume duct built by
#   python geometry/tims_geometry.py goharrokhi2019 --domain flume_goharrokhi --duct
# into out/v5/, then recycle its outlet profile until the flow is developed, at each velocity:
#   sbatch --export=ALL,DUCT=duct_gflume,VEL="0.15 0.30 0.60",TAG=goharrokhi_precursor precursor.slurm
# Writes out/${TAG}/U<U>_profile.csv for run_tims_case_v3.py transient --inlet-profile.
#SBATCH --job-name=tims-precursor
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "${SLURM_SUBMIT_DIR:-.}/site.sh"     # site settings: modules, Python environment, Fluent paths (site.example.sh)
export I_MPI_FABRICS=shm FLUENT_AFFINITY=0 I_MPI_PIN=off UCX_TLS=self,sm
export TIMS_FLUENT_ARGS=""
DUCT=${DUCT:-duct_gflume}
GEO=../geometry/out/v5/scale_1.00/${DUCT}
MESH=out/v5/${DUCT}_s1.00.msh.h5
N=${SLURM_NTASKS}
echo "precursor duct=${DUCT} velocities=${VEL} cores=${N} host=$(hostname)"
python run_tims_case_v3.py mesh --geometry "${GEO}.msh" --cores "${N}" --out out/v5
python run_tims_case_v3.py precursor --mesh "${MESH}" --zones "${GEO}_zones.json" --cores "${N}" \
  --velocities ${VEL} --out "out/${TAG:-goharrokhi_precursor}" ${ARGS:-}
ls -l "out/${TAG:-goharrokhi_precursor}"/*_profile.csv
