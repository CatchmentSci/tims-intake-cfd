# Site settings sourced by the Slurm templates (tims_job.slurm, tims_continue.slurm, screen.slurm,
# precursor.slurm). Copy to site.sh in this folder and edit for your cluster; site.sh is git-ignored.
# The values below are placeholders: module names, paths and versions differ between sites.

# Ansys Fluent 2025 R1 and Python 3.12
module purge
module load <your-ansys-fluent-2025R1-module>
module load <your-python-3.12-module>

# Python environment with the packages in requirements.txt (ansys-fluent-core etc.)
source <path-to-venv>/bin/activate

# PyFluent finds Fluent through AWP_ROOT<version>. Some site modules set only FLUENT_DIR, so set
# both explicitly. <ansys-install>/v251 is the folder that holds fluent/fluent25.1.0.
export AWP_ROOT251=<ansys-install>/v251
export PYFLUENT_FLUENT_ROOT="${FLUENT_DIR:-${AWP_ROOT251}/fluent}"

# RHEL 9 nodes may lack libnsl.so.1, which the Fluent host process needs. Ansys ships a copy;
# uncomment if Fluent fails to start with "libnsl.so.1: cannot open shared object file".
# export LD_PRELOAD="${AWP_ROOT251}/fluent/fluent25.1.0/addons/gpuapp/cromwell/amd/dependencies/lnamd64/libnsl.so.1${LD_PRELOAD:+:$LD_PRELOAD}"

# Optional: nodes the licence-retry resubmission should avoid (comma-separated), e.g. a known bad node
# export TIMS_EXCLUDE=
