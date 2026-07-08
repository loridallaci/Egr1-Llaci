#!/bin/bash
# ---------------------------------------------------------------------
# Deploy + submit Suppl Fig 2 panels a/b (cortex dev cell numbers + UMAPs)
# to HTCF. Copies the R script and the sbatch to a remote work dir and
# submits the job. Requires a LIVE ssh master connection to HTCF (HTCF
# needs interactive password+Duo, which the shared config multiplexes):
#
#   ssh -fN login.htcf.wustl.edu      # open once, approve Duo; persists 8h
#   bash deploy_SupplFig2_HTCF.sh     # then this runs with no re-auth
#
# Verify the master is up first:  ssh -O check login.htcf.wustl.edu
# ---------------------------------------------------------------------
set -euo pipefail

HOST="login.htcf.wustl.edu"
REMOTE_DIR="/home/lllaci/cortex_suppfig2"
HERE="$(cd "$(dirname "$0")" && pwd)"

echo ">> checking ssh master connection to $HOST ..."
if ! ssh -O check "$HOST" 2>/dev/null; then
  echo "!! no live master connection. Open one first (interactive, Duo):" >&2
  echo "     ssh -fN $HOST" >&2
  exit 1
fi

echo ">> creating remote dir $REMOTE_DIR"
ssh "$HOST" "mkdir -p '$REMOTE_DIR'"

echo ">> copying script + sbatch"
scp "$HERE/Figure_SupplFig2_cellnumbers_umaps.R" \
    "$HERE/run_SupplFig2_cellnumbers_umaps_HTCF.sbatch" \
    "$HOST:$REMOTE_DIR/"

echo ">> submitting job (all 5 stages)"
ssh "$HOST" "cd '$REMOTE_DIR' && sbatch run_SupplFig2_cellnumbers_umaps_HTCF.sbatch"

echo ">> submitted. Monitor with:  ssh $HOST 'squeue -u l.llaci'"
echo ">> outputs will land in:      $REMOTE_DIR/output/"
echo ">> pull them back with:       scp '$HOST:$REMOTE_DIR/output/SupplFig2_*' '$HERE/output/'"
