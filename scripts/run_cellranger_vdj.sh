#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 --reference /path/to/custom_vdj_ref --manifest config/cellranger_libraries.csv --outdir results/cellranger [--execute] [--cores 8] [--mem 32]"; }
REF='' ; MANIFEST='' ; OUTDIR='results/cellranger' ; CORES=8 ; MEM=32 ; EXEC=0
while (($#)); do case "$1" in
 --reference) REF=$2;shift 2;; --manifest) MANIFEST=$2;shift 2;; --outdir) OUTDIR=$2;shift 2;; --cores) CORES=$2;shift 2;; --mem) MEM=$2;shift 2;; --execute) EXEC=1;shift;; *) usage;exit 2;; esac;done
[[ -n "$REF" && -n "$MANIFEST" ]] || { usage;exit 2; }
if ((EXEC));then command -v cellranger >/dev/null || { echo 'cellranger not found' >&2;exit 1; }; [[ -d "$REF" ]] || { echo 'Reference directory not found' >&2;exit 1; };fi
mkdir -p "$OUTDIR"; OUTDIR=$(cd "$OUTDIR" && pwd)
while IFS=, read -r run fastqs sample_prefix; do
 [[ "$run" == 'run' || -z "$run" ]] && continue
 [[ "$run" =~ ^[a-zA-Z0-9_-]+$ ]] || { echo "Invalid run ID: $run" >&2;exit 2; }
 cmd=(cellranger vdj --id="$run" --reference="$REF" --fastqs="$fastqs" --sample="$sample_prefix" --chain=IG --localcores="$CORES" --localmem="$MEM")
 printf 'cd %q && ' "$OUTDIR"; printf '%q ' "${cmd[@]}"; printf '\n'
 if ((EXEC));then [[ -d "$fastqs" ]] || { echo "FASTQ directory missing: $fastqs" >&2;exit 1; };(cd "$OUTDIR" && "${cmd[@]}");fi
done < "$MANIFEST"
