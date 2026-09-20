#!/usr/bin/env bash
# Optional reruns with separately installed, appropriately licensed upstream tools.
# Usage: bash code/run_inference.sh epos|blockbuster|fitcoal
set -euo pipefail
cd "$(dirname "$0")/.."
mode="${1:-}"
mkdir -p results/inference results/inference/normalized_sfs
python3 - <<'PY'
from pathlib import Path
from code.compare_demography import CASES, sfs_values, validate_sfs
folder = Path('data')
validate_sfs(folder)
for _, name, _, _, _, _, _ in CASES:
    values = sfs_values(folder / name)
    (Path('results/inference/normalized_sfs') / name).write_text(
        ' '.join(str(int(x)) if x.is_integer() else str(x) for x in values) + '\n')
PY
case "$mode" in
  epos)
    : "${EPOS_BIN:?Set EPOS_BIN to an installed EPOS executable}"
    [[ -x "$EPOS_BIN" ]] || { echo "EPOS_BIN is not executable" >&2; exit 2; }
    for record in 'Elife:Elife.YRI.4usfs:5334461' 'Science:FitCoal.YRI.usfs:826650000' 'MBE:MBE.YRI.usfs:1350734279'; do
      IFS=: read -r name file length <<< "$record"
      "$EPOS_BIN" "results/inference/normalized_sfs/$file" -u 1.25e-8 -U -l "$length" \
        > "results/inference/${name}.epos.log" 2>&1
    done
    ;;
  blockbuster)
    : "${BLOCKBUSTER_SH:?Set BLOCKBUSTER_SH to an installed blockbuster.sh}"
    [[ -f "$BLOCKBUSTER_SH" ]] || { echo "BLOCKBUSTER_SH does not exist" >&2; exit 2; }
    for record in 'Elife:Elife.YRI.4usfs:5334461:1' 'Science:FitCoal.YRI.usfs:826650000:1' \
                  'MBE:MBE.YRI.usfs:1350734279:1' 'Zhen:GR.YRI.4fold.sfs:5767108:0'; do
      IFS=: read -r name file length oriented <<< "$record"
      mkdir -p "results/inference/blockbuster_$name"
      bash "$BLOCKBUSTER_SH" --sfs "$(pwd)/results/inference/normalized_sfs/$file" \
        --prefixe_directory "$(pwd)/results/inference/blockbuster_$name" \
        -L "$length" -m 1.25e-8 -g 29 -o "$oriented" \
        > "results/inference/${name}.blockbuster.log" 2>&1
    done
    ;;
  fitcoal)
    : "${FITCOAL_JAR:?Set FITCOAL_JAR to FitCoal.jar}"
    : "${FITCOAL_TABLES:?Set FITCOAL_TABLES to the compatible tables directory}"
    [[ -f "$FITCOAL_JAR" && -d "$FITCOAL_TABLES" ]] || { echo "FitCoal inputs missing" >&2; exit 2; }
    java -cp "$FITCOAL_JAR" FitCoal.calculate.RextFitCoal_SinglePopDecoder \
      -table "$FITCOAL_TABLES" -input "$(pwd)/results/inference/normalized_sfs/FitCoal.YRI.usfs" \
      -output "$(pwd)/results/inference/Science.YRI" \
      -mutationRate 0.0000125 -generationTime 29 -genomeLength 826650 \
      > results/inference/Science.fitcoal.log 2>&1
    ;;
  *)
    echo 'Usage: bash code/run_inference.sh epos|blockbuster|fitcoal' >&2
    exit 2
    ;;
esac
