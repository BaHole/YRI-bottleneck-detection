#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 code/compare_demography.py --data data --output results
