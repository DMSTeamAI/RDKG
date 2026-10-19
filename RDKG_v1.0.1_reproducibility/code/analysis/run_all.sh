#!/bin/bash
# Usage: cd code/analysis ; bash run_all.sh   (run from this folder; outputs are written here)
# Needs: python3 with pandas numpy scipy. Downloads the v1.0.0 release files and hp.obo into ./data
set -e
mkdir -p data && cd data
[ -f edges_no_trials.csv ] || curl -L -o edges_no_trials.csv https://github.com/wangjl99/RDKG/releases/download/v1.0.0/edges_no_trials.csv
[ -f nodes_no_trials.csv ] || curl -L -o nodes_no_trials.csv https://github.com/wangjl99/RDKG/releases/download/v1.0.0/nodes_no_trials.csv
[ -f hp.obo ] || curl -L -o hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/hp.obo
cd ..
export RDKG_DIR=$PWD/data
python3 a1.py; python3 a2.py; python3 a3.py; python3 a4.py; python3 a6.py; python3 a7.py; python3 a8.py
python3 a9.py; python3 a10.py eval; python3 a10.py MONDO:0009265   # DrugBank-only repurposing (replaces Marfan example)
python3 ../release/clean_release.py   # writes ./release_v1.0.1/
