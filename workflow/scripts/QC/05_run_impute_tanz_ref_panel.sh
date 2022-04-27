#!/bin/bash
OUT_DIR='../Imputed/Tanz/'
parallel -j 23 ./05_impute_tanz_ref_panel.sh {} ">" $OUT_DIR"chr"{}".impute.log" ::: {{1..22},"X"}
