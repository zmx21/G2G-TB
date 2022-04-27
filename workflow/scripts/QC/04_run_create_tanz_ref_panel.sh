#!/bin/bash
OUT_DIR="../../../WGS/Ref_Panel/"
parallel -j 23 ./04_create_tanz_ref_panel.sh {} ">" $OUT_DIR"chr"{}".create.tanz.log" ::: {{1..22},"X"}
