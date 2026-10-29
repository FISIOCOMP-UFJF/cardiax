#!/bin/bash
cd "$(dirname "$0")"
BIN=$(pwd)/../../build/app/cardiac_cycle
MESH=$(pwd)/malha_descarregada_ldrb_ms.xml
COMMON="-num_threads 1 -ksp_type preonly -pc_type lu -c ToRORdLand -dt 0.02 -dt_mech 1 -save_state 5"
export OMP_NUM_THREADS=1

mkdir -p regazzoni/full regazzoni/rest
(cd regazzoni/full && $BIN -f $MESH $COMMON -circ 1 -rr 20 -beats 1 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_circulation.dat regazzoni/full/
(cd regazzoni/rest && $BIN -f $MESH $COMMON -circ 1 -rr 20 -beats 1 -restore ../full/checkpoint_t_10.00ms.h5 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_circulation.dat regazzoni/rest/
h5diff regazzoni/full/checkpoint_t_15.00ms.h5 regazzoni/rest/checkpoint_t_15.00ms.h5
h5diff regazzoni/full/checkpoint_t_20.00ms.h5 regazzoni/rest/checkpoint_t_20.00ms.h5
tail -n 10 regazzoni/full/malha_descarregada_ldrb_ms.xml_circulation.dat > regazzoni/full/malha_descarregada_ldrb_ms.xml_circulation.dat.tail
tail -n 10 regazzoni/rest/malha_descarregada_ldrb_ms.xml_circulation.dat > regazzoni/rest/malha_descarregada_ldrb_ms.xml_circulation.dat.tail
diff regazzoni/full/malha_descarregada_ldrb_ms.xml_circulation.dat.tail regazzoni/rest/malha_descarregada_ldrb_ms.xml_circulation.dat.tail

mkdir -p regazzoni_mef/full regazzoni_mef/rest
(cd regazzoni_mef/full && $BIN -f $MESH $COMMON -circ 1 -rr 20 -beats 1 -lam 1 -lamrate 1 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_circulation.dat regazzoni_mef/full/
(cd regazzoni_mef/rest && $BIN -f $MESH $COMMON -circ 1 -rr 20 -beats 1 -lam 1 -lamrate 1 -restore ../full/checkpoint_t_10.00ms.h5 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_circulation.dat regazzoni_mef/rest/
h5diff regazzoni_mef/full/checkpoint_t_15.00ms.h5 regazzoni_mef/rest/checkpoint_t_15.00ms.h5
h5diff regazzoni_mef/full/checkpoint_t_20.00ms.h5 regazzoni_mef/rest/checkpoint_t_20.00ms.h5
tail -n 10 regazzoni_mef/full/malha_descarregada_ldrb_ms.xml_circulation.dat > regazzoni_mef/full/malha_descarregada_ldrb_ms.xml_circulation.dat.tail
tail -n 10 regazzoni_mef/rest/malha_descarregada_ldrb_ms.xml_circulation.dat > regazzoni_mef/rest/malha_descarregada_ldrb_ms.xml_circulation.dat.tail
diff regazzoni_mef/full/malha_descarregada_ldrb_ms.xml_circulation.dat.tail regazzoni_mef/rest/malha_descarregada_ldrb_ms.xml_circulation.dat.tail

mkdir -p piecewise/full piecewise/rest
(cd piecewise/full && $BIN -f $MESH $COMMON -pw 1 -pwinit 40 -pwedp 10 -rr 20 -beats 1 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_piecewise.dat piecewise/full/
(cd piecewise/rest && $BIN -f $MESH $COMMON -pw 1 -pwinit 40 -pwedp 10 -rr 20 -beats 1 -restore ../full/checkpoint_t_35.00ms.h5 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_piecewise.dat piecewise/rest/
for t in 40 45 50 55 60; do h5diff piecewise/full/checkpoint_t_$t.00ms.h5 piecewise/rest/checkpoint_t_$t.00ms.h5; done
tail -n 25 piecewise/full/malha_descarregada_ldrb_ms.xml_piecewise.dat > piecewise/full/malha_descarregada_ldrb_ms.xml_piecewise.dat.tail
tail -n 25 piecewise/rest/malha_descarregada_ldrb_ms.xml_piecewise.dat > piecewise/rest/malha_descarregada_ldrb_ms.xml_piecewise.dat.tail
diff piecewise/full/malha_descarregada_ldrb_ms.xml_piecewise.dat.tail piecewise/rest/malha_descarregada_ldrb_ms.xml_piecewise.dat.tail

mkdir -p legacy/full legacy/rest
(cd legacy/full && $BIN -f $MESH $COMMON -t 20 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_pvloop.txt legacy/full/
(cd legacy/rest && $BIN -f $MESH $COMMON -t 20 -restore ../full/checkpoint_t_10.00ms.h5 > log.txt 2>&1)
cp malha_descarregada_ldrb_ms.xml_pvloop.txt legacy/rest/
h5diff legacy/full/checkpoint_t_15.00ms.h5 legacy/rest/checkpoint_t_15.00ms.h5
h5diff legacy/full/checkpoint_t_20.00ms.h5 legacy/rest/checkpoint_t_20.00ms.h5
tail -n 10 legacy/full/malha_descarregada_ldrb_ms.xml_pvloop.txt > legacy/full/malha_descarregada_ldrb_ms.xml_pvloop.txt.tail
tail -n 10 legacy/rest/malha_descarregada_ldrb_ms.xml_pvloop.txt > legacy/rest/malha_descarregada_ldrb_ms.xml_pvloop.txt.tail
diff legacy/full/malha_descarregada_ldrb_ms.xml_pvloop.txt.tail legacy/rest/malha_descarregada_ldrb_ms.xml_pvloop.txt.tail
