NJOBS=20
NCORES=20
NEVENTSPERJOB=100
ENERGY=40
MINRAP=0.5
MAXRAP=4.5
CONFIG=PbPb_MB.cmnd
OUTDIR=data/PYTHIA_PbPb_MB

SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))

parallel -j $NCORES \
  'mkdir -p '"$OUTDIR"'/{1} && \
   ./GenPYTHIA --events '"$NEVENTSPERJOB"' --energy '"$ENERGY"' --ymin '"$MINRAP"' --ymax '"$MAXRAP"' --cmnd '"$CONFIG"' --out '"$OUTDIR"'/{1} > '"$OUTDIR"'/{1}/log_sim_seed{1}.txt' \
  ::: $(seq $SEEDSTART $SEEDEND)