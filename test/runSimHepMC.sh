if [[ -z $(which parallel) ]] ; then 
  echo "parallel is not available, install it first"
  exit 1
fi

NJOBS=20
NCORES=20
NEVENTSPERJOB=100
OUTDIR=real_config/PYTHIA_PbPb_MB
LAYOUT=na6pLayout_real.ini
GENERATOR=genHEPMCWithSpectTracking.C+
SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))
EVENTSTOSKIP=$(($SEEDSTART))

parallel -j $NCORES \
  'na6psim -n '"$NEVENTSPERJOB"' \
     -g '"$NA6PROOT_ROOT"'/share/test/'"$GENERATOR"'\(\"real_config/PYTHIA_PbPb_MB/{1}/GenPYTHIA_HepMC3.root\",false,0\) \
     --load-ini '"$LAYOUT"' \
     --disable-write-ini \
     --configKeyValues "keyval.output_dir='"$OUTDIR"'/{1}" \
     > '"$OUTDIR"'/{1}/log_sim_seed{1}.txt' \
  ::: $(seq $SEEDSTART $SEEDEND)