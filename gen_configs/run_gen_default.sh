#!/bin/bash -x

# a wrapper to generate hepmc file
# note: this can only generate events from a certain paramter choice

NEVENT=$1

# the MG process dir
MDIR=proc

# cd into current genpack's dir
cd "$( dirname "${BASH_SOURCE[0]}" )"

# sample one line from the paramter file if mg5_params.dat exists
if [ -f mg5_params.dat ]; then
    RANDOM_PARAMS=$(shuf -n 1 mg5_params.dat)
    IFS=' ' read -ra PARAMS <<< $RANDOM_PARAMS
fi

# step1: setup the MG dir for the first time
if [ ! -d $MDIR ]; then
    echo "Setup MG dir"
    $MG5_PATH/bin/mg5_aMC mg5_step1.dat
fi

# step2: generate event

## write mg5_step2.dat
cp -f mg5_step2_templ.dat mg5_step2.dat
## if mg5_params.dat exists, replace $1, $2, ... by the customized params
if [ -f mg5_params.dat ]; then
    for i in `seq 1 ${#PARAMS[@]}`; do
        sed -i "s/\$${i}/${PARAMS[$((i-1))]}/g" mg5_step2.dat
    done
fi

sed -i "s/\$NEVENT/$NEVENT/g" mg5_step2.dat # initialize nevent
sed -i "s/\$SEED/$RANDOM/g" mg5_step2.dat # initialize seed, important!

# if mg5_step2_run_card_templ exists, copy it to the MG dir
if [ -f mg5_step2_run_card_templ.dat ]; then
    cp -f mg5_step2_run_card_templ.dat $MDIR/Cards/run_card.dat
fi

## generate MG events
rm -rf $MDIR/Events/*
cat mg5_step2.dat | $MDIR/bin/generate_events pilotrun
if [ -f $MDIR/Events/pilotrun_decayed_1/unweighted_events.lhe.gz ]; then
    echo "Run MG with MadSpin. Use the MadSpin generated LHE"
    mv -f $MDIR/Events/pilotrun_decayed_1/unweighted_events.lhe.gz .
else
    mv -f $MDIR/Events/pilotrun/unweighted_events.lhe.gz .
fi

# run pythia
rm -f events.hepmc
LD_LIBRARY_PATH=$MG5_PATH/HEPTools/lib:$LD_LIBRARY_PATH $MG5_PATH/HEPTools/MG5aMC_PY8_interface/MG5aMC_PY8_interface py8.dat
