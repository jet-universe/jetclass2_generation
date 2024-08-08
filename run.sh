#!/bin/bash -x

PROC=$1
NEVENT=$2
NEVENT_GEN=$3
JOBNUM=$4

# Setup environment

# ============ basic configuration ============
MG5_PATH=/PATH/TO/MG5_aMC_v2_9_18
DELPHES_PATH=/PATH/TO/Delphes-3.5.0
OUTPUT_PATH=/PATH/TO/OUTPUT_DIR

## some env variables are required by the softwares
LHAPDFCONFIG=/PATH/TO/lhapdf-config
LHAPDF_DATA_PATH=/PATH/TO/share/LHAPDF
PYTHIA8DATA=/PATH/TO/MG5_aMC_v2_9_18/HEPTools/pythia8/share/Pythia8/xmldoc

## fixed configuration
GENCFG_PATH=$(realpath gen_configs)
DELPHES_CARD_PATH=$(realpath delphes_cards)/delphes_card_CMS_JetClassII_onlyFatJet.tcl
OUTPUT_PATH=$(realpath $OUTPUT_PATH)

# =============================================

# Create workdir
RANDSTR=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 10; echo)
WORKDIR=$(realpath $OUTPUT_PATH)/workdir_$(date +%y%m%d-%H%M%S)_${RANDSTR}_$(echo "$PROC" | sed 's/\//_/g')_$JOBNUM
mkdir -p $WORKDIR

cd $WORKDIR

generate_delphes(){
    # all GEN production logic inside the gen_configs folder
    # generate GEN events
    rm -f events.hepmc
    MG5_PATH=$MG5_PATH ./run_gen.sh $NEVENT_GEN

    # run delphes
    ln -s $DELPHES_PATH/MinBias_100k.pileup .
    rm -f events_delphes.root
    $DELPHES_PATH/DelphesHepMC2 $DELPHES_CARD_PATH events_delphes.root events.hepmc
    return $?
}

# generate delphes, in a batch of NEVENT_GEN
nbatch=$((NEVENT / NEVENT_GEN))

for ((i=0; i<nbatch; i++)); do

    echo "Batch: $i"

    # copy genpack if not exist
    cd $WORKDIR
    if [ ! -d "proc_base" ]; then
        mkdir proc_base
        cp -r $GENCFG_PATH/$PROC/* proc_base/
        # if genpack does not have a run_gen.sh, use the default
        if [ ! -f "proc_base/run_gen.sh" ]; then
            cp $GENCFG_PATH/run_gen_default.sh proc_base/run_gen.sh
        fi
    fi
    cd $WORKDIR/proc_base

    generate_delphes

    # if return code is 0
    if [ $? -eq 0 ]; then
        # successful
        mv events_delphes.root $WORKDIR/events_delphes_$i.root
    fi
    cd $WORKDIR

    # intermediate file merging for every 100 batches
    if [ $(((i+1) % 100)) -eq 0 ]; then
        hadd -f $WORKDIR/merged_events_delphes_$i.root $WORKDIR/events_delphes_*.root
        if [ $? -eq 0 ]; then
            rm -f $WORKDIR/events_delphes_*.root
        fi
    fi

done

# combine all root
if [ $nbatch -eq 1 ]; then
    mv $WORKDIR/events_delphes_0.root events_delphes.root
else
    hadd -f events_delphes.root $WORKDIR/*.root
fi
mkdir -p $OUTPUT_PATH/$PROC

# transfer the file
mv -f events_delphes.root $OUTPUT_PATH/$PROC/events_delphes_$JOBNUM.root

# remove workspace
rm -rf $WORKDIR

echo -e "\033[1mJob done. Generated $NEVENT events for $PROC.\033[0m"
echo -e "\033[1mDelphes file path: $OUTPUT_PATH/$PROC/events_delphes_$JOBNUM.root\033[0m"
