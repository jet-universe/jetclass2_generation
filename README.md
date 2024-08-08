# JetClass-II dataset generation

This repo provides the generation scripts for the [JetClass-II dataset](https://huggingface.co/datasets/jet-universe/jetclass2).

## Setup environment

Please ensure the following software is configured:
- [MadGraph5_aMC@NLO](https://launchpad.net/mg5amcnlo) and [Pythia8](https://pythia.org): Download MadGraph, and install Pythia within it. Download the [`2HDM`](https://feynrules.irmp.ucl.ac.be/wiki/2HDM) model within MadGraph.
- [LHAPDF](https://lhapdf.hepforge.org): Install LHAPDF or configure an existing installation so that it can be interfaced with MadGraph.
- [Delphes](http://cp3.irmp.ucl.ac.be/projects/delphes): Install Delphes or configure an existing installation. Link the pileup file `/eos/cms/store/group/upgrade/delphes/PhaseII/MinBias_100k.pileup` within this directory.

The script `run.sh` produces events via MadGraph + Pythia + Delphes in one go. Please configure the environment variables in the script.

## Run MadGraph + Pythia + Delphes in one go
 
```bash
# Start the production by executing:
# ./run.sh [process_name] [num_tot_events] [num_events_per_gen_step] [job_num]

# For example:
./run.sh jetclass2/train_higgs2p 100 50 0
```

<details>
  <summary>Expand to see full setup for JetClass-II production.</summary>

| Process | Command |
| --- | --- |
| **`Res2P`** (neutral $X$) | `./run.sh jetclass2/train_higgs2p 100000 100 [job_num]`, `[job_num]` ranges from 0–191 |
| **`Res2P`** (charged $X$) | `./run.sh jetclass2/train_higgspm2p 100000 100 [job_num]`, `[job_num]` ranges from 0–143 |
| **`Res34P`** | `./run.sh jetclass2/train_higgs4p 100000 100 [job_num]`, `[job_num]` ranges from 0–1439 |
| **`QCD`:** | `./run.sh jetclass2/train_qcd 100000 100 [job_num]`, `[job_num]` ranges from 0–239 |

</details>

## MadGraph + Pythia step

The MadGraph + Pythia job is configured by `[process_name]`, which takes a subdirectory path in [`gen_configs`](gen_configs). The following files are present to define the job:
```yaml
├── mg5_step1.dat         # Definition of the MadGraph process
├── mg5_step2_templ.dat   # Configuration for launching a MadGraph process. $1, $2, etc., serve as placeholders for parameters from a given list
├── mg5_params.dat        # Parameter list to replace $1, $2, etc.
└── py8.dat               # Pythia8 data file
```

## Delphes step

The Delphes step is configured by a Delphes card. Three Delphes cards are provided in `delphes_cards`, which are derived from the CMS detector conditions in Run 2/3, adding track smearing as done in the JetClass-I production [[JetClass Delphes card]](https://github.com/jet-universe/jetclass_generation/blob/main/delphes_card.tcl), and emulating pileup with $<\mu>=$50, mitigating the PU with the PUPPI algorithm.

The three JetClass-II Delphes cards have gradually increasing content and computational costs:

| Delphes card | Description |
| --- | --- |
| **[`delphes_card_CMS_JetClassII_onlyFatJet.tcl`](delphes_cards/delphes_card_CMS_JetClassII_onlyFatJet.tcl)** | This data card only produces large-*R* jets `JetPUPPIAK8`, `JetPUPPIAK15` and the corresponding GEN-jets `GenJetAK8`, `GenJetAK15`. |
| **[`delphes_card_CMS_JetClassII_lite.tcl`](delphes_cards/delphes_card_CMS_JetClassII_lite.tcl)** | This data card also perserves small-*R* jets `JetPUPPI` (with *R*=0.4) and other event-level objects, including `Electron`, `Muon`, `Photon`, `PuppiMissingET`, etc. |
| **[`delphes_card_CMS_JetClassII.tcl`](delphes_cards/delphes_card_CMS_JetClassII.tcl)** | This data card includes consistent information as in the nominal CMS card, also including the CHS jets with different radii: `Jet`, `JetAK8`, and `JetAK15`. |

## Produce ntuples from Delphes

Enter `delphes_analyzers`, compile, and run the `makeNtuples.C` macro to convert a Delphes file into ntuples. The ntuple is jet-based, with its branch definition available [here](https://github.com/jet-universe/sophon?tab=readme-ov-file#variable-details).

The following example works on EL9 machines.

```bash
# Setup environment

source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:/cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.1pre09-9fe9c/x86_64-el9-gcc13-opt/include

# Compile and run macro
root -b -q 'makeNtuples.C++("events_delphes_higgs2p_example.root", "out.root", "JetPUPPIAK8", "GenJetAK8", true)'

## defination: 
##  void makeNtuples(TString inputFile, TString outputFile, TString jetBranch = "JetPUPPIAK8", TString genjetBranch = "GenJetAK8", bool assignQCDLabel = false, bool debug = false)
```
