# JetClass-II dataset generation

This repository contains the event-generation and Delphes-production scripts used for the [JetClass-II dataset](https://huggingface.co/datasets/jet-universe/jetclass2).

## Repository layout

```text
.
├── run.sh                  # End-to-end MadGraph/Pythia/Delphes production
├── gen_configs/            # Process-specific generation configurations
├── delphes_cards/          # Delphes detector cards
└── delphes_analyzers/      # ROOT macro to convert Delphes ROOT files to ntuples
```

## Software needed on lxplus alma9

The current scripts expect:

- `MadGraph5_aMC@NLO` with `pythia8`, `hepmc`, and `lhapdf6` installed inside `HEPTools`
- `LHAPDF6` PDF data available to MadGraph/Pythia
- `Delphes`
- `ROOT` for the ntuple-making step

The examples below use EL9-compatible software on lxplus.

### 1. Load a ROOT/compiler environment

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
```

### 2. Install MadGraph + internal Pythia8/HepMC/LHAPDF

```bash
cd /path/where/you/want/the/tools, e.g. /afs/cern.ch/work/m/mpresill/private/tools
wget https://launchpadlibrarian.net/828125194/MG5_aMC_v3.5.12.tar.gz
tar -xzf MG5_aMC_v3.5.12.tar.gz
cd MG5_aMC_v3_5_12

cat > install_mg5.dat <<'MG5'
install pythia8
install hepmc
install lhapdf6
install model 2HDM
exit
MG5

./bin/mg5_aMC install_mg5.dat
```

This is installing a series of packaged within Madgraph software, which will be placed (for my examples) in `/afs/cern.ch/work/m/mpresill/private/tools/MG5_aMC_v3_5_12/HEPTools/`

#### Required patch for MadSpin

MadGraph's standalone f2py Makefile templates call f2py with `--include-paths=<dir>`
(single token). The f2py CLI shipped with numpy only recognizes the two-token form
`--include-paths <dir>`; with the `=` form the path is silently dropped, and MadSpin's
standalone spin-correlated decay matrix elements (needed whenever `madspin=ON` with
`spinmode onshell`, e.g. the `jetclass2/train_zz` config) fail to compile with
`error: unknown file type '' (from '--include-paths=...')`. Apply the fix once, right
after installing:

```bash
./fix_mg5_madspin_f2py.sh /afs/cern.ch/work/m/mpresill/private/tools/MG5_aMC_v3_5_12
```

(from the repository root; the script is idempotent and safe to re-run, e.g. after
re-downloading MG5).

Notes (skip if not expert):

- The existing JetClass-II resonance configs use the `2HDM` model.
- The provided `pp -> ZZ` prototype added in this repository uses the Standard Model (`import model sm`) and does not need any extra UFO model.
- If you prefer a site-wide LHAPDF installation, point `LHAPDFCONFIG` and `LHAPDF_DATA_PATH` in `run.sh` to that installation instead of the internal MG5 copy.

### 3. Install Delphes

```bash
cd /path/where/you/want/the/tools , e.g. /afs/cern.ch/work/m/mpresill/private/tools 
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz
tar -xzf Delphes-3.5.0.tar.gz
cd Delphes-3.5.0
make -j"$(nproc)"
ln -sf /eos/cms/store/group/upgrade/delphes/PhaseII/MinBias_100k.pileup MinBias_100k.pileup
```

The production script expects `MinBias_100k.pileup` to be available directly under `${DELPHES_PATH}`, so if you don't get errors from the last command it's just fine.

## Configure `run.sh`

From the repository root, edit the paths near the top of `run.sh`.
This is the example frmo my installation path:

```bash
G5_PATH=/afs/cern.ch/work/m/mpresill/private/tools/MG5_aMC_v3_5_12
DELPHES_PATH=/afs/cern.ch/work/m/mpresill/private/tools/Delphes-3.5.0
OUTPUT_PATH=/afs/cern.ch/work/m/mpresill/private/pilot-polarized-pnet/jet_class_test/jetclass2_generation/output
LHAPDFCONFIG=/afs/cern.ch/work/m/mpresill/private/tools/MG5_aMC_v3_5_12/HEPTools/lhapdf6_py3/bin/lhapdf-config
LHAPDF_DATA_PATH=/afs/cern.ch/work/m/mpresill/private/tools/MG5_aMC_v3_5_12/HEPTools/lhapdf6_py3/share/LHAPDF
PYTHIA8DATA=/afs/cern.ch/work/m/mpresill/private/tools/MG5_aMC_v3_5_12/HEPTools/pythia8/share/Pythia8/xmldoc```
````

Optional override:

```bash
export DELPHES_CARD_PATH=/absolute/path/to/custom_delphes_card.tcl
```

If `DELPHES_CARD_PATH` is not set, the script uses `delphes_cards/delphes_card_CMS_JetClassII_onlyFatJet.tcl`.

## Run MadGraph + Pythia + Delphes

From the repository root:

```bash
cd jetclass2_generation
./run.sh [process_name] [num_tot_events] [num_events_per_gen_step] [job_num]
```

Examples:

```bash
cd jetclass2_generation
./run.sh jetclass2/train_higgs2p 100 50 0
./run.sh jetclass2/train_zz 100 50 0
```

### Run 3 quick recipe (`jetclass2/train_zz`)

The `jetclass2/train_zz` config in this repository is already set for:

- 13.6 TeV proton-proton collisions (`ebeam1=6800`, `ebeam2=6800`)
- NNPDF3.1 LO (`lhaid=315000`)
- MadSpin Z decays (configured in `mg5_step2_madspin_card_templ.dat`)

Minimal commands:

```bash
cd jetclass2_generation
# edit run.sh once: MG5_PATH, DELPHES_PATH, OUTPUT_PATH, LHAPDFCONFIG, LHAPDF_DATA_PATH, PYTHIA8DATA
./run.sh jetclass2/train_zz 100000 100 0
```

Important runtime notes:

- `num_tot_events` must be divisible by `num_events_per_gen_step`.
- Output ROOT files are written to `${OUTPUT_PATH}/[process_name]/events_delphes_[job_num].root`.
- The current Higgs-like configs sample one random line from `mg5_params.dat` for each generation batch.
- `run.sh` creates a temporary work directory under `${OUTPUT_PATH}` and removes it after a successful job.

### Included process configurations

| Process | Command example | Notes |
| --- | --- | --- |
| `jetclass2/train_higgs2p` | `./run.sh jetclass2/train_higgs2p 100000 100 [job_num]` | Neutral resonance pair prototype |
| `jetclass2/train_higgspm2p` | `./run.sh jetclass2/train_higgspm2p 100000 100 [job_num]` | Charged resonance pair prototype |
| `jetclass2/train_higgs4p` | `./run.sh jetclass2/train_higgs4p 100000 100 [job_num]` | 4-parton resonance prototype |
| `jetclass2/train_qcd` | `./run.sh jetclass2/train_qcd 100000 100 [job_num]` | Pythia8-only QCD production |
| `jetclass2/train_zz` | `./run.sh jetclass2/train_zz 100000 100 [job_num]` | Standard Model `pp -> ZZ` prototype at 13.6 TeV with MadSpin Z decays |
| `jetclass2/train_zwpLhad` | `./run.sh jetclass2/train_zwLhad 100000 100 [job_num]` | Standard Model `pp -> Z W+L`, W+ longitudinal decaying hadronically and Z leptonically prototype at 13.6 TeV with MadSpin |

## How process configs are structured

A MadGraph+Pythia config directory may contain:

```text
mg5_step1.dat         # MadGraph process definition
mg5_step2_templ.dat   # MadGraph launch commands, with $NEVENT/$SEED placeholders
mg5_params.dat        # Optional parameter scan file; one line is sampled per batch
py8.dat               # Pythia8 card used by MG5aMC_PY8_interface
```

The default wrapper is `gen_configs/run_gen_default.sh`.

## New `pp -> ZZ` prototype

The prototype configuration is located in:

```text
gen_configs/jetclass2/train_zz
```

It uses:

- `import model sm`
- `generate p p > z z`
- `set ebeam1 6800` and `set ebeam2 6800` for 13.6 TeV proton-proton collisions
- `set lhaid 315000` (NNPDF3.1 LO via LHAPDF)
- `madspin=ON` with a dedicated MadSpin card in `mg5_step2_madspin_card_templ.dat`
- Z decays configured in MadSpin with leptonic and hadronic channels (`z -> l+ l-` and `z -> q q~`) rather than in the Pythia card

This makes it a simple starting point for LHC diboson studies on top of the existing JetClass-II workflow.

## Produce ntuples from Delphes output

The Delphes analyzer lives in `delphes_analyzers`.

```bash
cd jetclass2_generation/delphes_analyzers
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:/cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.1pre09-9fe9c/x86_64-el9-gcc13-opt/include

root -b -q 'makeNtuples.C++("events_delphes_higgs2p_example.root", "out.root", "JetPUPPIAK8", "GenJetAK8", true)'
```

Macro signature:

```cpp
void makeNtuples(TString inputFile,
                 TString outputFile,
                 TString jetBranch = "JetPUPPIAK8",
                 TString genjetBranch = "GenJetAK8",
                 bool assignQCDLabel = false,
                 bool debug = false)
```
