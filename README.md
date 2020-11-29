# XGBClassifier

XGBClassifier： a tool for predicting single amino acid variant pathogenicity.

## Example
```bash
python3 main.py example example/O96017.pdb O96017 example/O96017_sav.txt MAESTRO_OSX_x64  p2rank_2.2 CHEK2.fasta NP_009125.1 example/O96017_fathmm.txt
```

## Installation
Download
```bash
main.py
evolutionaryAnalysis.jl
XGBoostClassifier.joblib
```
Required Input Files:

```bash
SAV file
PDB file
Orthologs file in FASTA format
FATHMM output file
```
Required Packages

Python Packages

```bash
Bio
sklearn
xgboost
prody
rhapsody
skbio
pyrosetta
```
Other Packages
```bash
foldx(rotabase.txt)
freesasa: command line tool
maestro: stand-alone software
P2Rank: stand-alone software
plmc-master: stand-alone software
clustal-omega: command line tool
paup: command line tool
```
## Usage
Args
```bash
[1] working Dir: relative path
[2] PDF file: full path
[3] uniprot accession
[4] SAV file (All SAVs should be from the SAME CHAIN; make sure all the SAVs are in the PDB)
[5] MAESTRO directory
[6] P2Rank directory
[7] Orthologs file (including the human amino acid)
[8] Protein accession
[9] FATHMM output file
```
