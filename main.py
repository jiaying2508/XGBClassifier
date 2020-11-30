#run XGBoost Classifier to predict SAV disease pathogenicity
'''
python3 main.py example example/O96017.pdb O96017 example/O96017_sav.txt MAESTRO_OSX_x64  p2rank_2.2 CHEK2.fasta NP_009125.1 example/O96017_fathmm.txt
foldx in the same dir as main.py
freesasa is command line tool

args:
    [1] working Dir: relative path
    [2] old pdb: full path
    [3] uniprot accession
    [4] change_list_file/ FROM THE SAME CHAIN ON PDB FILE / make sure all the SAVs are in the PDB
    [5] MAESTRO dir
    [6] P2Rank dir
    [7] orthologs sequence file
    [8] protein accession
    [9] fathmm output file
'''
import os
import re
import sys

import shutil
from shutil import copyfile
from skbio import TreeNode
from io import StringIO
from Bio import Phylo
from Bio import SeqIO
import numpy as np
from scipy import stats
import random
import math

import prody
from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.toolbox import cleanATOM
from pyrosetta.teaching import *
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta.rosetta.protocols.docking import *

import pandas as pd
from evcouplings.couplings import CouplingsModel
from evcouplings.mutate import predict_mutation_table, single_mutant_matrix
import rhapsody as rd

currDir = os.getcwd()
#print(currDir)

workingDir = sys.argv[1]
proteinAccession = sys.argv[8]

################################################################################
#codon dictionary
codon3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
codon1to3 = {}
for key in codon3to1.keys():
    codon1to3[codon3to1[key]] = key

################################################################################
'''
clean PDB file
'''
################################################################################
#os.system('')
savfile = sys.argv[4]
sav_list = []
sav_1codon_list = []
sav_dict = {}
sav_1codon_dict = {}
rawDict = {}
posDict = {}
evDict = {}
#process mutation
file = open(savfile)
wt = ''
pos = ''
var = ''
for line in file:
    line = line.rstrip()
    m = re.search('([A-Za-z]+)(\d+)([A-Za-z]+)', line)
    wt = m.group(1).upper()
    if len(wt) == 1:
        wt = codon1to3[wt]
    pos = m.group(2)
    var = m.group(3).upper()
    if len(var) == 1:
        var = codon1to3[var]

    sav_list.append([wt, pos, var])
    wt_1codon = codon3to1[wt]
    var_1codon = codon3to1[var]
    sav_1codon_list.append([wt_1codon, pos, var_1codon])
    k = wt + pos + var
    sav_dict[k] = []
    rawDict[line] = k
    if pos in posDict.keys():
        posDict[pos].append(k)
    else:
        posDict[pos] = [k]
    kk = wt_1codon + pos + var_1codon
    sav_1codon_dict[kk] = k
    kkk = pos + var_1codon
    evDict[kkk] = k
file.close()

#print(sav_list)
#print(sav_1codon_list)

oldfile = open(sys.argv[2])
uniprot = sys.argv[3]

found = False
chain = ''
for line in oldfile:
    line=line.strip().split()
    #print (line[3],line[5], wt,res)
    try:
        if line[3]==wt and line[5]==pos:
            found=True
            chain=line[4]
            break
    except:
        pass
oldfile.close()
if not found:
    print('No Matched Position in PDB File.')
    exit()
oldfile=open(sys.argv[2],'r')
pdbfile = '{}/{}_clean.pdb'.format(workingDir,uniprot)
newfile=open(pdbfile,'w')
for line in oldfile:
    line=line.strip()
    temp=line
    temp=temp.split()
    try:
        if temp[4]==chain:
            newfile.write(line)
            newfile.write('\n')
    except:
        pass
oldfile.close()
newfile.close()

mutant_file = open('{}/individual_list_{}.txt'.format(workingDir, uniprot), 'w')
for sav in sav_1codon_list:
    mutant_file.write('{}{}{}{};\n'.format(sav[0],chain,sav[1],sav[2]))
mutant_file.close()

################################################################################
'''
Run evolutionAnalysis.jl
'''
################################################################################
os.system('julia evolutionAnalysis.jl {}/{} {} {}'.format(currDir, workingDir, sys.argv[7], proteinAccession))

################################################################################
'''
find variation number and EVMutation
'''
################################################################################
def write_fasta(file, accession, sequence):
    file.write(">{}\n".format(accession))
    if len(sequence) <= 70:
        file.write("{}\n".format(sequence))
    else:
        line_num_exact = len(sequence) / 70
        line_num = math.floor(line_num_exact)
        for i in range(0,line_num):
            start = i * 70
            stop = i * 70 + 70
            file.write(sequence[start:stop])
            file.write("\n")


        start = line_num * 70
        stop = len(sequence)
        file.write(sequence[start:stop])
        file.write("\n")

outputDir = '{}/{}'.format(currDir, workingDir)
################################################################################
with open('{}/{}_aligned.fasta'.format(outputDir,proteinAccession)) as f:
    alist = [line.rstrip() for line in f]
seqDict = {}

accession = ''
seq = ''
for line in alist:
    if '>' in line:
        if accession != '':
            seqDict[accession] = seq
        accession = line.replace('>', '')
        accession = accession.replace('_', '')
        accession = accession.replace('.', '')
        seq = ''
    else:
        seq = seq + line
seqDict[accession] = seq

homoAccession = proteinAccession
homoAccession = homoAccession.replace('_', '')
homoAccession = homoAccession.replace('.', '')
homoSeq = seqDict[homoAccession]
seqLength = len(homoSeq)

#################################################################################
#Phylogenetic tree
################################################################################
Phylo.convert('{}/{}_trees.nex'.format(outputDir, proteinAccession), 'nexus', '{}/{}_tree.tree'.format(outputDir, proteinAccession), 'newick')

f = open('{}/{}_tree.tree'.format(outputDir, proteinAccession), 'r')
tree = f.readline()
tree = re.sub(r':\d.\d+', '', tree)
tree = TreeNode.read(StringIO(tree))

####################################################################
#find set of characters given a list of taxa l and position number i
####################################################################
def findSet(l, i, seqDict):
    ll = []
    for j in l:
        ll.append(str(seqDict[j][i]))
    ll = list(set(ll))
    return ll

def updateVN(node, child, variation_number, seqDict, length):
    allClades = re.findall(r'[a-zA-Z0-9]+', str(node))
    for c in child:
        cClades = re.findall(r'[a-zA-Z0-9]+', str(c))
        subClades = list(set(allClades) - set(cClades))
        for i in range(length):
            cSet = findSet(cClades, i, seqDict)
            subSet = findSet(subClades, i, seqDict)
            for item in cSet:
                if item not in subSet:
                    variation_number[i] = variation_number[i] + 1
    return variation_number

def generateVN(tree, seqDict, seqLength):
    variation_number = np.zeros((seqLength,), dtype=int)
    queue = []
    queue.append(tree)

    while len(queue) > 0:
        node = queue.pop()
        if node.is_tip():
            continue
        else:
            child = node.children
            variation_number = updateVN(node, child, variation_number, seqDict, seqLength)
            for c in child:
                queue.append(c)
    return variation_number

variation_number = generateVN(tree, seqDict, seqLength)

homoIndexList = []
f_vn = []
for i in range(len(homoSeq)):
    if str(homoSeq[i]) != '-':
        homoIndexList.append(i)
        f_vn.append(variation_number[i])

outputFile = open("{}/vn_{}.txt".format(outputDir,proteinAccession), 'w')

vn_max = max(f_vn)
vn_min = min(f_vn)

vn_score = 0
for i in range(len(homoIndexList)):
    j = i + 1
    vn = variation_number[homoIndexList[i]]
    vn = (vn - vn_min) / (vn_max - vn_min)
    outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))
    if str(j) in posDict.keys():
        for k in posDict[str(j)]:
            sav_dict[k].append(vn)
outputFile.close()

################################################################################
outputFile = open("{}/evmutation_{}.txt".format(outputDir,proteinAccession), "w")
seq = seqDict[homoAccession]
seq1 = ''
for i in homoIndexList:
    seq1 = seq1 + seq[i]
write_fasta(outputFile, homoAccession, seq1)

for key in seqDict.keys():
    if key == homoAccession:
        continue
    seq = seqDict[key]

    seq1 = ''
    for i in homoIndexList:
        seq1 = seq1 + seq[i]
    write_fasta(outputFile, key, seq1)
outputFile.close()
################################################################################
'''
use rhapsody to find polyphen2
'''
################################################################################
test_SAVs = []#['{} {} {} {}'.format(uniprot, pos, wt_1codon, var_1codon)]
for sav in sav_1codon_list:
    test_SAVs.append('{} {} {} {}'.format(uniprot, sav[1], sav[0], sav[2]))
rh = rd.Rhapsody()
rh.queryPolyPhen2(test_SAVs)

polyphen2_file = 'pph2-full.txt'
file = open(polyphen2_file)
d_score = ''
s_score = ''
for line in file:
    line = line.rstrip()
    if len(line) == 0:
        continue
    line = line.split('\t')
    for i in range(len(line)):
        line[i] = line[i].replace(' ', '')
    if uniprot in line[0]:
        k = line[2] + line[1] + line[3]
        if k in sav_1codon_dict.keys():
            sav_dict[sav_1codon_dict[k]].append(float(line[22]))
            sav_dict[sav_1codon_dict[k]].append(float(line[23]))

file.close()
################################################################################
'''
calculate evmutation
'''
################################################################################
existingFiles = os.listdir(workingDir)
evmutationfile = '{}.params'.format(proteinAccession)

ll = (len(homoIndexList) - 1) * 0.2
#print(ll)

if evmutationfile in existingFiles:
    pass
else:
    os.system('plmc-master/bin/plmc -o {}/{}.params -le {} -lh 0.01 -m 100 {}/evmutation_{}.txt'.format(workingDir, proteinAccession, ll, workingDir, proteinAccession))

c = CouplingsModel('{}/{}.params'.format(workingDir, proteinAccession))
singles = single_mutant_matrix(c, output_column='effect_prediction_epistatic')

print('performaing EVMutation Analysis')
pd.DataFrame(singles).to_csv('{}/evmutation_{}_result.txt'.format(workingDir, proteinAccession))

evfile = '{}/evmutation_{}_result.txt'.format(workingDir, proteinAccession)
file = open(evfile)
ev_score = ''
for line in file:
    line = line.rstrip()
    if len(line) == 0:
        continue
    line = line.split(',')
    #print(line)
    k = line[3] + line[5]
    if k in evDict.keys():
        sav_dict[evDict[k]].append(float(line[8]))
file.close()
################################################################################
'''
process fathmm score
'''
################################################################################
fathmm_file = sys.argv[9]
file = open(fathmm_file)
for line in file:
    line = line.rstrip()
    line = line.split('\t')
    #print(line)
    if line[3] in sav_1codon_dict.keys():
        sav_dict[sav_1codon_dict[line[3]]].append(float(line[5]))
file.close()
################################################################################
'''
run FoldX
'''
#print(pdbfile)
################################################################################
existingFiles = os.listdir(workingDir)
repairPDBfile = '{}_clean_Repair.pdb'.format(uniprot)
copyfile(pdbfile, '{}/{}_clean.pdb'.format(currDir, uniprot))
if repairPDBfile in existingFiles:
    print('###Repaired PDB file exists###')
else:
    print('###Using FoldX to generate Repaired PDB file###')
    os.system('./foldx --command=RepairPDB --pdb={} --output-dir={}'.format(pdbfile,workingDir))
os.remove('{}/{}_clean.pdb'.format(currDir, uniprot))
print('###Using FoldX to calculate free energy###')
foldx_out = '{}/Dif_{}_clean_Repair.fxout'.format(workingDir,uniprot)

try:
    os.remove(foldx_out)
except:
    pass

copyfile('{}/{}'.format(workingDir, repairPDBfile), '{}/{}'.format(currDir, repairPDBfile))
os.system('./foldx --command=BuildModel --pdb={}/{} --mutant-file={}/individual_list_{}.txt --output-dir={} --out-pdb=false'.format(workingDir,repairPDBfile,workingDir,uniprot, workingDir))
os.remove('{}/{}'.format(currDir, repairPDBfile))
file = open(foldx_out)
foldx_free_energy = ''
index = 0
for line in file:
    line = line.rstrip()
    if len(line) == 0:
        continue
    line = line.split()
    if uniprot in line[0]:
        sav = sav_list[index]
        k = sav[0] + sav[1] + sav[2]
        sav_dict[k].append(float(line[1]))
        index += 1


file.close()
################################################################################
'''
run freeSASA
'''
################################################################################
freeSASAfile = '{}/freesasa_{}.txt'.format(workingDir,uniprot)
os.system('freesasa --format=seq --shrake-rupley -n 200 --probe-radius 1.2 --n-threads 4 {} > {}'.format(pdbfile, freeSASAfile))
file = open(freeSASAfile)
freesasa_score = ''
count = 0
for line in file:
    line = line.rstrip()
    if len(line) == 0:
        continue
    line = line.split()
    if count == 0:
        count = 1
        continue
    if line[2] in posDict.keys():
        for k in posDict[line[2]]:
            sav_dict[k].append(float(line[5]))

file.close()
################################################################################
'''
run MAESTRO
'''
################################################################################
maestroDir = sys.argv[5]
for sav in sav_1codon_list:
    mut = '{}{}{{{}}}'.format(sav[0],sav[1],sav[2])
    maestro_file = '{}/{}_maestro.txt'.format(workingDir, uniprot)
    os.system('{}/maestro {}/config.xml {} --evalmut=\'{}\' --bu > {}'.format(maestroDir, maestroDir, pdbfile, mut,maestro_file))
    file = open(maestro_file)
    maestro_score = ''
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
        line = line.split()
        if sav[1] in line[2]:
            k = sav[0] + sav[1] + sav[2]
            sav_dict[sav_1codon_dict[k]].append(float(line[3]))
            break
    file.close()
################################################################################
'''
run prody
ANM, Mechanical stiffness, Effectiveness, Sensitivity
'''
################################################################################
#ANM
filename = pdbfile
prot=prody.parsePDB(filename)
anm,sel = prody.calcANM(prot)
flucts=prody.calcSqFlucts(anm[0])

pdbIndexDict = {}
file = open(pdbfile)
index = -1
currPos = 0
count = 0
for line in file:
    line = line.strip()
    line = line.split()
    if count == 0:
        index = 0
        currPos = int(line[5])
        pdbIndexDict[currPos] = index
        count = 1
        continue
    if currPos == int(line[5]):
        continue
    index += 1
    currPos = int(line[5])
    pdbIndexDict[currPos] = index
file.close()

for pos in posDict.keys():
    for k in posDict[pos]:
        sav_dict[k].append(flucts[pdbIndexDict[int(pos)]])

#mechanical_stiffness
gfp=prody.parsePDB(filename)
calphas=gfp.ca
anm = prody.ANM('ANM analysis')
anm.buildHessian(calphas, cutoff=13.0)
anm.calcModes(n_modes='all')
stiffness = prody.calcMechStiff(anm, calphas)
ms=[]
for row in stiffness:
    ms.append(str(sum(row)/len(row)))

for pos in posDict.keys():
    for k in posDict[pos]:
        sav_dict[k].append(float(ms[pdbIndexDict[int(pos)]]))

#Effectieness
ampar_ca = prody.parsePDB(filename, subset='ca')

anm_ampar = prody.ANM('AMPAR')
anm_ampar.buildHessian(ampar_ca)
anm_ampar.calcModes()
prs_mat, eff, sen = prody.calcPerturbResponse(anm_ampar)
eff=eff.tolist()
sen=sen.tolist()
for pos in posDict.keys():
    eff1=eff[pdbIndexDict[int(pos)]]
    sen1=sen[pdbIndexDict[int(pos)]]
    for k in posDict[pos]:
        sav_dict[k].append(eff1)
        sav_dict[k].append(sen1)
################################################################################
'''
pyrosetta
'''
################################################################################
def pyrosettaScore(chain, pos, var_1codon):
    init()
    pose=pose_from_pdb(filename)
    scorefxn=get_fa_scorefxn()
    scorefxn.set_weight(fa_atr,0)
    scorefxn.set_weight(fa_rep,0)
    scorefxn.set_weight(fa_intra_rep,0)
    scorefxn.set_weight(fa_sol,0)
    scorefxn.set_weight(lk_ball_wtd,0)
    scorefxn.set_weight(fa_intra_sol,0)
    scorefxn.set_weight(fa_elec,0)
    scorefxn.set_weight(hbond_lr_bb,0)
    scorefxn.set_weight(hbond_sr_bb,0)
    scorefxn.set_weight(hbond_bb_sc,0)
    scorefxn.set_weight(hbond_sc,0)
    scorefxn.set_weight(dslf_fa13,0)
    scorefxn.set_weight(rama_prepro,0)
    scorefxn.set_weight(p_aa_pp,0)
    scorefxn.set_weight(fa_dun,0)
    scorefxn.set_weight(omega,0)
    scorefxn.set_weight(pro_close,0)
    scorefxn.set_weight(yhh_planarity,0)
    scorefxn.set_weight(fa_intra_sol_xover4,0)
    wt_score=scorefxn(pose)
    #print (scorefxn.show(pose))
    mutate_residue(pose,pose.pdb_info().pdb2pose(chain,int(pos)),var_1codon)
    mutant_score=scorefxn(pose)
    diff = mutant_score - wt_score
    return mutant_score, diff

for pos in posDict.keys():
    for k in posDict[pos]:
        kk = codon3to1[k[-3:]]
        m,d = pyrosettaScore(chain, pos, kk)
        sav_dict[k].append(m)
        sav_dict[k].append(d)

################################################################################
'''
P2Rank
'''
################################################################################
p2rankDir = sys.argv[6]
p2rank_file = '{}/{}_clean.pdb_residues.csv'.format(workingDir, uniprot)

try:
    os.remove(p2rank_file)
except:
    pass
os.system('{}/prank predict -f {}'.format(p2rankDir, pdbfile))
shutil.move('{}/test_output/predict_{}_clean/{}_clean.pdb_residues.csv'.format(p2rankDir, uniprot, uniprot), workingDir)

file = open(p2rank_file)
p2rank_score = 0
for line in file:
    line = line.rstrip()
    if len(line) == 0:
        conitnue
    line = line.replace(' ', '')
    line = line.split(',')
    if line[1] in posDict.keys():
        p2rank_score = 1
        if float(line[5]) < 0.5:
            p2rank_score = 0
        for k in posDict[line[1]]:
            sav_dict[k].append(p2rank_score)

file.close()
################################################################################
'''
generate data array
'''
################################################################################
from joblib import dump, load
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import StandardScaler
from tpot.builtins import StackingEstimator
from sklearn.preprocessing import FunctionTransformer
from copy import copy
from matplotlib import pyplot
import math
from xgboost import XGBClassifier

xgboost = load('XGBoostClassifier.joblib')
print('####################################################################')
print('#')
print('#Prediction written to file: {}/{}_xgboost_prediction.csv'.format(workingDir,uniprot))
print('#')
print('####################################################################')

outputFile = open('{}/{}_xgboost_prediction.csv'.format(workingDir,uniprot), 'w')
outputFile.write('uniprot,accession,prediction_proba,prediction,variation_nuumber,dScore,Score1,evmutation,FATHMM,FoldX,SASA,maestro,ANM,MS,effectiveness,sensitivity,pyrosetta_mutant,pyrosetta_diff,active_site\n')
for sav in sav_dict.keys():
    data = sav_dict[sav]
    #print(data)
    data = np.array([data])
    proba = xgboost.predict_proba(data)[0][1]
    p = xgboost.predict(data)[0]
    outputFile.write('{},{},{},{}'.format(uniprot, proteinAccession,proba,p))
    for k in sav_dict[sav]:
        outputFile.write(',{}'.format(k))
    outputFile.write('\n')

outputFile.close()
