################################################################################
#required packages:
#   clustal omega
#   paup
#   python: skbio, Bio, scipy
#
#ARGS:
#1: input/output dir
#2: homolog sequence file
#3: target sequence fasta identification
#4: protein or dna; default dna
#
#Example:
#julia evolutionAnalysis.jl FULL_PATH/example APP_refseq_protein.fasta NP_000475.1 protein
#
################################################################################

using Random
using CSV
################################################################################
#genbank_dict:
#    store accession as key and organism as value
################################################################################
dir = ARGS[1]
outputDir = dir

function replaceAA(sequence)
    seq = ""
    for char in sequence
        if occursin(char,"ARNDCQEGHILKMFPSTWYV")
            seq = seq * char
        else
            seq = seq * "X"
        end
    end
    return seq
end

#seq = "AACJ"
#println(replaceAA(seq))

function write_fasta(file, accession, sequence)
    write(file, ">$accession\n")
    if length(sequence) <= 70
        write(file, "$sequence\n")
    else
        line_num_exact = length(sequence) / 70
        line_num = floor(Int, line_num_exact)
        #println(line_num)
        for i in 0:line_num - 1
            start = i * 70 + 1
            stop = i * 70 + 70
            write(file, sequence[start:stop], "\n")
        end

        start = line_num * 70 + 1
        stop = length(sequence)
        write(file, sequence[start:stop], "\n")
        write(file, "\n")
    end
end

################################################################################
#read fasta file
################################################################################
in_file = open("$dir/$(ARGS[2])", "r")
println("reading $(ARGS[2]) file")
accession = ""
sequence = ""
fastaDict = Dict()
for line in eachline(in_file)
    if length(line) == 0
        continue
    end
    global accession
    global sequence
    #println(line)
    if occursin(">", line)
        #println(line)
        if accession == ""
            accession = match(r">([A-Za-z0-9_]+.[0-9]+)", line)[1]
            continue
        end
        sequence = replaceAA(sequence)
        fastaDict[accession] = sequence
        sequence = ""
        accession = match(r">([A-Za-z0-9_]+.[0-9]+)", line)[1]
    else
        sequence = sequence * line
    end
end
sequence = replaceAA(sequence)
fastaDict[accession] = sequence
close(in_file)

################################################################################
#write fasta file
################################################################################

file = open("$dir/$(ARGS[3]).fasta", "w")
for key in keys(fastaDict)
    write_fasta(file, key, fastaDict[key])
end
close(file)

################################################################################
#run cluster omega
################################################################################
files = Base.readdir("$dir")
println("running CLUSTAL OMEGA for multiple sequence alignment")
#run(`clustalo -i $outputDir/$(ARGS[3]).fasta -o $outputDir/$(ARGS[3])_aligned.fasta --auto -v --force`)

if "$(ARGS[3])_aligned.fasta" in files
    println("$(ARGS[3])_aligned.fasta found")
else
    #println("running CLUSTAL OMEGA for multiple sequence alignment")
    run(`clustalo -i $dir/$(ARGS[3]).fasta -o $dir/$(ARGS[3])_aligned.fasta --auto -v --force`)
end

################################################################################
#convert fasta identifier name
################################################################################
file = open("$dir/$(ARGS[3])_aligned.fasta", "r")
seqDict = Dict()

accession = ""
seq = ""
for line in eachline(file)
    global accession
    global seq
    if isempty(line)
        continue
    end
    if occursin(">", line)
        if length(accession) != 0
            seqDict[accession] = seq
            seq = ""
        end
        #l = split(line, " ")
        #println(l[1])
        accession = match(r">([A-Za-z0-9_]+.[0-9]+)", line)[1]
        accession = replace(accession, "_" => "")
        accession = replace(accession, ">" => "")
	    accession = replace(accession, "." => "")
        #println(accession)
    else
        #println(line)
        seq = seq * line
    end
end
seqDict[accession] = seq
nchar = length(seq)
close(file)
################################################################################
#convert fasta to nexus file
################################################################################
function writeNexus(fileName, seqDict, sequencetype)
    nexus_file = open(fileName, "w")
    println("generating $nexus_file nexus file")
    write(nexus_file, "#NEXUS\n\n")
    write(nexus_file, "BEGIN TAXA;\n")
    write(nexus_file, "\tTITLE Taxa;\n")
    write(nexus_file, "\tDIMENSIONS NTAX=$(length(seqDict));\n")

    taxa = join(sort(collect(keys(seqDict))), " ")
    dimension = 0
    for key in keys(seqDict)
        dimension = length(seqDict[key])
        break
    end

    write(nexus_file, "\tTAXLABELS\n")
    write(nexus_file, "\t\t$taxa\n\t;\n\n")
    write(nexus_file, "END;\n\n")

    write(nexus_file, "BEGIN CHARACTERS;\n")
    write(nexus_file, "\tTITLE Character_Matrix;\n")
    write(nexus_file, "\tDIMENSIONS NCHAR=$dimension;\n")
    write(nexus_file, "\tFORMAT DATATYPE = $(sequencetype) GAP = - MISSING = ?;\n")
    write(nexus_file, "\tMATRIX\n")

    for key in sort(collect(keys(seqDict)))
        write(nexus_file, "\t$key $(seqDict[key])\n")
    end

    write(nexus_file, "\n;\nEND;")
    close(nexus_file)
end

nexusFile = "$dir/$(ARGS[3]).nex"
seqtype = "protein"

#println(ARGS[4])
#println(seqtype)
#exit()
writeNexus(nexusFile, seqDict, seqtype)

################################################################################
#run paup
################################################################################
function getTreeCMD(nexusFileName, outputFile, dir)
    #=
    DEFINE SETTINGS
    =#

    #DESIGNATED OUTGROUP
    #defineOutgroup = "reference /only"

    #LOCATION OF PAUP
    #paupApp = "/gpfs/runtime/opt/paup/4.0a166/bin/paup"

    #TREE SEARCH METHOD
    #tree search method for simultaneous analysis tree & support tests
    treeSearchMethod = "hsearch nreps=1000 swap=tbr multrees=no"

    #BOOTSTRAP TREE SEARCH METHOD
    #   specify method for bootstrap tree searhc
    #   bootstrapMethod = "search=bandb"
    bootstrapMethod = "search=heuristic nreps=50"

    createTreeCmdFile = open(outputFile, "w")
    write(createTreeCmdFile, "#NEXUS\n\n")
    write(createTreeCmdFile, "set warnReset = no;\n")
    write(createTreeCmdFile, "set increase = auto;\n")
    write(createTreeCmdFile, "set datastorage=full;\n")
    write(createTreeCmdFile, "set criterion=parsimony;\n")
    #perform tree search & calculate bootsraps for SA tree
    #write(createTreeCmdFile, "[SIMULTANEOUS ANALYSIS]\n")
    write(createTreeCmdFile, "execute $nexusFileName;\n")
    #write(createTreeCmdFile, "outgroup $defineOutgroup;\n")
    write(createTreeCmdFile, "$treeSearchMethod;\n")
    write(createTreeCmdFile, "filter best;\n")
    treeFile = replace(outputFile, "_getTree.cmd" => "_trees.nex")
    write(createTreeCmdFile, "savetrees file=$treeFile format=nexus replace=yes root=yes;\n")
    #write(createTreeCmdFile, "contree /strict=yes treeFile=SA_Consensus.nex replace=yes;\n")
    #write(createTreeCmdFile, "execute SA_Consensus.nex;\n")
    #write(createTreeCmdFile, "savetrees file=SA_Consensus.nex format=nexus replace=yes root=yes;\n")
    #write(createTreeCmdFile, "bootstrap $bootstrapMethod;\n")
    #write(createTreeCmdFile, "savetrees from=1 to=1 file=SA_Bootstrap.nex savebootp=nodelabels maxdecimals=0 format=nexus replace=yes;\n\n")

    write(createTreeCmdFile, "quit warnTsave=no;\n")
    close(createTreeCmdFile)
end
getTreeCMD(nexusFile, "$dir/$(ARGS[3])_getTree.cmd", dir)

#=
run paup
=#
run(`paup $dir/$(ARGS[3])_getTree.cmd`)

println("Generating Variation Number")
#run(`python3 caos.py $(ARGS[3]) $(ARGS[1])`)

exit()

#=
run mrbayes
=#
nexusFile = "$(ARGS[1])_mrbayes.nex"
function generateMrbayes(fileName, seqDict, sequencetype, nchar)
    nexus_file = open(fileName, "w")
    println("generating $nexus_file nexus file")
    write(nexus_file, "#NEXUS\n\n")
    write(nexus_file, "begin data;\n")
    write(nexus_file, "\tdimensions ntax=$(length(seqDict)) nchar=$nchar;\n")
    write(nexus_file, "\tformat datatype=$(sequencetype) interleave=no gap=-;\n")
    write(nexus_file, "\tmatrix\n")

    for key in sort(collect(keys(seqDict)))
        write(nexus_file, "$key $(seqDict[key])\n")
    end

    write(nexus_file, "\n;\nend;")
    close(nexus_file)

    mrbayesFile = open("mrbayes_batch.txt", "w")
    if sequencetype == "dna"
        write(mrbayesFile, "set autoclose=yes nowarn=yes\n")
        write(mrbayesFile, "execute $nexusFile\n")
        write(mrbayesFile, "lset nst=6 rates=gamma\n")
        write(mrbayesFile, "mcmc ngen=1000000 samplefreq=100 printfrq=500000 diagnfreq=10000\n")
        write(mrbayesFile, "sumt\n")
        write(mrbayesFile, "sump\n")
        write(mrbayesFile, "quit\n")

    else
        write(mrbayesFile, "set autoclose=yes nowarn=yes\n")
        write(mrbayesFile, "execute $nexusFile\n")
        write(mrbayesFile, "prset aamodelpr = mixed\n")
        write(mrbayesFile, "mcmc nchains = 1 ngen = 500000\n")
        write(mrbayesFile, "sumt\n")
        write(mrbayesFile, "sump\n")
        write(mrbayesFile, "quit\n")
    end
    close(mrbayesFile)

end
generateMrbayes(nexusFile, seqDict, sequencetype, nchar)
#run(`mb mrbayes_batch.txt`)
exit()
