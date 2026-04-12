<img src="img/pathogen-Havey-igh-twitter-header1.jpg" width="100%" />

# Microbial Genomics In-Class Activity

### DGP-485 Data Science for Biomedical Researchers
**April 9, 2025**

*Egon A. Ozer, MD PhD (<e-ozer@northwestern.edu>)*  
 
---

### Introduction 
In this exercise we will assemble the whole genome sequence of an isolate of the bacterium *Staphyloccus aureus* from Illumina short reads, assess the assembly quality, identify antibiotic resistance markers, and determine variants relative to a reference genome sequence.

We will perform all these tasks on Quest using Conda / Mamba. For a nice introduction to Conda, see [here](https://www.dataschool.io/intro-to-conda-environments/) or [here](https://docs.conda.io/projects/conda/en/latest/index.html) for more detail. For more information about using Mamba or Conda on Quest specifically, take a look at [this page](https://services.northwestern.edu/TDClient/30/Portal/KB/ArticleDet?ID=2064).

> <img src="img/warn.png" width="25"> You can also install Mamba or Conda on your personal computer (all of the analyses in this exercise can easily be done on a laptop). I recommend Mamba or Micromamba as the installation of software packages into environments with Mamba tends to be much faster and smoother than with Conda or Anaconda. See [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for more information on installing Mamba or Micromamba. If you install Conda and/or Mamba on you computer you do not need to load it as a module.  

### Software to install: 

#### Remote file broswer

File browsers allow you to access and browse remote filesystems (like Quest) and drag and drop files or folders back and forth. 

**Options:**

<img src="img/cyberduck-icon-384.png" width="40"/>[Cyberduck](https://cyberduck.io/) (This one is pretty, but will ask you for money)

<img src="img/95988_filezilla_icon.png" width="40">[FileZilla](https://filezilla-project.org/) (This one is free)

### 1. Logging into Quest and activating the environment  

To log into quest, open your terminal program and type

```
ssh <NETID>@login.quest.northwestern.edu
```

Where `<NETID>` is replaced with your NetID (such as `abc123` for example).

Then enter your NetID password and hit Enter.

**a. Create a working folder**

1. To get to the project folder, use the `cd` (change directory) command followed by the path to the exercise folder

```
cd /projects/e30682/MicrobialWGS
```

2. In this folder, create a folder for yourself with the `mkdir` (make directory) command.

Again: `<NETID>` should be replaced with your NetID

```
mkdir <NETID>
```

3. Move into your folder:

```
cd <NETID>
```

4. Ensure you are in your folder by using the `pwd` (present working directory) command

```
pwd
```

The output of the command should be `/projects/e30682/MicrobialWGS/<NETID>`


**b. Load the Mamba module on Quest and activate the pre-made environment** 

1. This first command loads Quest's `mamba` module.

```
module load mamba
```  

2. Now activate the pre-made mamba environment called `assembly_env`

```
conda activate /projects/e30682/MicrobialWGS/conda_envs/assembly_env
```

> <img src="img/warn.png" width="25"> If this is the first time you've used Mamba or Conda on Quest, you'll probably get an error message here. The following commands will get it working for you, after which you can rerun the command *b2* above. You'll only need to run these two commands once and Conda/Mamba should work for you every time you use them on Quest going forward.
> ```
> conda init bash
> source ~/.bashrc
> ```

> <img src="img/warn.png" width="25"> If you want to recreate this environment in your home directory or on your own computer, you can copy the `assembly_environment.yaml` file from the `/projects/e30682/MicrobialWGS/conda_envs` directory and create the environment using the `mamba env create -f assembly_environment.yaml` command. 


### 2. Perform quality trimming of the sequencing reads

Before we assemble, we'll do some read trimming. This step removes any Illumina adapter sequences that may have made it into the reads due to read-through of short library fragments. These adapter sequences can sometimes get incorporated into your assembly and cause misassemblies. Better to take the time to remove them prior to assembly. 

We are going to use [fastp](https://github.com/OpenGene/fastp) to perform the trimming step. Fastp can perform a number of functions including adapter removal as well remove low-quality sequences from the reads. 

_Commands_ 

```Shell
fastp \
    --in1 /projects/e30682/MicrobialWGS/reads/SA_1.fastq.gz \
    --in2 /projects/e30682/MicrobialWGS/reads/SA_2.fastq.gz \
    --out1 SA_trimmed_paired_1.fastq.gz --unpaired1 SA_trimmed_unpaired_1.fastq.gz \
    --out2 SA_trimmed_paired_2.fastq.gz --unpaired2 SA_trimmed_unpaired_2.fastq.gz \
    -h SA_fastp.html -j SA_fastp.json \
    -w 1
```    
_Settings Used_

Setting | Descripton
--- | ---
`--in1 / --in2` | Input read files
`--out1 / --out2` |  Paired output read files
`--unpaired1 / --unpaired2` | Singleton read files (only one of the paired reads passed filters)  
`-h` | Filtering and trimming report, html format
`-j` | Filtering and trimming reoprt, json format
`-w` | Number of parallel threads to use
   

See [fastp manual](https://github.com/OpenGene/fastp/blob/master/README.md) for more detail on settings and other options.

_Outputs_

Files | Description
--- | ---
`SA_trimmed_paired_1.fastq.gz` & `_2.fastq.gz` | Paired reads remaining after trimming
`SA_trimmed_unpaired_1.fastq.gz` & `_2.fastq.gz` | Unpaired reads remaining after trimming 
`SA_fastp.html` | Filtering and trimming report. Can view with a web browser

## 3. Assemble with SPAdes

Now we'll generate a _de novo_ whole genome assembly from our trimmed reads. For this we'll use the assembler [SPAdes](https://cab.spbu.ru/software/spades/) ([Github site](https://github.com/ablab/spades)).  


> <img src="img/warn.png" width="25"> This step takes too many resources and is too slow to do in the login node. We'll submit a Slurm batch script for this instead.

```
nano spades.sh
```

Copy and paste the below code into the file you just opened on Quest using nano

*spades.sh*

```Shell
#!/bin/bash
#SBATCH --job-name="spades"
#SBATCH -A e30682
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o spades.out
#SBATCH -e spades.err

eval "$(conda shell.bash hook)"
conda activate /projects/e30682/MicrobialWGS/conda_envs/assembly_env

spades.py \
    -o SA_spades \
    -1 SA_trimmed_paired_1.fastq.gz \
    -2 SA_trimmed_paired_2.fastq.gz \
    --isolate \
    -t 12 --cov-cutoff auto \
    | tee SA_spades_output.txt

```

To close and save the document, first type `Ctrl + x` to exit, then `y`, then hit `Enter` to save the file with the same name.

To submit the job request, type: 

```
sbatch spades.sh
```

To check if your job is queued, running, or finished, you can use the following command:

```
squeue --me
```

**Outputs**

All of the output files can be found in the `SA_spades` folder. There are a lot of files in there, so we're just going to pick a few of the most relevant to describe.

Files | Description
--- | ---
`contigs.fasta` | Assembly contig sequences.
`scaffolds.fasta` | Scaffold sequnences. Scaffolds consist of contigs joined by N's at sites with paired reads support.
`spades.log` | Log file. Useful for troubleshooting poor or failed assemblies.
`assembly_graph.fastg` | Graph of assembly showing contig connections. Can be useful for assembly quality assessment using [Bandage](https://rrwick.github.io/Bandage/). 


## 4. Assembly quality assessment - part 1: CheckM

First we'll use [CheckM](https://ecogenomics.github.io/CheckM/) to estimate assembled genome completeness and contamination. Again this is a slow one so we'll submit it as a Slurm job.

```
nano checkm.sh
```

*checkm.sh*

```Shell
#!/bin/bash
#SBATCH --job-name="checkm"
#SBATCH -A e30682
#SBATCH -p short
#SBATCH -t 00:15:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o checkm.out
#SBATCH -e checkm.err

eval "$(conda shell.bash hook)"
conda activate /projects/e30682/MicrobialWGS/conda_envs/assembly_env

echo -e "SA\tSA_spades/contigs.fasta" > SA.checkm_source.txt
export CHECKM_DATA_PATH=/projects/e30682/MicrobialWGS/ref/checkm_database

checkm lineage_wf \
  -t 12 \
  --tab_table \
  -f SA.checkm_results.txt \
  SA.checkm_source.txt \
  SA_checkm
```

```
sbatch checkm.sh
```

Results will be summarized in the `SA.checkm_results.txt` file. The ones to pay attention to are the *Completeness*, *Contamination*, and *Strain heterogeneity* values. 

From the [CheckM wiki](https://github.com/Ecogenomics/CheckM/wiki/Reported-Statistics#qa):

* completeness: estimated completeness of genome as determined from the presence/absence of marker genes and the expected collocalization of these genes.
* contamination: estimated contamination of genome as determined by the presence of multi-copy marker genes and the expected collocalization of these genes.
* strain heterogeneity: estimated strain heterogeneity as determined from the number of multi-copy marker pairs which exceed a specified amino acid identity threshold (default = 90%). High strain heterogeneity suggests the majority of reported contamination is from one or more closely related organisms (i.e. potentially the same species), while low strain heterogeneity suggests the majority of contamination is from more phylogenetically diverse sources.

## 5. Assembly quality assessment - part 2: Quast

Running CheckM can take up to 10 minutes to complete, so while that is running we'll perform some other measures of assembly quality. Assembly quality assement includes, but is not limited to, determination of the total number of contigs or scaffolds, the total length of the sequence, and the GC content. 

* Low numbers of contigs usually indicate a pretty good assembly where many unambiguous connections could be made. High numbers of contigs may indicate low read coverage or possible genome or read contamination from another source (though not always)
* Total contig length (sum of the lengths of all contigs) should be compared to the expected genome size of the species or organism. [NCBI Genome](https://www.ncbi.nlm.nih.gov/genome/) is a good source of expected genome sizes. 
* Percent GC content (i.e. the percentage of the sequence that is either cytosine or guanine bases) should also be compared to the expected genome GC content. 
* N50 is length of the contig such that contigs that length or longer account for >= 50% of the total assembly size. A small N50 value may indicate a more fragmented assembly.

_Commands_

```Shell
quast SA_spades/contigs.fasta \
  -r /projects/e30682/MicrobialWGS/ref/NCTC_8325.fasta \
  -t 1
```

_Settings_

_Settings Used_

Setting | Descripton
--- | ---
`SA_spades/contigs.fasta` | Path to the input assembly file
`-r` |  Reference genome sequence file (OPTIONAL, see * below)
`-t` | Number of parallel threads to use

\* Using a reference genome for comparision with the `-r` setting is optional, but can be used to determine possible misassemblies. Just be warned that if your reference genome sequence is not very closely related, you may end up with a number of false-positive missassembly warnings. _Caveat emptor_. 

_Outputs_

By default, outputs will be put into a directory called "quast_results". Here are a few of the output files worth noting:

Files | Description
--- | ---
`report.html` | A report file that can be opened in a web browser like Chrome or Safari.
`report.pdf` | A pdf version of the assembly report
`report.tsv` | A tab-separated version of the report table that can be viewed in a spreadsheet program like Excel or in a text editor

## 6. Identify antimicrobial resistance genes and mutations with AMRFinder

[AMRFinder](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) is a tool and database developed by NCBI to identify antimicrobial resistance (AMR) genes, resistance-associated point mutations, and select other classes of genes from assembled nucleotide sequence. 

_Commands_

```Shell
amrfinder -n SA_spades/contigs.fasta \
  -O Staphylococcus_aureus \
  --plus \
  --threads 1 \
  --output SA_amrfinder.txt
```

_Settings Used_

Setting | Descripton
--- | ---
`-n` | Path to the input assembly nucleotide sequence file
`-O` | Taxonomy group of the organism 
`--plus` | Add virulence genes to the report
`--threads` | Number of parallel threads to use
`--output` | File to output results to 

_Outputs_

Files | Description
--- | ---
`SA_amrfinder.txt` | A report file that can be downloaded and opened in Excel

[More detail about output format](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format)

## 7. Align reads to reference genome and identify variants

[Snippy](https://github.com/tseemann/snippy) is a nice "all-in-one" pipeline for generating alignments and using those alignments to determine variants. If you supply a genbank file as the reference sequence, the program will also annnotate the variants, i.e. call synonymous, non-synonymous, or frameshift. Will also output consensus sequences. 

**Commands**

```
nano snippy.sh
```

*snippy.sh*

```Shell
#!/bin/bash
#SBATCH --job-name="snippy"
#SBATCH -A e30682
#SBATCH -p short
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o snippy.out
#SBATCH -e snippy.err

eval "$(conda shell.bash hook)"
conda activate /projects/e30682/MicrobialWGS/conda_envs/alignment_env

snippy \
    --outdir SA_snippy \
    --reference /projects/e30682/MicrobialWGS/ref/NCTC_8325.gbk \
    --R1 /projects/e30682/MicrobialWGS/reads/SA_1.fastq.gz \
    --R2 /projects/e30682/MicrobialWGS/reads/SA_2.fastq.gz\
    --cpus 12
```

```
sbatch snippy.sh
```

**Settings**

Setting | Description
--- | ---
`--outdir` | Name of the directory the output files will go into. You'll get an error if this directory already exists
`--reference` | Genome sequence to use as an alignment reference. Can be fasta or genbank file.
`--R1` & `--R2` | Paired-end read files
`--cpus` | Number of compute cores to use.

For more detail on settings, visit <https://github.com/tseemann/snippy> or run `snippy --help`

**Outputs**

_Adapted from <https://github.com/tseemann/snippy>. Follow the link to see the full list of output files._

Extension | Description
----------|--------------
.tab | A simple [tab-separated](http://en.wikipedia.org/wiki/Tab-separated_values) summary of all the variants
.html | A [HTML](http://en.wikipedia.org/wiki/HTML) version of the .tab file for viewing in a web browser like Chrome or Safari
.bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Includes unmapped, multimapping reads. Excludes duplicates.
.bam.bai | Index for the .bam file
.aligned.fa | A version of the reference but with `-` at position with `depth=0` and `N` for `0 < depth < --mincov` (**Note: snippy manual says this file "does not have variants," but in the current version 4.6.0 it does incorporate the single nucleotide variants**)
.consensus.fa | A version of the reference genome with *all* variants instantiated (substitutions and indels), but no masking of positions under minimum depth
.consensus.subs.fa | A version of the reference genome with *only substitution* variants instantiated, but no masking of positions under minimum depth
.log | A log file with the commands run and their outputs

**When you are all done...**
Deactivate the conda environment. 

```
conda deactivate
```

---

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="img/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
