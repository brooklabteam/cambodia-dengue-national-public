# cambodia-dengue-national-public

# cambodia-dengue-national

This repo contains directions and scripts for our Cambodia dengue paper, which can be accessed as a preprint at [Brook et al. 2024](https://www.medrxiv.org/content/10.1101/2022.06.08.22276171v5). All final figures can be accessed in the "final-figures" folder, and all data used to generate them can be accessed in the "data" subfolder. Within the "figure-development" folder, you will additionally find subfolders corresponding to each figure ("Fig1", "Fig2", "Fig3", "Fig4", "Fig5") in the manuscript which contain the scripts used to (a) generate any data included in the folder and (b) actually produce the figure itself.


Below, you will find detailed directions of the broad workflow used in the  preparation of our Bayesian timetrees in Figure .

---

## Workflow to Generate Bayesian Timetrees

---

### Gathering background sequences from GenBank

We first sought to supplement our existing Cambodia sequences with those available in GenBank [(NCBI Virus)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). We collected full or partial genome nucleotide sequences >10000 bp in length, selecting ALL existing sequences for Cambodia available for Dengue 1 up to a collection date of December 31, 2022 (tax id 11053; 197 sequences, including the 63 we contributed) and Dengue 2 (tax id 11060; 176 sequences, including the 120 we contributed). We supplemented these with genomes collected from other major Southeast Asian countries (nine), which were Laos, Myanmar, Malaysia, Thailand, Vietnam, Brunei, Indonesia, the Philippines, and Singapore. All countries were represented in both the DENV-1 and DENV-2 datasets. To avoid overrepresenting certain countries outside of Cambodia, we limited sequence selection to a maximum of three genomes collected per year from each available year per country, beginning in 2000, the earliest year in which any Cambodian sequences were recovered, and ending in 2022.

The subfolder, "gather-sequences", includes the original NCBI datafile used to select sequences ("allSEasiadengue.csv"), the code used to query this list and select the sub-select the background sequences ("query-genbank.R"), and the output file that lists the accession number, genotype, date, and country of origin for all sequences beyond those generated in our own sampling, which are used in this paper ("All_Seq_SE_Asia.csv"). The subfolder also includes two text files used to download these sequences as a batch .fasta file from GenBank ("DENV1-NCBI.txt" and "DENV2-NCBI.txt"), as well as the resulting .fasta batch files from those downloads ("GenBankDENV1.fasta" and "GenBankDENV1.fasta").

---

### Including our own sequences

We generated our own sequences for this analysis which we submitted to GenBank between 2021 and 2024 (DENV1 (63): accession numbers PP470671-PP470673, PP470602, OK159935-OK159976, OL411495-OL411499, and OL412140, OL412678,OL412703, OQ678009, OQ678018, OQ678060,  OQ678061, OQ678099, OQ678244, OQ678252, OQ683856 OQ683880; DENV2 (120): OL412740, OL414717-OL414730, OL414732-OL414747, OL414749-OL414765, OL420733, OL435143, OP999336, OQ674509, OQ678012, OQ678017, OQ678059, OQ678101, OQ678102, OQ683881, PP396279-PP396280, PP406368-PP406376, PP406447-PP406455, PP411227-PP411235, PP411249-PP411257, PP439574-PP439582, PP446688-PP446694, PP469533-PP469534, PP470598-PP470603; DENV4 (9):  MZ976858-MZ976860, OQ678015,PP396732, PP447197-PP447200. To generate these sequences, we (a) performed metagenomic Next Generation Sequencing on samples from patients reporting to clinic with undiagnosed fevers, then ran them through the [CZID](czid.org) pipeline to query the origin of each fever. For all samples positive for dengue, which generated reliable contigs in the metagenomic pipeline, we then sought to generate consensus genomes using the [ARTIC network's Nextflow consensus genome pipeline](https://github.com/connor-lab/ncov2019-artic-nf), mapping each sequence to the closest GenBank accession number hit in the original metagenomic run of IDseq, using a cutoff of 10 reads per nucleotide site to make a consensus call (sites with <10 reads were called as "N").

Using this method, we generated 63 full genome sequences of DENV-1 and 120 full genome sequences of DENV-2. Sequences were verified visually by comparing raw reads to consensus genomes in Geneious, as well as by checking their congruence with conensus genomes using a similar approach to the Nextflow pipeline implemented in CZID. All contributed genomes were >10000 bps in length and a maximum of 90 Ns (corresponding to <1% of the dengue virus genome). We additionally contributed 9 DENV-4 genomes to GenBank, though data were too sparse to analyze in earnest.

---

### Alignment and model selection

We include the combined files of both NIH and background SE Asian sequences from GenBank in the subfolder "alignment" as files "GenBankDENV1_prename.fasta" and "GenBankDENV2_prename.fasta". Sequence headers were renamed in the script "rename-fasta.R" to be compatible with BEAST, such that the date of each sequence could be easily discerned after the underscore ("_") in each title. These files ("allDENV1_beast.fasta" and "allDENV2_beast.fasta") were then sent to [MAFFT](https://mafft.cbrc.jp/alignment/server/) for alignment, prior to tree-building in BEAST. The resulting alignment files are contained in the alignment folder ("DENV1aligned.fasta" and "DENV2aligned.fasta").

After alignment, we evaluated the optimal nucleotide substitution model for each set of genomes, using [ModelTest-NG](https://github.com/ddarriba/modeltest). Output files from both DENV-1 and DENV-2 runs of ModelTest-NG are available for viewing in the subfolder, "modeltest-ng-output" within the "alignment" folder. Both found the strongest support (AIC and AICc) for GTR+I+G4 nucleotide substitution models. 

---

### Building a phylogenetic tree using BEAST2

The BEAST and BEAST2 communities maintain a number of extremely helpful tutorials that you should practice to get BEAST up and running on your home computer. See [here](https://taming-the-beast.org/tutorials/). 

After model selection, we next built a Bayesian phylogenetic tree for each set of dengue genomes in BEAST2, using the GTR+I+G4 nucleotide substitution model recommended by ModelTest-NG. The first step in this process requires generation of a .xml file for input into BEAST, using the program BEAUti. (Note: The GTR+I+G4 substitution model is easily specified in BEAUti; however, for less common forms, see this great blog post [here](https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/). Additionally, it is possible to generate .xml files outside of BEAUti, for example in a pipeline scenario, though this approach was not needed for this phylogeny).

To prepare the .xml file, we used the following parameters in tab inputs at the top of the screen in BEAUti:
 - Tip Dates: We used the date of sample collection as the "Tip Date." For any sample from GenBank which only listed year of collection, we set the tip date to July 31 of the year of collection. Each alignment was uploaded to BEAUti with sequence names arranged so as to easily decipher the date.
 - Site Model: Following output from ModelTest-NG, we selected a "Gamma Site Model" with Gamma category 4 and an estimated 0.001 for the proportion invariant. The exact parameter entries can be seen here:
 
![](BEAST-tree/guide-pics/BEAST-site-model.png)

- **Clock Model**: After [Salje et al. 2017](https://science.sciencemag.org/content/355/6331/1302.abstract), we adopted a strict molecular clock. Similar coalescent times were observed with relaxed clock models as well. We set the clock rate at 7.9 x 10^-4 substitutions/site/year, the median reported value in [Sall et al. 2020](https://journals.asm.org/doi/full/10.1128/JVI.01738-09).

![](BEAST-tree/guide-pics/BEAST-clock-model.png)

- **Priors**: We used a Coalescent Bayesian Skyline model. Nucleotide substitution rates were left at default, assuming a Gamma distribution, Alpha=0.05, Beta=10, initial value of 1, spanning a range of 0 to infinity. The following table specifies all other priors (Offset = 0 and Miniordimension= 1 for all priors):

| Prior  | Distribution  | Parameter 1  | Parameter 2  | Mean In Real Space |  Initial Value | Lower | Upper | Dimension | Estimate
|---|---|---|---|---|---|---|---|---|---|
| bPopSizes  |  -- | --  |  -- |--   | 380  |0 |  -- | 1  | check |
| clockRate  |  logNormal | 0.001  |  1.25 |check   | 0  |-- |  -- | 1  | check |
| freqParameter  | Uniform   | 0  | 1  | --  | 0.25  | 0  | 1  |  4  |check |
| gammaShape  | exponential  | 1  | --  | --  | 1  |  -- |  -- |  1 | check  |
| proportionInvariant  | Uniform  |  0 |  1 | --  |  0.001 | 0  |  1 | 1  | check |

The prior entry page can be seen here:

![](BEAST-tree/guide-pics/BEAST-priors.png)

- **MCMC**: We used an MCMC chain length of 100,000,000 iterations and set tracelog and treelog every 10,000 iterations. All other MCMC metrics were left at default. 

- The original DENV1.xml and DENV2.xml files can be found in the "beast-denv" subfolder of  "alignment" folder. This were then run using BEAST2 on the UChicago and UC Berkeley computing cluster to produce the files found in the "BEAST-output" folder.

---

### Visualizing Bayesian timetree

The initial 10% of MCMC iterations were removed as burn-in. Parameter convergence was assessed visually using Tracer v1.6. Once verified, we used TreeAnnotator to average across the BEAST tree output, then visualized the resulting tree in R. Scripts for preparation of the timetree witnessed in our paper are available in folder "Figures" and subfolders "Fig1" and "Fig3". 
