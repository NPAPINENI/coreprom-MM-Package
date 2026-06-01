## A command-line pipeline to generate DNABERT-CorePromMM sequence properties and shape features  from a labeled promoter TSV.
#### Developed by the Davuluri Lab, Stony Brook University.
### Given a TSV file with DNA sequences and labels, this pipeline:
1. Generates 6-mer tokenized sequences (sliding window = 1)
2. Splits data into train / test / dev (70 / 15 / 15, stratified)
3. Computes 10 sequence composition features (A/C/G/T fractions, CpG1/2/3, etc.)
4. Computes 14 DNA shape features via DNAshapeR (MGW, Roll, ProT, HelT, EP, Rise, Shift, Slide, Tilt, Buckle, Opening, Shear, Stagger, Stretch)
5. Writes 6 output TSV files (2 per split)

### Requirements
Python >= 3.8
Anaconda or Miniconda
R >= 4.0 installed and accessible via the R command

### Step 1 — Install the package
python -m pip install git+https://github.com/NPAPINENI/coreprom-MM-Package.git

### Step 2 — Install R dependencies
git clone https://github.com/NPAPINENI/coreprom-MM-Package.git
Rscript coreprom-MM-Package/install_r_deps.R

### Step 4 — Run
coreprom-props Promoter input.tsv ./output/  (input tsv file should have columns: DNA Sequence and label)
