# DOGE - Discove Organelle GEnome Variants

**A reference-free, CLI / GUI and multi-platform compatible and portable pipeline for high-quality organelle SNP calling**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-GPL-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

![DOGE icon](https://github.com/Gipsy-The-Sheller/DOGE/blob/main/icon.svg){:height="50px"}

DOGE is a cross-platform program aiming to handle SNP calling of organelle genomes and other genomes with complex structural variants (e.g. reversions and rearrangements), using ProgressiveMauve as its core alignment tool. It uses locally collinear blocks (LCBs) as basic SNP calling units and inherently supports auto partitioning based on LCBs. Moreover, DOGE can be run on Linux, Mac OS X and **Windows**, as its dependencies, ProgressiveMauve, MAFFT and TrimAl, are all available on these systems. Apart from command line mode, DOGE has a built-in GUI based on PyQt5. It can not only be started up in command line but also be embedded to other PyQt5-based desktop bioinformatic platforms as a plugin.

[TOC]

## Why DOGE?

Unlike reference-dependant SNP calling pipelines (e.g. Snippy), DOGE employs a reference-free approach, which enables itself to take full use of sequences and call SNPs for phylogenies without reference-caused biases. Also, DOGE identifies collinearity of whole sequences, taking full use of phylogenetic signals from homologous none-coding regions, not only gene regions.

below is a humorous comic illustrating advantages of DOGE.

![Why DOGE](https://github.com/Gipsy-The-Sheller/DOGE/blob/main/artworks/WhyDOGE.png)

## How it works?

A standard process of DOGE pipeline includes 3 necessary steps and 3 optional steps. Below is a list illustrating how it work.

1. **Whole-genome alignment with Mauve** This is the main step of DOGE. After the alignment, all LCBs are splitted from ProgressiveMauve's output xmfa file.
2. **Refining alignment with MAFFT (optional)** This step uses MAFFT to re-align each LCB. Note that we noticed **this step may introduce errors to LCBs with many gap regions** but it is kept in the pipeline and by default, skipped in the analysis. However, you can use `-r True` or `--refine_alignment True` to  use this step.
3. **Filtering out low-coverage LCBs (optional)** This step discards LCBs shared in genomes which count less than a threashold. By default, the threashold is 3, a minimum number to produce parsimony-informative sites for phylogenies. You can also customize the threashold.
4. **Trimming LCBs (optional)** This step uses TrimAl (with --automated1 parameter to fit in different LCB circumstances) to trim gap regions of LCBs.
5. **SNP calling** This step calls SNPs from processed LCB regions and forms a fasta supermatrix.
6. **Partitioning** This step allocates each SNP site with a partition based on its belonging LCB and generate:
   - a NEXUS file containing the SNP matrix and the partition scheme;
   - a NEXUS file with only partition charsets.

below is an illustration of DOGE process.

![How DOGE Works](https://github.com/Gipsy-The-Sheller/DOGE/blob/main/artworks/HowDOGEworks.png)

## Quick start

You can get DOGE by 2 means: *download the source code* or *get a released package* (only available for Windows users now)

### Configure runtime environment for source code

As a python script, DOGE has the following dependencies: `BioPython` and `PyQt5`. You may install these by pip:

```bash
pip install biopython PyQt5
```

The default paths to ProgressiveMauve, MAFFT and TrimAl are respectively:

```
./mauve/progressiveMauve.exe
./mafft/mafft.bat
./trimal/bin/trimal.exe
```

You can also install them to other directories and customize their paths by using the arguments `--mauve`, `--mafft` and `--trimal`.

### Run DOGE in CLI mode

DOGE CLI requires 2 arguments: The `--input` directory containing all genomes' FASTA files, and the `--output` directory for its results. For advanced usage of DOGE, you can get a detailed argument list by typing `--help`:

```
DOGE: Discover Organelle GEnome Variants

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input directory with genome FASTA files [required]
  -o OUTPUT, --output OUTPUT
                        Output directory for results [required]
  --mauve MAUVE         Path to Mauve executable
  --mafft MAFFT         Path to MAFFT executable
  --trimal TRIMAL       Path to TrimAl executable
  -l LOWEST_COVERAGE, --lowest_coverage LOWEST_COVERAGE
                        Minimum number of sequences in an LCB file
  -r REFINE_ALIGNMENT, --refine_alignment REFINE_ALIGNMENT
                        True / False: Refine MAUVE alignment of LCBs with MAFFT (Warning: perhaps harmful  
                        for distantly related taxa!)
  -f FILTER, --filter FILTER
                        True / False: Filter out low-coverage LCBs
  -t TRIM, --trim TRIM  True / False: Trim gap regions with TrimAl before SNP calling
  -w, --window          Launch GUI with pre-filled parameters
```

### Run DOGE in GUI mode

DOGE has a PyQt5-based GUI that can be started from command line:

```
python doge.py -w
```

or

```
python doge.py --window
```

After this a GUI window of DOGE will appear. If you have already typed any other arguments, then they'll be automatically filled in the GUI and will be highlighted green.

![DOGE GUI](https://github.com/Gipsy-The-Sheller/DOGE/blob/main/artworks/DOGE-GUI.png)

### Integrate DOGE to other platforms

The core GUI of DOGE is placed into a `QWidget` class, which allows it to be integrated into other platforms developed with PyQt5. Any scripts can get a copy of DOGE class by `DOGE.Plugin.run()`. DOGE has been integrated into `YRTools`, a new plugin-based, and native-Chinese supported bioinformatic platform developed by us, which will come out soon:

![DOGE in YRTools](https://github.com/Gipsy-The-Sheller/DOGE/blob/main/artworks/DOGE-YRTools.png)

## Citations

If you use DOGE, please cite its Github repository or its poster, and used dependencies:

1. https://github.com/Gipsy-The-Sheller/DOGE or
2. Xu, Z.-J. DOGE: Discover Organelle GEnome Variants. DOI: http://dx.doi.org/10.13140/RG.2.2.24360.64001
3. Darling AE, Mau B, Perna NT. progressiveMauve: multiple genome alignment with gene gain, loss and rearrangement. PLoS One. 2010 Jun 25;5(6):e11147. doi: 10.1371/journal.pone.0011147.

If you use MAFFT in the pipeline, please cite:
- Katoh, K., Standley, D.M., 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol. Biol. Evol. 30, 772-780.

If you use TrimAl in the pipeline, please cite:
- Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics. 25: 1972-1973. doi: 10.1093/bioinformatics/btp348.

Hopefully you may find DOGE useful for your work.