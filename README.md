# ![nf-funPHEhr](docs/images/nf-core-funphehr_logo_light.png#gh-light-mode-only) ![funPHEhr](docs/images/nf-core-funphehr_logo_dark.png#gh-dark-mode-only)


## Introduction

**nf-funphehr** is a bioinformatics pipeline that can be used to analyse Nanopore sequencing data obtained from fungal isolates. It takes a samplesheet and fastq files as input, performs QC, trimming, assembly, assembly QC and annotation(in dev). This is the public version of the developing pipeline being trialled at SA Pathology. 

![nf-funphehr metro map](docs/images/nf-funphehr_metro_map.png)
1. Trim adapters ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Trim low quality and short reads (['chopper`](https://github.com/wdecoster/chopper))
3. Present QC of reads post trimming (['Nanoplot'](https://github.com/wdecoster/NanoPlot))
4. Screen reads for contamination (['kraken2'](https://github.com/DerrickWood/kraken2))
5. Denovo Assembly
   - flye (['flye'](https://github.com/fenderglass/Flye))
   OR 
   miniasm -> minimap2 -> racon 
   (['miniasm'](https://github.com/lh3/miniasm))
   (['minimap2'](https://github.com/lh3/minimap2))
   (['racon'](https://github.com/lbcb-sci/racon))

6. Polishing reads (optional) (['medaka'], (https://github.com/nanoporetech/medaka))
7. extract of ITS1-5.8S-ITS2 region (['ITSx'](https://microbiology.se/software/itsx/))
8. Assembly assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
9. Assembly quality (['BUSCO'](https://busco.ezlab.org/))
4. Present metrics from run ([`MultiQC`](http://multiqc.info/))

## Usage
:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

** Please not nextflow.config file may need to be have paths updated for databases for kraken2, busco and for annotation steps

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
ID,LongFastq,GenomeSize, species
231862455,./data/S1_long_fastq.gz,14.0m,"Candida albicans"
231495562,./data/S1_long_fastq.gz,26.0m,"Candida parapsilosis"

```

Each row represents a fastq file (single-end) of long reads.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run https://github.com/SAP-PHE-Bioinformatics/FunPHEhr/ \
   -profile <docker/apptainer/singularity/...> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::


## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/funphehr/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/funphehr/output).

## Credits

nf-funphehr was originally written by jacob.

Thank you to all and acknowledgement to all the authors of tools used. 
* __[chopper](https://github.com/wdecoster/chopper)__  
Rust implementation of NanoFilt+NanoLyse, both originally written in Python. This tool, intended for long read sequencing such as PacBio or ONT, filters and trims a fastq file.
_Wouter De Coster, Rosa Rademakers, [NanoPack2: population-scale evaluation of long-read sequencing data](https://doi.org/10.1093/bioinformatics/btad311, Bioinformatics, (2023)_  

* __[Flye](https://github.com/fenderglass/Flye)__  
De novo assembler for single molecule sequencing reads using repeat graphs  
_Kolmogorov, M, Yuan, J, Lin, Y, Pevzner, P, [Assembly of Long Error-Prone Reads Using Repeat Graphs](https://doi.org/10.1038/s41587-019-0072-8), Nature Biotechnology, (2019)_  

* __[Minimap2](https://github.com/lh3/minimap2)__  
A versatile pairwise aligner for genomic and spliced nucleotide sequences  
_Li, H [Minimap2: pairwise alignment for nucleotide sequences.](https://doi.org/10.1093/bioinformatics/bty191) Bioinformatics, 34:3094-3100. (2018)_  

* __[Miniasm](https://github.com/lh3/miniasm)__  
Ultrafast de novo assembly for long noisy reads (though having no consensus step)  
_Li, H [Miniasm: Ultrafast de novo assembly for long noisy reads](https://github.com/lh3/miniasm_  

* __[Medaka](https://github.com/nanoporetech/medaka)__  
Sequence correction provided by ONT Research  
_Li, H [Medaka: Sequence correction provided by ONT Research](https://github.com/nanoporetech/medaka)_  

* __[Porechop](https://github.com/rrwick/Porechop)__  
Adapter trimmer for Oxford Nanopore reads  
_Wick, RR, Judd, LM, Gorrie, CL, Holt, KE, [Completing bacterial genome assemblies with multiplex MinION sequencing.](https://doi.org/10.1099/mgen.0.000132) Microb Genom. 3(10):e000132 (2017)_  

* __[BUSCO](https://busco.ezlab.org/)__  
Adapter trimmer for Oxford Nanopore reads  
_Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, [C BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes](https://doi.org/10.1093/molbev/msab199) Molecular Biology and Evolution (2021)_  


* __[Racon](https://github.com/lbcb-sci/racon)__  
Ultrafast consensus module for raw de novo genome assembly of long uncorrected reads  
_Vaser, R, Sović, I, Nagarajan, N, Šikić, M, [Fast and accurate de novo genome assembly from long uncorrected reads.](http://dx.doi.org/10.1101/gr.214270.116) Genome Res. 27, 737–746 (2017)._  

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#funphehr` channel](https://nfcore.slack.com/channels/funphehr) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

This pipeline was built using nextflow and following nf-core template. 

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
