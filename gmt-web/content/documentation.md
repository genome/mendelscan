## Documentation ##
----

MendelScan is a command-line program with multiple subcommands (e.g. score, rhro, and sibd). Each subcommand has a unique set of inputs and outputs. For the list of available subcommands, enter:

<pre class="terminal">
java -jar MendelScan.jar --help
</pre>

## Available Subcommands ##
----

These subcommands are currently supported:

<pre class="terminal">
java -jar MendelScan.jar score &#35; Prioritize a VCF
java -jar MendelScan.jar rhro  &#35; Perform RHRO analysis
java -jar MendelScan.jar sibd  &#35; Perform SIBD analysis
</pre>

For detailed usage information, enter the subcommand followed by -h or --help, e.g.:

<pre class="terminal">
java -jar MendelScan.jar score -h
</pre>

For those familiar with Java, the auto-generated [Javadoc][] documentation may be useful as well.

[Javadoc]: http://genome.github.io/mendelscan/

### score: Variant Scoring and Prioritization ###

The `score` command of MendelScan takes 4 inputs:

1. A pedigree file in [PED format][PED] that indicates the name, gender, and affectation status of the samples in the VCF. Samples in the VCF but not in the PED file will be treated as affected females.
2. A VCF file that has been annotated with dbSNP information (a task that can be completed with the current dbSNP VCF file and the [joinx][] utility).
3. Variant annotation information in Variant Effect Predictor (VEP) format, ideally with canonical, hgnc, Polyphen, SIFT, and Condel options.
4. Gene expression for the tissue(s) of interest. This should be a one-column text file with HUGO symbols ordered according to their expression level (highest to lowest). This is optional but highly recommended; many gene expression datasets are freely available.

MendelScan calculates four individual scores (segregation, population, annotation, and expression) for each variant. Each score is a value between 0 and 1 reflecting the likelihood that the variant could be disease-causing. The default settings will prioritize novel/rare protein-altering variants in highly expressed genes that segregate in autsomal-dominant fashion. Many of the scoring parameters can be adjusted to suit different kinds of studies. An overall score, taken as the product of the four scores, reflects the relative priority (higher = more likely to cause disease) of variants based on these criteria.

The output file contains each variant along with the overall and individual scores, as well as annotation, population, expression, and segregation data that were used to compute them. A VCF output option is also available; it places an similar but abbreviated information in the INFO field.

### rhro: Rare Heterozygote Rule Out ###

The `rhro` subcommand of MendelScan takes three inputs:

1. A pedigree file in [PED format][PED] that indicates the name, gender, and affectation status of the samples in the VCF. Samples in the VCF but not in the PED file will be treated as affected females.
2. A VCF file that has been annotated with dbSNP information (a task that can be completed with the current dbSNP VCF file and the [joinx][] utility).
3. A BED file of chromosome centromere coordinates (optional but recommended).

The RHRO method identifies candidate regions consistent with autosomal dominant inheritance based on the idea that a disease-causing haplotype will manifest regions of rare heterozygous variants shared by all affecteds, and an absence of homozygous differences between affected pairs (which would indicate that a pair had no haplotype in common).

There are two output files from this command. One contains all informative variants (rare heterozygotes shared by affecteds, or variant positions with homozygous differences between affected pairs). The second output is a window of RHRO regions that are consistent with autosomal dominant inheritance given the inputs and assumptions described here.

### sibd: Shared Identity-by-Descent ###

The `sibd` subcommand of MendelScan uses BEAGLE FastIBD results to identify regions of maximum identity-by-descent (IBD) among affected pairs. It requires the user to run BEAGLE FastIBD on the sequencing data (which requires conversion of the VCF to BEAGLE format and a "markers" file). This should be done on a per-chromosome basis. Then, the following files should be provided as inputs to MendelScan for each chromosome:

1. A pedigree file in [PED format][PED] that indicates the name, gender, and affectation status of the samples in the VCF. Samples in the VCF but not in the PED file will be treated as affected females.
2. The BEAGLE markers file for the chromosome at hand, which typically includes four columns: physical position (chrom:position), map position (morgans), allele1, and allele2.
3. The BEAGLE FastIBD output file (*.fibd) for the chromosome in uncompressed format. It should have five columns: sample1, sample2, index1, index2, and score. The index fields correspond to the markers file; MendelScan will convert these to genomic coordinates and print them to the output file.

MendelScan breaks the chromosome into windows of a user-specified resolution (default: 100,000 bp) and, for each window, determines the number of affected pairs that shared an IBD segment in that window. Typically, windows in which >90% of possible affected pairs were IBD suggests a candidate haplotype. All windows are output to a second output file (if specified) or STDOUT.

### Example

Included in the repository is an example data set using [1000 Genomes][] data.  You extract that data and run the following example:

<pre class="terminal">
$ tar -zxvf example_data.tar.gz
$ cd example_data
$ java -jar MendelScan.jar score variants.vcf \
    --vep-file annotation.vep \
    --ped-file family.ped \
    --gene-file gene-expression.txt \
    --output-file mendelscan.tsv \
    --output-vcf mendelscan.vcf
Reading input from variants.vcf
Loading sample information from family.ped...
1 males, 2 cases, 1 controls
Loading gene expression information from gene-expression.txt...
Expression rank loaded for 38545 genes
Loading VEP from annotation.vep...
11181 variants had VEP annotation
Scoring variants under dominant disease model
3 samples in VCF (2 affected, 1 unaffected, 1 male)
11181 variants in VCF file
11181 matched with VEP annotation
12846   variants_common
337 variants_known
18  variants_mutation
97  variants_novel
1359    variants_rare
466 variants_uncommon
</pre>

[1000 Genomes]: http://www.1000genomes.org
