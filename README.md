# MendelScan

## Overview

Please see the official site for MendelScan, on [GMT][], for detailed information.

[GMT]: http://gmt.genome.wustl.edu/mendelscan/

## Releases

Please see the [releases][] page for official releases.

## Build

You can build the Java source by running `ant compile` or a JAR file by running `ant dist`.

[releases]: https://github.com/genome/mendelscan/releases

## Run

You can run the JAR file on the commandline like any Java JAR file, e.g. `java -jar ./MendelScan.jar`.

## Source Documentation

Javadoc is published on the GitHub [page][gh-pages].

[gh-pages]: http://genome.github.io/mendelscan/

## Example

Included in the repository is an example data set using [1000 Genomes][] data.

    tar -zxvf example_data.tar.gz
    cd example_data
    java -jar MendelScan.jar score variants.vcf \
        --vep-file annotation.vep --ped-file family.ped --gene-file gene-expression.txt \
        --output-file mendelscan.tsv --output-vcf mendelscan.vcf

[1000 Genomes]: http://www.1000genomes.org
