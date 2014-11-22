exome_depth_pipeline
====================

A Bpipe pipeline for easily running ExomeDepth.

Overview
========

ExomeDepth is a fantastic tool for identifying CNVs in exome or targeted sequencing data.

This small pipeline just makes it a bit easier to run by implementing a set of defaults
that should work for most situations. With this pipeline you just need to supply
the following:

  * a set of BAM files - one sample per BAM file 
  * a BED file of the target region of your exome or targeted sequencing data
  * the human genome reference FASTA to which the BAM files were aligned

Install
=======

1.  To run the pipeline you need to install Bpipe:

http://bpipe.org


2. Once you have done that, you can check out this pipeline anywhere:

git clone git@github.com:ssadedin/exome_depth_pipeline.git

3. You also need to install ExomeDepth - see the instructions for that which
come with ExomeDepth.

Running
=======

An example of running the pipeline is like this:


    bpipe run -p target_bed=target_regions.bed \
              -p ref=gatk.ucsc.hg19.fasta \
              exome_depth_pipeline.groovy *.bam

Of course, you need to replace the target regions, reference fasta and BAM files with your own.

On a small target region this takes about 2 minutes to run on a moderately powerful server. If you do not have 
enough RAM to run all the chromosomes in parallel, you can limit throughput with -n <num threads> as an argument
to Bpipe.
