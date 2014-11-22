// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// ExomeDepth Pipeline - The core pipeline stage is in here.
//
//////////////////////////////////////////////////////////////////
exome_depth = {

    doc title : "Exome Depth Per-Chromosome Analysis",
        desc : """
               Runs ExomeDepth on a single chromosome for a set of BAM files.
               <p>
               NOTE: each BAM file is assumed to contain a single sample. Multisample BAM
                     files are not supported.
               """

    requires target_bed : "BED file containing target regions to analyse (eg: exome design covered regions)",
             chr : "Chromosome to process",
             ref : "Genome reference in FASTA format, indexed using samtools faindex"

    var transition_probability : "0.0001"

    msg "Using Exome Depth transition probability = $transition_probability"

    

    R({"""

        library(ExomeDepth)
        library(Rsamtools)

        read.bed = function(f) {
          read.table(f, col.names=c("chr","start","end","id"), fill=1)
        }

        # Reference sequence
        hg19.fasta = "$ref"

        # Read the target / covered region
        print("Reading target regions for $chr from $target_bed")
        dsd.covered = read.bed(pipe("grep '$chr[^0-9]' $target_bed"))

        # ExomeDepth wants the columns named differently
        dsd.covered = data.frame(
            chromosome=dsd.covered\$chr, 
            start=dsd.covered$start, 
            end=dsd.covered$end,  
            name=paste(dsd.covered\$chr,dsd.covered$start,dsd.covered\$end,sep="-")
        )

        # BAM files are the primary input - convert to R vector
        bam.files = c('${inputs.bam.join("','")}') 

        all.samples = sapply(bam.files, function(bam.file) {
            # Scan the bam header and parse out the sample name from the first read group
            read.group.info = strsplit(scanBamHeader(file.name)[[1]]$text[["@RG"]],":")
            names(read.group.info) = sapply(read.group.info, function(field) field[[1]]) 
            return(read.group.info$SM[[2]])
        })

        print(sprintf("Processing %d samples",length(all.samples)))

        # Finally we can call ExomeDepth
        bam.counts <- getBamCounts(bed.frame = dsd.covered,
                                  bam.files = bam.files,
                                  include.chr = F,
                                  referenceFasta = hg19.fasta)

        print("Successfully counted reads in BAM files")

        # Note: at this point bam.counts has column names reflecting the 
        # file names => convert to actual sample names which is more convenient
        colnames(bam.counts) = c("GC", all.samples)

        # Problem: sample names starting with numbers get mangled. So convert them back, 
        # but ignore the first column which is actually the GC percentage
        all.samples = colnames(bam.counts)[-1]

        for(test.sample in all.samples) {

            print(sprintf("Processing sample %s", test.sample))

            reference.samples = all.samples[-match(test.sample, all.samples)]

            bam.counts.df = as.data.frame(bam.counts[,reference.samples])[,-1:-6]

            test.sample.counts = bam.counts[,test.sample][[1]]

            print(sprintf("Selecting reference set for %s ...", test.sample ))
            dsd.reference = select.reference.set(
                                     test.counts = bam.counts[,test.sample][[1]],
                                     reference.counts = as.matrix(as.data.frame(bam.counts[,reference.samples])[,-1:-6]),
                                     bin.length = dsd.covered$end - dsd.covered$start
                                    )

            # Get counts just for the reference set
            dsd.reference.counts = apply(bam.counts.df[,dsd.reference\$reference.choice,drop=F],1,sum)

            print(sprintf("Creating ExomeDepth object ..."))
            ed = new("ExomeDepth",
                          test = test.sample.counts,
                          reference = dsd.reference.counts,
                          formula = "cbind(test, reference) ~ 1")


            print(sprintf("Calling CNVs ..."))
            found.cnvs = CallCNVs(x = ed,
                                    transition.probability = $transition_probability,
                                    chromosome = dsd.covered$chromosome,
                                    start = dsd.covered$start,
                                    end = dsd.covered$end,
                                    name = dsd.covered$name)


            results = found.cnvs@CNV.calls
            results$sample = rep(test.sample, nrow(results))

            print(sprintf("Writing %d results to $output.tsv ...", nrow(results))
            write.table(file="$output.tsv", 
                        x=results,
                        row.names=F,
                        append=T)
        }

        print("Finished")

    """},'exome_depth')
}

@produce("exome_depth.cnvs.tsv")
merge_ed = {

    doc "Merges results from multiple ExomeDepth analyses together"

    exec """
        cat $inputs.exome_depth.tsv | grep -v '^"sample"' | awk '{ if(NR==1 || \$1 != "\\"start.p\\"") print \$0 }' > $output.tsv
    """
}
