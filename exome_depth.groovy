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

    
    produce("exome_depth.${chr}.tsv") {

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
            chr.covered = read.bed(pipe("grep '$chr[^0-9]' $target_bed"))

            # ExomeDepth wants the columns named differently
            chr.covered = data.frame(
                chromosome=chr.covered\$chr, 
                start=chr.covered$start, 
                end=chr.covered$end,  
                name=paste(chr.covered\$chr,chr.covered$start,chr.covered\$end,sep="-")
            )

            # BAM files are the primary input - convert to R vector
            bam.files = c('${inputs.bam.join("','")}') 

            all.samples = sapply(bam.files, function(file.name) {
                # Scan the bam header and parse out the sample name from the first read group
                read.group.info = strsplit(scanBamHeader(file.name)[[1]]$text[["@RG"]],":")
                names(read.group.info) = sapply(read.group.info, function(field) field[[1]]) 
                return(read.group.info$SM[[2]])
            })

            print(sprintf("Processing %d samples",length(all.samples)))

            # Finally we can call ExomeDepth
            bam.counts <- getBamCounts(bed.frame = chr.covered,
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
                reference = select.reference.set(
                                         test.counts = bam.counts[,test.sample][[1]],
                                         reference.counts = as.matrix(as.data.frame(bam.counts[,reference.samples])[,-1:-6]),
                                         bin.length = chr.covered$end - chr.covered$start
                                        )

                # Get counts just for the reference set
                reference.counts = apply(bam.counts.df[,reference\$reference.choice,drop=F],1,sum)

                print(sprintf("Creating ExomeDepth object ..."))
                ed = new("ExomeDepth",
                              test = test.sample.counts,
                              reference = reference.counts,
                              formula = "cbind(test, reference) ~ 1")


                print(sprintf("Calling CNVs ..."))
                found.cnvs = CallCNVs(x = ed,
                                        transition.probability = $transition_probability,
                                        chromosome = chr.covered$chromosome,
                                        start = chr.covered$start,
                                        end = chr.covered$end,
                                        name = chr.covered$name)


                results = found.cnvs@CNV.calls
                results$sample = rep(test.sample, nrow(results))

                print(sprintf("Writing %d results to $output.tsv ...", nrow(results)))
                write.table(file="$output.tsv", 
                            x=results,
                            row.names=F,
                            append=T)
            }

            print("Finished")

        """},'exome_depth')
    }
}

@produce("exome_depth.cnvs.tsv")
merge_ed = {

    doc "Merges results from multiple ExomeDepth analyses together"

    exec """
        cat $inputs.tsv | grep -v '^"sample"' | awk '{ if(NR==1 || \$1 != "\\"start.p\\"") print \$0 }' > $output.tsv
    """
}
