import os

samples = ['Ctrl1', 'Ctrl2', 'Ctrl3', 'NASH1', 'NASH2', 'NASH3']
groups = ["Ctrl", "NASH"]
marks = ['H3K27ac', 'H3K27me3', 'H3K4me3']

bam_dir = "chip_bam"
out_dir = "5.4"

rule all:
    input:
        expand(os.path.join(out_dir, "{sample}.sorted.bdg"), sample=samples),
        expand(os.path.join(out_dir, "{sample}.sorted.bw"), sample=samples),
        expand(os.path.join(out_dir, "{mark}.{group}.matrix.gz"), mark=marks, group=groups),

rule sort_bdg:
    input:
        chrsize = "GCF_012559485.2.chrom.sizes",
        bdg = os.path.join(bam_dir, "{sample}_MergedCG_5x.sort.bedgraph"),
    output:
        os.path.join(out_dir, "{sample}.sorted.bdg"),
    shell:
        '''
        awk -F '\\t' '(NR==FNR){{a[$1];next}}($1 in a){{OFS="\\t";print $0}}' {input.chrsize} {input.bdg} |sort -k1,1 -k2,2n - > {output}
        '''

rule bdg_to_bw:
    input:
        chrsize = "GCF_012559485.2.chrom.sizes",
        bdg = os.path.join(out_dir, "{sample}.sorted.bdg"),
    output:
        os.path.join(out_dir, "{sample}.sorted.bw"),
    shell:
        '''
        bedGraphToBigWig {input.bdg} {input.chrsize} {output}
        '''

rule computeMatrix:
    input:
        bw = os.path.join(out_dir, "{mark}.{group}.bw"),
        gene = "GCF_012559485.2_MFA1912RKSv2.ncbiRefSeq.gene.bed"
    output:
        matrix = os.path.join(out_dir, "{mark}.{group}.matrix.gz"),
    params:
        computeMatrix = "computeMatrix", 
    shell:
        """
        {params.computeMatrix} scale-regions\
            --scoreFileName {input.bw} \
            --regionsFileName {input.gene} \
            --samplesLabel {wildcards.mark} \
            -o {output.matrix}\
            -p {resources.tasks}\
            -b 2000\
            -a 0\
            -m 5000\
            --missingDataAsZero
        """
