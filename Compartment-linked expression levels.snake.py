compartment_resolutions = [
    100000, 50000, 25000
]

rule all:
    input:
        expand("{resolution}.ab_gene_stat.svg", resolution=compartment_resolutions),

rule ab_gene_stat:
    input:
        compartment = "NASH.{resolution}.cis.vecs.tsv",
        name = "gene_transcript.txt",
        rna = "gene_TPM.xls",
        gene_bed = "MFA1912RKSv2.ncbiRefSeq.transcript.bed",
    output:
        comparment_gene_file = "{resolution}.comparment_gene.txt", 
        tb = "{resolution}.ab_gene_stat.txt",
        fig = "{resolution}.ab_gene_stat.png",
    params:
        Rscript = "Rscript",
        python = "python3.9",
    shell: 
        """
        bedtools intersect -a {input.compartment} -b {input.gene_bed} -wa -wb|awk -F"\\t" '{{OFS="\\t";print $5,$11, $1"-"$2"-"$3}}'|uniq > {output.comparment_gene_file} 
        {params.python} work.py {output.comparment_gene_file}  {output.tb}
        {params.Rscript} box.R {output.tb} {output.fig}
        """

rule plot:
    input:
        "{resolution}.ab_gene_stat.txt",
    output:
        "{resolution}.ab_gene_stat.svg",
    shell:
        """
        Rscript 8.2-2.box.R {input} {output}
        """
