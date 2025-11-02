forms = ["stable", "expand", "shrink", "shift"]
resolutions = [10000, 25000]
out_dir = "results/"
rule all:
    input:
        expand(os.path.join(out_dir, "{form}.{resolution}.genes"), form=forms, resolution=resolutions),

rule work:
    input:
        gene_bed = "MFA1912RKSv2.ncbiRefSeq.gene.bed", 
        tad_domain = "{form}.{resolution}.domain.bed",
    output:
        gene_file = os.path.join(out_dir, "{form}.{resolution}.genes"),
    shell:
        """
        bedtools intersect -a {input.gene_bed} -b {input.tad_domain} -wa|awk '{{print $4}}'|sort|uniq > {output.gene_file}
        """    
