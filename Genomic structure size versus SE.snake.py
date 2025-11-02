bed_dir = "results"
prefix = [f"random_{i}" for i in range(100)] + ["SE"]
work_dir = "results"
resolutions = [5000, 10000]
resolutions_tad = [10000, 25000]

rule all:
    input:
        expand(os.path.join(work_dir, "{prefix}.{resolution}.bedpe"), prefix=prefix, resolution=resolutions),
        expand(os.path.join(work_dir, "{prefix}.{resolution}.bed"), prefix=prefix, resolution=resolutions_tad),

rule loop:
    input:
        loop = "NASH.{resolution}.hicpeaks.txt", 
        bed = os.path.join(bed_dir, "{prefix}.bed"),
    output:
        os.path.join(work_dir, "{prefix}.{resolution}.bedpe"), 
    shell:
        """
        bedtools pairtobed -a {input.loop} -b {input.bed} > {output}
        """

rule tad:
    input:
        tad = "NASH.{resolution}.domain.bed", 
        bed = os.path.join(bed_dir, "{prefix}.bed"),
    output:
        os.path.join(work_dir, "{prefix}.{resolution}.bed"), 
    shell:
        """
        bedtools intersect -a {input.tad} -b {input.bed} -wa|uniq > {output}
        """
