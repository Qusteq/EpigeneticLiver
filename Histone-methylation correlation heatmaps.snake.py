marks = ['H3K27ac', 'H3K4me3', 'H3K27me3']
bw_dir = "../5.4"
work_dir = "./6.4-5"
kinds = ["Low", "Intermediate", "High"]

rule all:
	input:
		expand(os.path.join(work_dir, "{mark}.methylation.heatmap.svg"), mark=marks), 

rule matrix_by_methy:
    input:
        bw = os.path.join(bw_dir, "{mark}.NASH.bw"),
        bed = expand(os.path.join(work_dir, "{kind}-methylation.bed"), kind=kinds),
        mat = os.path.join(work_dir, "{mark}.methylation.mat.gz"),
    output:
        png = os.path.join(work_dir, "{mark}.methylation.heatmap.svg"),
    shell:
        """
		plotHeatmap -m {input.mat} -out {output.png} \
			--sortRegions descend --sortUsing median \
			--regionsLabel 'Low' 'Intermediate' 'High' \
			-x '' --dpi 500 
        """
