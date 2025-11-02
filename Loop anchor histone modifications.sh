computeMatrix reference-point\
    --referencePoint center\
    --scoreFileName H3K4me3.Ctrl.bw  H3K4me3.NASH.bw H3K27ac.Ctrl.bw  H3K27ac.NASH.bw H3K27me3.Ctrl.bw  H3K27me3.NASH.bw \
    --regionsFileName NASH.uniqe.regions.bed \
    -o heat.matrix.gz \
    -p 30 \
    -b 3000\
    -a 3000\
    --missingDataAsZero

plotHeatmap -m heat.matrix.gz \
     -out heat.png \
     --colorMap Reds Reds Oranges Oranges Blues Blues\
     --whatToShow 'heatmap and colorbar' \
     --kmeans 4 \
    --outFileSortedRegions region.bed \
    --zMin 0 0 0 0 0 0 \
    --zMax 8 8 5 5 3 3\
     --dpi 100 
