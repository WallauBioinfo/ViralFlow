process coveragePlot {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

    input:

      tuple val(sample_id), path(bam), path(bai)
    
    output:
        path("*coveragePlot*"), optional: true
        path("coveragePlot_result.txt"), optional: true, emit: result, hidden: true
    script:
    
    //bam = bam_file[0].toString()
    depth = params.depth
    html = "${sample_id}_coveragePlot.html"
    png = "${sample_id}_coveragePlot.png"
    svg = "${sample_id}_coveragePlot.svg"
    
    """
    #!/usr/bin/env python
    # ----- import libraries -------------------------------------------------
    import pysam
    import subprocess

    #Loading the BAM file
    bam_file = pysam.AlignmentFile('${bam}', 'rb')

    #Checking the number of mapped reads
    n_mapped_reads = bam_file.count()
    
    if n_mapped_reads > 0:
      
      references = bam_file.references
      if len(references) == 0:
        raise ValueError("No references found in BAM header.")
        
      genomecode = references[0]
    
      subprocess.run([f"bamdash -b ${bam} -r {genomecode} -c ${depth} -e svg"], shell=True)
      subprocess.run([f"mv {genomecode}_plot.svg ${svg}"], shell=True)

      subprocess.run([f"bamdash -b ${bam} -r {genomecode} -c ${depth} -e png"], shell=True)
      subprocess.run([f"mv {genomecode}_plot.png ${png}"], shell=True)

      subprocess.run([f"bamdash -b ${bam} -r {genomecode} -c ${depth}"], shell=True)
      subprocess.run([f"mv {genomecode}_plot.html ${html}"], shell=True)
   
    else:
      result = "No mapped reads were found in the sorted BAM file for sample ${sample_id}. The coverage plot will not be generated for it."
      with open('coveragePlot_result.txt', 'w') as f:
        f.write(result)
      """
}

process snpPlot {
    publishDir "${params.outDir}/${sample_id}_results/", mode: "copy"

    input:

      tuple val(sample_id), path("${sample_id}.depth${params.depth}.fa.algn")

      
    output:
        path("*snpPlot*")

    script:
    
    plot_name = "${sample_id}_snpPlot"
    

    """

    snipit ${sample_id}.depth${params.depth}.fa.algn -o ${plot_name} --solid-background 
    snipit ${sample_id}.depth${params.depth}.fa.algn -o ${plot_name} -f svg --solid-background

    """
}
