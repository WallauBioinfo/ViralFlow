process coveragePlot {
    tag "${meta.id}"
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: "copy", pattern: "*_coveragePlot.*"

    input:

      tuple val(meta), path(bam), path(bai)
      // bamdash -c: the depth a position must exceed to count as recovered in
      // the plot's stats. The caller passes its mode's consensus threshold, so
      // the figure describes the positions the consensus keeps; see GENPLOTS.
      val(depth)

    output:
        // Named, not "*coveragePlot*": that glob also caught
        // coveragePlot_result.txt and published it as a plot.
        path("${meta.id}_coveragePlot.{html,png,svg}"), optional: true, emit: plots
        path("coveragePlot_result.txt"), optional: true, emit: result
    script:

    // bin/coverage_plot.py: a missing HTML plot fails the task, a missing PNG
    // or SVG is written to coveragePlot_result.txt, which GENPLOTS logs as a
    // warning. The script's docstring has the reasoning.
    """
    coverage_plot.py \\
        --bam ${bam} \\
        --sample ${meta.id} \\
        --threshold ${depth} \\
        --result-file coveragePlot_result.txt
    """
}

process snpPlot {
    tag "${meta.id}"
    publishDir { "${params.outDir}/${meta.id}_results/" }, mode: "copy"

    input:

      tuple val(meta), path("${meta.id}.depth${params.depth}.fa.algn")


    output:
        path("*snpPlot*")

    script:

    plot_name = "${meta.id}_snpPlot"


    """

    snipit ${meta.id}.depth${params.depth}.fa.algn -o ${plot_name} --solid-background
    snipit ${meta.id}.depth${params.depth}.fa.algn -o ${plot_name} -f svg --solid-background

    """
}
