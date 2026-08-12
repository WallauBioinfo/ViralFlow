process compileOutputs{
  publishDir "${workflow.outputDir}/COMPILED_OUTPUT/", mode: "copy"
  label "singlethread"
  
  input:
    val(go)
    val(virus_tag)

  output:
    path("*")
  script:
    """
    compileOutput.py -dD ${workflow.outputDir} \
                            -oD ./ \
                            --depth ${params.depth} \
                            -virus_tag ${virus_tag}
    """
}
