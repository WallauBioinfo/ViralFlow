//enable dsl 2
nextflow.enable.dsl = 2
include {run_porechop} from '../modules/runPorechop.nf'
include {run_minimap2} from '../modules/runMinimap2.nf'
include {run_amplicon_clip} from '../modules/runAmpliconClip.nf'
include {run_bam_utils} from '../modules/runBamUtils.nf'
include {run_clair3} from '../modules/runClair3.nf'
include {run_bcftools; run_bcftools_consensus} from '../modules/runBcftools.nf'

workflow NANOPORE {
    take:
        reads_ch // tuple (meta, fastq)
        ref // path to reference genome

    main:
    
    // remove adapters (porechop)
    run_porechop(reads_ch)
    reads_ch = run_porechop.out

    // do alignment (minimap2) 
    run_minimap2(reads_ch, ref)
    bams_ch = run_minimap2.out // tuple (meta, sorted_bam, bai)
    
    // do variant calling (clair3)
    run_clair3(bams_ch, ref, params.clair3_chunk_size, params.clair3_qual, params.mapping_quality, params.clair3_model)

    // normlalize indes and filter variants (bcftools)
    run_bcftools(run_clair3.out, ref, params.af_threshold)
    
    run_bcftools.out.map{meta, vcf -> tuple(meta.id, meta, vcf)}
        .set { vcf_ch } // tuple (id, meta, vcf)

    bams_ch.map { meta, sorted_bam -> tuple(meta.id, sorted_bam)} // tuple (id, sorted_bam)
    .join(vcf_ch) // tuple (id, sorted_bam, meta, vcf)
    .map { _id, sorted_bam, meta, vcf ->  tuple(meta, vcf, sorted_bam[0])}
    .set { vcf_bam_ch } // tuple (meta, vcf, sorted_bam)

    // call consensus sequence (bcftools consensus)
    run_bcftools_consensus(vcf_bam_ch, ref, params.np_min_depth)

    emit:
        bams_ch = run_minimap2.out //meta, sorted_bam
}

workflow {
    // Define the input files
    reads_ch = parse_mnf(params.mnf) // tuple (meta, fastq)

    log.info("${params.base_container} ${params.mnf} ${params.outDir}")
    // run workflow
    NANOPORE(reads_ch, params.ref)
}


def parse_mnf(mnf) {
    def mnf_rows = Channel.fromPath(mnf)
            .splitCsv(header: true, sep: ',')
            .map { row -> 
                    // set meta
                    def meta = [id: row.sample_id]

                    // declare channel shape
                    tuple(meta, row.fastq)
                 }

    return mnf_rows
}