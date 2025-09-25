process FASTQC {
   //conda "./environment.yaml"
   //container "docker.io/staphb/fastqc"
   publishDir "$params.publish_dir/fastqc", mode: 'copy'

   input:
   tuple val(sample), path(reads)

   output:
   path '*_fastqc.{zip,html}', emit: report

   script:
   """
   fastqc --quiet --threads $params.cpus ${reads}
   """
}
