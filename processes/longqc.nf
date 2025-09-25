/*process LONGQC {
   //conda "./environment.yaml"
   //container "docker.io/cymbopogon/longqc"
   publishDir "$params.publish_dir/longqc", mode: 'copy'

   input:
   tuple val(sample), path(reads)

   output:
   tuple val(sample), path("longqc/*"), emit: longqc

   script:
   """
   longqc sampleqc -x $params.seq_type -o ${sample} -p $params.cpus ${reads}
   """
}
*/
