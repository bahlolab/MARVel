version 1.0

workflow MitoHPC_Workflow {
  input {
    File cram_input
    File crai_input
    String sample_name
    String mitohpc_docker = "dpuiu1/mitohpc"
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb = 50
    Int? preemptible_tries_hpc
	
  }

  call RunMitoHPC {
    input:
      cram_input = cram_input,
      crai_input = crai_input,
      sample_name = sample_name,
      mitohpc_docker = mitohpc_docker,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      preemptible_tries_hpc = preemptible_tries_hpc
  }

  output {
    Array[File] all_outputs = RunMitoHPC.all_outputs
  }
}

task RunMitoHPC {
  input {
    File cram_input
    File crai_input
    String sample_name
    String mitohpc_docker
    Int cpu
    Int memory_gb
    Int disk_gb
    Int? preemptible_tries_hpc
  }

  command <<<
    set -euo pipefail

    mkdir -p bams out
    mv ~{cram_input} bams/
    mv ~{crai_input} bams/

    . /MitoHPC/scripts/init.sh

    # Export variables 
    export HP_ADIR="bams"
    export HP_ODIR="out"
    export HP_L=

    find "$HP_ADIR/" \( -name '*.bam' -o -name '*.cram' \) -readable | ls2in.pl -out \$HP_ODIR | sort -V > \$HP_IN
    /MitoHPC/scripts/run.sh | bash
  >>>

  runtime {
    docker: mitohpc_docker
    cpu: cpu
    memory: "~{memory_gb}G"
    disks: "local-disk ~{disk_gb} HDD",
    preemptible: select_first([preemptible_tries_hpc, 5])
  }

  output {
    Array[File] all_outputs = glob("out{,/**}/*")
  }
}
