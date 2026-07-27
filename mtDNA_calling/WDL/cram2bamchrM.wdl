version 1.0

workflow CramToChrMBam {
  input {
    String sample_id             
    File input_cram              
    File input_cram_index        
    File reference_fasta         
    File reference_fasta_index   
    Int n_cpu = 4
    Int memory_gb = 16           
    Int disk_gb = 25             
  }

  call ConvertCramToChrMBam {
    input:
      sample_id = sample_id,
      input_cram = input_cram,
      input_cram_index = input_cram_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      n_cpu = n_cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File sorted_bam = ConvertCramToChrMBam.sorted_bam
    File sorted_bam_index = ConvertCramToChrMBam.sorted_bam_index
  }
}

task ConvertCramToChrMBam {
  input {
    String sample_id
    File input_cram
    File input_cram_index        
    File reference_fasta
    File reference_fasta_index
    Int n_cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail

    
    if [ ! -f "~{input_cram_index}" ]; then
        echo "ERROR: CRAI index not found!" 1>&2
        exit 1
    fi

    export TMPDIR=$(pwd)/tmp
    mkdir -p $TMPDIR

    samtools view -T ~{reference_fasta} -b ~{input_cram} chrM > ~{sample_id}.chrM.bam

    samtools sort -@ ~{n_cpu} -T $TMPDIR/~{sample_id}.tmp -o ~{sample_id}.chrM.sorted.bam ~{sample_id}.chrM.bam

    samtools index ~{sample_id}.chrM.sorted.bam

    rm ~{sample_id}.chrM.bam
    rm -rf $TMPDIR
  >>>

  output {
    File sorted_bam = "~{sample_id}.chrM.sorted.bam"
    File sorted_bam_index = "~{sample_id}.chrM.sorted.bam.bai"
  }

  runtime {
    cpu: n_cpu
    memory: "~{memory_gb}G"
    disks: "local-disk ~{disk_gb} HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
  }
}