version 1.0

workflow MityWorkflow {
  input {
    String sample_id
    File input_cram
    File input_crai                
    String reference = "hg38"
    Int cpu = 4                    
    Int memory_gb = 16             
    Int disk_gb = 50
  }

  call RunMity {
    input:
      sample_id = sample_id,
      input_cram = input_cram,
      input_crai = input_crai,
      reference = reference,
      cpu = cpu,
      memory_gb = memory_gb
  }

  output {
    File mity_vcf = RunMity.mity_vcf
    File mity_tbi = RunMity.mity_tbi
    File mity_report = RunMity.mity_report
  }
}

task RunMity {
  input {
    String sample_id
    File input_cram
    File input_crai
    String reference
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail

    echo "Running mity call for ~{sample_id}..."
    mkdir -p out

    mity call \
        --reference ~{reference} \
        --prefix ~{sample_id} \
        --output-dir out \
        --normalise \
        ~{input_cram} || { echo 'mity call failed'; exit 1; }

    echo "Running mity report for ~{sample_id}..."
    mity report \
        --prefix ~{sample_id} \
        --min_vaf 0.01 \
        --contig chrM \
        --output-dir out \
        out/~{sample_id}.normalise.vcf.gz || { echo 'mity report failed'; exit 1; }

    echo "mity completed for ~{sample_id}."
  >>>

  output {
    File mity_vcf = "out/~{sample_id}.normalise.vcf.gz"
    File mity_report = "out/~{sample_id}.mity.report.xlsx"  # adjust if report output differs
    File mity_tbi = "out/~{sample_id}.normalise.vcf.gz.tbi"
  }

  runtime {
    docker: "drmjc/mity"
    cpu: cpu
    memory: "~{memory_gb}G"
    disks: "local-disk ~{disk_gb} HDD"
  }
}