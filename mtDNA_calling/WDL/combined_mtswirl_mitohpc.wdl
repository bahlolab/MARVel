version 1.0

import "https://raw.githubusercontent.com/WEHI-ResearchComputing/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/AlignAndCallR1_v2_5_Single.wdl" as AlignAndCallR1_Single
import "https://raw.githubusercontent.com/WEHI-ResearchComputing/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/AlignAndCallR2_v2_5_Single.wdl" as AlignAndCallR2_Single
import "https://raw.githubusercontent.com/WEHI-ResearchComputing/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/LiftoverTools_v2_5_Single.wdl" as LiftoverTools_Single
import "https://raw.githubusercontent.com/WEHI-ResearchComputing/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/ProduceSelfReferenceFiles_v2_5_Single.wdl" as ProduceSelfReferenceFiles_Single
import "https://raw.githubusercontent.com/WEHI-ResearchComputing/mtSwirl/master/WDL/v2.5_MongoSwirl_Single/MongoTasks_v2_5_Single.wdl" as MongoTasks_Single

workflow MtToolsAndConversion {

  meta {
    description: "Takes in an hg38 bam or cram and outputs VCF of SNP/Indel calls on the mitochondria."
    allowNestedInputs: true
  }

  input {
    ## common input for both tools ##
    File wgs_aligned_input_bam_or_cram
    File wgs_aligned_input_bam_or_cram_index
    String sample_name

    ## mtSwirl input ##
    File mt_interval_list
    File nuc_interval_list

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File blacklisted_sites
    File blacklisted_sites_index

    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list

    File HailLiftover
    File FaRenamingScript
    File CheckVariantBoundsScript
    File CheckHomOverlapScript
    File JsonTools

    Boolean force_manual_download
    String? requester_pays_project
    String? m2_extra_args
    String? m2_filter_extra_args
    String? printreads_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float? verifyBamID
    Boolean compress_output_vcf = false
    Boolean compute_numt_coverage = false
    Boolean use_haplotype_caller_nucdna = true
    Int haplotype_caller_nucdna_dp_lower_bound = 10

    # Some CRAMs (e.g., AoU) contain an XQ tag that doesn't play well with gatk RevertSam.
    # Enable this flag to avoid using this tag by skipping the restore hardclips step.
    Boolean skip_restore_hardclips = false

    #Docker and version arguments
    String gatk_version = "4.2.6.0"
    File? gatk_override
    String? gatk_docker_override
    String ucsc_docker
    String genomes_cloud_docker
    String haplochecker_docker
    String gatk_samtools_docker

    #Optional runtime arguments
    Int? printreads_mem
    Int? lift_coverage_mem
    Int? n_cpu_subsetbam
    Int? n_cpu_m2_hc_lift
    Int? n_cpu_bwa
    Int? preemptible_tries
    Int? preemptible_tries_hpc 

    ## mitoHPC inputs ##
    String mitohpc_docker = "dpuiu1/mitohpc"
    Int hpc_cpu = 4
    Int hpc_memory_gb = 16
    Int hpc_disk_gb = 50

    ## cram2chrM conversion inputs ##  
    Int cram2bam_cpu = 1
    Int cram2bam_mem = 8           
    Int cram2bam_disk = 50
  }

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS hg38 bam or cram"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    mt_interval_list: "Picard style interval list file, with header and single interval representing chrM, eg chrM 1 16569 + ., and putative NUMT intervals."
  }

  String self_ref_suffix = ".self.ref"
  
  ## run cram chrM extraction ##
  call ConvertCramToChrMCram {
    input:
      sample_id = sample_name,
      input_cram = wgs_aligned_input_bam_or_cram,
      input_cram_index = wgs_aligned_input_bam_or_cram_index,
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      cram2bam_cpu = cram2bam_cpu,
      cram2bam_mem = cram2bam_mem,
      cram2bam_disk = cram2bam_disk,
      preemptible_tries_hpc = preemptible_tries_hpc
  }
  ## run mitoHPC ##
  call RunMitoHPC {
    input:
      cram_input = wgs_aligned_input_bam_or_cram,
      crai_input = wgs_aligned_input_bam_or_cram_index,
      sample_name = sample_name,
      mitohpc_docker = mitohpc_docker,
      cpu = hpc_cpu,
      memory_gb = hpc_memory_gb,
      disk_gb = hpc_disk_gb,
      preemptible_tries_hpc = preemptible_tries_hpc
  }
## run mtSwirl ##
  call MongoTasks_Single.MongoSubsetBamToChrMAndRevert as SubsetBamToChrMAndRevert { # modify + include the chrM + NUMT
    input:
      input_bam = wgs_aligned_input_bam_or_cram,
      input_bai = wgs_aligned_input_bam_or_cram_index,
      sample_name = sample_name,

      mt_interval_list = mt_interval_list,
      nuc_interval_list = nuc_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      requester_pays_project = requester_pays_project,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_samtools_docker,
      gatk_version = gatk_version,
      printreads_extra_args = printreads_extra_args,
      force_manual_download = force_manual_download,
      read_length = max_read_length,
      skip_restore_hardclips = skip_restore_hardclips,
      coverage_cap = 100000,
      mem = printreads_mem,
      n_cpu = n_cpu_subsetbam,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCallR1_Single.AlignAndCallR1 as AlignAndCallR1 {
    input:
      input_bam = SubsetBamToChrMAndRevert.output_bam,
      input_bai = SubsetBamToChrMAndRevert.output_bai,
      sample_name = sample_name,

      mt_interval_list = mt_interval_list,
      nuc_interval_list = nuc_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_mean_coverage = SubsetBamToChrMAndRevert.mean_coverage,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      m2_extra_args = m2_extra_args,
      m2_filter_extra_args = m2_filter_extra_args,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      verifyBamID = verifyBamID,
      compress_output_vcf = compress_output_vcf,
      max_read_length = max_read_length,
      use_haplotype_caller_nucdna = use_haplotype_caller_nucdna,
      hc_dp_lower_bound = haplotype_caller_nucdna_dp_lower_bound,
      preemptible_tries = preemptible_tries,
      haplochecker_docker = haplochecker_docker,
      n_cpu = n_cpu_m2_hc_lift
  }

  call ProduceSelfReferenceFiles_Single.ProduceSelfReferenceFiles as ProduceSelfRefFiles {
    input:
      sample_name = sample_name,
      suffix = self_ref_suffix,
      mtdna_variants = AlignAndCallR1.split_vcf,
      nuc_variants = AlignAndCallR1.split_nuc_vcf,
      
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_dict = mt_dict,
      mt_interval_list = mt_interval_list,
      non_control_region_interval_list = non_control_region_interval_list,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      nuc_interval_list = nuc_interval_list,
      reference_name = "reference",
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      n_shift = 8000,
      compute_numt_coverage = compute_numt_coverage,
      FaRenamingScript = FaRenamingScript,
      CheckVariantBoundsScript = CheckVariantBoundsScript,
      CheckHomOverlapScript = CheckHomOverlapScript,
      genomes_cloud_docker = genomes_cloud_docker,
      ucsc_docker = ucsc_docker,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCallR2_Single.AlignAndCallR2 as AlignAndCallR2 {
    input:
      unmapped_bam = SubsetBamToChrMAndRevert.unmapped_bam,
      sample_name = sample_name,
      suffix = self_ref_suffix,

      mt_interval_list = ProduceSelfRefFiles.mt_interval_list_self,

      mt_self = ProduceSelfRefFiles.mt_self,
      mt_self_index = ProduceSelfRefFiles.mt_self_index,
      mt_self_dict = ProduceSelfRefFiles.mt_self_dict,
      self_cat = ProduceSelfRefFiles.mt_andNuc_self,
      self_cat_index = ProduceSelfRefFiles.mt_andNuc_self_index,
      self_cat_dict = ProduceSelfRefFiles.mt_andNuc_self_dict,
      mt_self_shifted = ProduceSelfRefFiles.mt_shifted_self,
      mt_self_shifted_index = ProduceSelfRefFiles.mt_shifted_self_index,
      mt_self_shifted_dict = ProduceSelfRefFiles.mt_shifted_self_dict,
      self_shifted_cat = ProduceSelfRefFiles.mt_andNuc_shifted_self,
      self_shifted_cat_index = ProduceSelfRefFiles.mt_andNuc_shifted_self_index,
      self_shifted_cat_dict = ProduceSelfRefFiles.mt_andNuc_shifted_self_dict,
      shift_back_chain = ProduceSelfRefFiles.self_shift_back_chain,

      force_call_vcf = ProduceSelfRefFiles.force_call_vcf,
      force_call_vcf_idx = ProduceSelfRefFiles.force_call_vcf_idx,
      force_call_vcf_shifted = ProduceSelfRefFiles.force_call_vcf_shifted,
      force_call_vcf_shifted_idx = ProduceSelfRefFiles.force_call_vcf_shifted_idx,

      non_control_interval = ProduceSelfRefFiles.non_control_interval_self,
      control_shifted = ProduceSelfRefFiles.control_shifted_self,
      blacklisted_sites = ProduceSelfRefFiles.blacklisted_sites_self,
      blacklisted_sites_index = ProduceSelfRefFiles.blacklisted_sites_index_self,
      
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      gatk_version = gatk_version,
      m2_extra_args = select_first([m2_extra_args," "]),
      m2_filter_extra_args = m2_filter_extra_args,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      verifyBamID = verifyBamID,
      compress_output_vcf = compress_output_vcf,
      max_read_length = max_read_length,
      preemptible_tries = preemptible_tries,
      hasContamination = AlignAndCallR1.hasContamination,
      contamination_major = AlignAndCallR1.contamination_major,
      contamination_minor = AlignAndCallR1.contamination_minor,
      n_cpu_bwa = n_cpu_bwa,
      n_cpu = n_cpu_m2_hc_lift
  }

  call MongoTasks_Single.MongoLiftoverVCFAndGetCoverage as LiftOverAfterSelf {
    input:
      sample_name = sample_name,
      original_filtered_vcf = ProduceSelfRefFiles.ref_homoplasmies_vcf,
      new_self_ref_vcf = AlignAndCallR2.split_vcf,
      reversed_hom_ref_vcf = ProduceSelfRefFiles.force_call_vcf_filters,

      mt_self = ProduceSelfRefFiles.mt_self,
      mt_self_index = ProduceSelfRefFiles.mt_self_index,
      mt_self_dict = ProduceSelfRefFiles.mt_self_dict,
      mt_self_shifted = ProduceSelfRefFiles.mt_shifted_self,
      mt_self_shifted_index = ProduceSelfRefFiles.mt_shifted_self_index,
      mt_self_shifted_dict = ProduceSelfRefFiles.mt_shifted_self_dict,
      chain_self_to_ref = ProduceSelfRefFiles.self_to_ref_chain,
      chain_ref_to_self = ProduceSelfRefFiles.ref_to_self_chain,

      input_bam_regular_ref = AlignAndCallR2.mt_aligned_bam,
      input_bam_regular_ref_index = AlignAndCallR2.mt_aligned_bai,
      input_bam_shifted_ref = AlignAndCallR2.mt_aligned_shifted_bam,
      input_bam_shifted_ref_index = AlignAndCallR2.mt_aligned_shifted_bai,
      self_control_region_shifted_reference_interval_list = ProduceSelfRefFiles.control_shifted_self,
      self_non_control_region_interval_list = ProduceSelfRefFiles.non_control_interval_self,
      
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      HailLiftover = HailLiftover,
      self_suffix = self_ref_suffix,

      n_cpu = n_cpu_m2_hc_lift,
      genomes_cloud_docker = genomes_cloud_docker,
      preemptible_tries = preemptible_tries
  }

  call MongoTasks_Single.MongoLiftoverSelfAndCollectOutputs as LiftoverSelfCoverage {
    input:
      sample_name = sample_name,
      self_ref_table = LiftOverAfterSelf.self_coverage_table,
      chain = ProduceSelfRefFiles.self_to_ref_chain,
      homoplasmic_deletions_coverage = LiftOverAfterSelf.gap_coverage,
      
      liftover_table = LiftOverAfterSelf.liftoverStats,
      mean_coverage = AlignAndCallR2.mean_coverage,
      median_coverage = AlignAndCallR2.median_coverage,
      major_haplogroup = AlignAndCallR1.major_haplogroup,
      contamination = AlignAndCallR1.contamination,
      nuc_variants_pass = AlignAndCallR1.nuc_variants_pass,
      n_reads_unpaired_dropped = SubsetBamToChrMAndRevert.reads_dropped,
      nuc_variants_dropped = ProduceSelfRefFiles.nuc_variants_dropped,
      mtdna_consensus_overlaps = ProduceSelfRefFiles.mtdna_consensus_overlaps,
      nuc_consensus_overlaps = ProduceSelfRefFiles.nuc_consensus_overlaps,
      ucsc_docker = ucsc_docker,
      preemptible_tries = preemptible_tries
  }

  if (compute_numt_coverage) {
    call NucCoverageAtEveryBase {
      input:
        input_bam_regular_ref = SubsetBamToChrMAndRevert.output_bam,
        input_bam_regular_ref_index = SubsetBamToChrMAndRevert.output_bai,
        input_bam_self_ref = AlignAndCallR2.nuc_mt_aligned_bam,
        input_bam_self_ref_index = AlignAndCallR2.nuc_mt_aligned_bai,
        input_bam_self_ref_shifted = AlignAndCallR2.nuc_mt_shifted_aligned_bam,
        input_bam_self_ref_shifted_index = AlignAndCallR2.nuc_mt_shifted_aligned_bai,
        nuc_interval_list = nuc_interval_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        self_nuc_interval_list = ProduceSelfRefFiles.nuc_interval_list_self,
        self_fasta = ProduceSelfRefFiles.mt_andNuc_self,
        self_fasta_index = ProduceSelfRefFiles.mt_andNuc_self_index,
        self_dict = ProduceSelfRefFiles.mt_andNuc_self_dict,
        self_nuc_interval_list_shifted = ProduceSelfRefFiles.nuc_interval_list_shifted_self,
        self_shifted_fasta = ProduceSelfRefFiles.mt_andNuc_shifted_self,
        self_shifted_fasta_index = ProduceSelfRefFiles.mt_andNuc_shifted_self_index,
        self_shifted_dict = ProduceSelfRefFiles.mt_andNuc_shifted_self_dict,
        preemptible_tries = preemptible_tries
    }

    call LiftoverTools_Single.LiftOverAndJoinCoverage as LiftOverAndJoinCoverage {
      input:
        ref_table = NucCoverageAtEveryBase.table_old_ref,
        self_table = NucCoverageAtEveryBase.table_new_self,
        self_table_shifted = NucCoverageAtEveryBase.table_new_self_shifted,
        chain = ProduceSelfRefFiles.nuc_self_to_ref_chain,
        mem = lift_coverage_mem,
        ucsc_docker = ucsc_docker,
        preemptible_tries = preemptible_tries
    }
  }
call TarMitoHPCOutput {
  input:
    all_outputs = RunMitoHPC.all_outputs,
    sample_id = sample_name
}

call TarMtswirlOutputs {
  input:
    final_vcf = LiftOverAfterSelf.liftover_r2_final_vcf,
    final_base_level_coverage_metrics = LiftoverSelfCoverage.reference_coverage,
    sample_id = sample_name
}
  output {
    File subset_bam = SubsetBamToChrMAndRevert.output_bam
    File subset_bai = SubsetBamToChrMAndRevert.output_bai
    File idxstats_metrics = SubsetBamToChrMAndRevert.idxstats_metrics
    File yield_metrics = SubsetBamToChrMAndRevert.yield_metrics
    File flagstat = SubsetBamToChrMAndRevert.flagstat

    File r1_vcf = AlignAndCallR1.out_vcf
    File r1_vcf_index = AlignAndCallR1.out_vcf_index
    File r1_nuc_vcf = AlignAndCallR1.nuc_vcf
    File r1_nuc_vcf_index = AlignAndCallR1.nuc_vcf_index
    File r1_nuc_vcf_unfiltered = AlignAndCallR1.nuc_vcf_unfiltered
    File r1_split_vcf = AlignAndCallR1.split_vcf
    File r1_split_vcf_index = AlignAndCallR1.split_vcf_index

    File self_mt_aligned_bam = AlignAndCallR2.mt_aligned_bam
    File self_mt_aligned_bai = AlignAndCallR2.mt_aligned_bai
    File self_ref_vcf = AlignAndCallR2.out_vcf
    File self_ref_vcf_index = AlignAndCallR2.out_vcf_idx
    File self_ref_split_vcf = AlignAndCallR2.split_vcf
    File self_ref_split_vcf_index = AlignAndCallR2.split_vcf_idx
    File self_base_level_coverage_metrics = LiftOverAfterSelf.self_coverage_table
    File self_reference_fasta = ProduceSelfRefFiles.mt_self
    File reference_to_self_ref_chain = ProduceSelfRefFiles.ref_to_self_chain
    File self_control_region_shifted = ProduceSelfRefFiles.control_shifted_self
    File self_non_control_region = ProduceSelfRefFiles.non_control_interval_self

    File liftover_fix_pipeline_log = LiftOverAfterSelf.liftover_r2_log
    File stats_outputs = LiftoverSelfCoverage.table

    File final_vcf = LiftOverAfterSelf.liftover_r2_final_vcf
    File final_rejected_vcf = LiftOverAfterSelf.liftover_r2_rejected_vcf
    File final_base_level_coverage_metrics = LiftoverSelfCoverage.reference_coverage
    File? numt_base_level_coverage = LiftOverAndJoinCoverage.reference_coverage
    
    File input_vcf_for_haplochecker = AlignAndCallR1.input_vcf_for_haplochecker
    File duplicate_metrics = AlignAndCallR2.duplicate_metrics
    File coverage_metrics = AlignAndCallR2.coverage_metrics
    File theoretical_sensitivity_metrics = AlignAndCallR2.theoretical_sensitivity_metrics
    File contamination_metrics = AlignAndCallR1.contamination_metrics

    # liftover stats – the full set are outputted in the stats outputs file above
    Int success_liftover_variants = LiftOverAfterSelf.n_liftover_r2_pass
    Int failed_liftover_variants = LiftOverAfterSelf.n_liftover_r2_failed
    Int fixed_liftover_variants = LiftOverAfterSelf.n_liftover_r2_fixed
    Int n_liftover_r2_left_shift = LiftOverAfterSelf.n_liftover_r2_left_shift
    Int n_liftover_r2_injected_from_success = LiftOverAfterSelf.n_liftover_r2_injected_from_success
    Int n_liftover_r2_ref_insertion_new_haplo = LiftOverAfterSelf.n_liftover_r2_ref_insertion_new_haplo
    Int n_liftover_r2_failed_het_dele_span_insertion_boundary = LiftOverAfterSelf.n_liftover_r2_failed_het_dele_span_insertion_boundary
    Int n_liftover_r2_failed_new_dupes_leftshift = LiftOverAfterSelf.n_liftover_r2_failed_new_dupes_leftshift
    Int n_liftover_r2_het_ins_sharing_lhs_hom_dele = LiftOverAfterSelf.n_liftover_r2_het_ins_sharing_lhs_hom_dele
    Int n_liftover_r2_spanning_complex = LiftOverAfterSelf.n_liftover_r2_spanning_complex
    Int n_liftover_r2_spanningfixrhs_sharedlhs = LiftOverAfterSelf.n_liftover_r2_spanningfixrhs_sharedlhs
    Int n_liftover_r2_spanningfixlhs_upstream = LiftOverAfterSelf.n_liftover_r2_spanningfixlhs_upstream
    Int n_liftover_r2_repaired_success = LiftOverAfterSelf.n_liftover_r2_repaired_success

    # other stats and constants; also in stats_outputs
    Int mean_coverage = AlignAndCallR2.mean_coverage
    Float median_coverage = AlignAndCallR2.median_coverage
    String major_haplogroup = AlignAndCallR1.major_haplogroup
    Float contamination = AlignAndCallR1.contamination
    Int nuc_variants_pass = AlignAndCallR1.nuc_variants_pass
    Int n_reads_unpaired_dropped = SubsetBamToChrMAndRevert.reads_dropped
    Int nuc_variants_dropped = ProduceSelfRefFiles.nuc_variants_dropped
    Int mtdna_consensus_overlaps = ProduceSelfRefFiles.mtdna_consensus_overlaps
    Int nuc_consensus_overlaps = ProduceSelfRefFiles.nuc_consensus_overlaps

    # mitoHPC output
    Array[File] all_outputs = RunMitoHPC.all_outputs

    ## cram extract chrM output ##
    File chrM_cram = ConvertCramToChrMCram.chrM_cram
    File chrM_crai = ConvertCramToChrMCram.chrM_crai

    ## tarball output (tar.gz) ##
    File mitohpc_tar = TarMitoHPCOutput.mitohpc_tar
    File mtswirl_tar = TarMtswirlOutputs.mtswirl_tar

  }
}

task NucCoverageAtEveryBase {
  input {
    File input_bam_regular_ref
    File input_bam_regular_ref_index
    File input_bam_self_ref
    File input_bam_self_ref_index
    File input_bam_self_ref_shifted
    File input_bam_self_ref_shifted_index
    File nuc_interval_list
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? self_nuc_interval_list
    File self_fasta
    File self_fasta_index
    File self_dict
    File? self_nuc_interval_list_shifted
    File self_shifted_fasta
    File self_shifted_fasta_index
    File self_shifted_dict

    Int? preemptible_tries
  }
  Int disk_size = ceil(size(input_bam_regular_ref, "GB") + size(input_bam_self_ref, "GB") + size(ref_fasta, "GB") * 2) + 20

  meta {
    description: "Mainly for QC to understand how NUMT coverage changes with remapping."
  }
  command <<<
    set -e

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=~{input_bam_regular_ref} \
      R=~{ref_fasta} \
      PER_BASE_COVERAGE=reference_pre_realignment.tsv \
      O=reference_pre_realignment.metrics \
      TI=~{nuc_interval_list} \
      BI=~{nuc_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=~{input_bam_self_ref} \
      R=~{self_fasta} \
      PER_BASE_COVERAGE=self_post_realignment.tsv \
      O=self_post_realignment.metrics \
      TI=~{self_nuc_interval_list} \
      BI=~{self_nuc_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=~{input_bam_self_ref_shifted} \
      R=~{self_shifted_fasta} \
      PER_BASE_COVERAGE=self_post_realignment_shifted.tsv \
      O=self_post_realignment_shifted.metrics \
      TI=~{self_nuc_interval_list_shifted} \
      BI=~{self_nuc_interval_list_shifted} \
      COVMAX=20000 \
      SAMPLE_SIZE=1
  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "1200 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File table_old_ref = "reference_pre_realignment.tsv"
    File table_new_self = "self_post_realignment.tsv"
    File table_new_self_shifted = "self_post_realignment_shifted.tsv"
  }
}

## Run MitoHPC ##
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
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: select_first([preemptible_tries_hpc, 5])
  }

  output {
    Array[File] all_outputs = glob("out{,/**}/*") 
  }
}

task ConvertCramToChrMCram {
  input {
    String sample_id
    File input_cram
    File input_cram_index        
    File reference_fasta
    File reference_fasta_index
    Int cram2bam_cpu
    Int cram2bam_mem
    Int cram2bam_disk
    Int? preemptible_tries_hpc
  }

  command <<<
    set -euo pipefail

    if [ ! -f "~{input_cram_index}" ]; then
        echo "ERROR: CRAI index not found!" 1>&2
        exit 1
    fi

   samtools view -@ ~{cram2bam_cpu} -T ~{reference_fasta} -C -o ~{sample_id}.chrM.cram ~{input_cram} chrM
   samtools index ~{sample_id}.chrM.cram
  >>>

  output {
    File chrM_cram = "~{sample_id}.chrM.cram"
    File chrM_crai = "~{sample_id}.chrM.cram.crai"
  }

  runtime {
    cpu: cram2bam_cpu
    memory: "~{cram2bam_mem}G"
    disks: "local-disk ~{cram2bam_disk} HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    preemptible: select_first([preemptible_tries_hpc, 5])
  }
}

task TarMitoHPCOutput {
  input {
    Array[File] all_outputs
    String sample_id
  }

  command <<<
    set -euo pipefail
    tar -czvf ~{sample_id}.mitohpc_outputs.tar.gz ~{sep=" " all_outputs}
  >>>

  output {
    File mitohpc_tar = "~{sample_id}.mitohpc_outputs.tar.gz"
  }

  runtime {
    docker: "ubuntu:22.04"
    memory: "4G"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}

task TarMtswirlOutputs {
  input {
    File final_vcf
    File final_base_level_coverage_metrics
    String sample_id
  }

  command <<<
    set -euo pipefail
    tar -czvf ~{sample_id}.mtswirl_outputs.tar.gz ~{final_vcf} ~{final_base_level_coverage_metrics}
  >>>

  output {
    File mtswirl_tar = "~{sample_id}.mtswirl_outputs.tar.gz"
  }

  runtime {
    docker: "ubuntu:22.04"
    memory: "4G"
    cpu: 1
    disks: "local-disk 10 HDD"
  }
}
