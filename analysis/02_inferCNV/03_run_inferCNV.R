library(infercnv)
infercnv_obj = CreateInfercnvObject(
                    raw_counts_matrix="/raid1/heyao/project/scHCC-tumor/analyses/02_inferCNV/data_for_run/counts_stroma_20210421.tsv",
                    annotations_file="/raid1/heyao/project/scHCC-tumor/analyses/02_inferCNV/data_for_run/annotation_stroma_20210421.tsv",
                    ref_group_names = NULL,
                    delim="\t",
                    gene_order_file="/raid1/heyao/project/scHCC-tumor/analyses/02_inferCNV/data_for_run/gene_info_stroma.tsv",
                    )




infercnv_obj = infercnv::run(infercnv_obj,
                            #analysis_mode = 'subclusters',
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                    
                             out_dir="/raid1/heyao/project/scHCC-tumor/analyses/02_inferCNV/out_inferCNV/hepatoctyes_sd1.75_stroma",  # dir is auto-created for storing outputs
                             scale_dat = FALSE,
                             cluster_by_groups=FALSE,   # cluster
                             denoise=TRUE,
                             HMM=FALSE,
                             HMM_type = 'i6', 
                             output_format='pdf',
                             tumor_subcluster_pval = 0.05,
                             sd_amplifier=1.75,  # sets midpoint for logistic
                             num_threads = 32,
                             BayesMaxPNormal = 0.8, 
                             noise_logistic=TRUE, # turns gradient filtering on
                             tumor_subcluster_partition_method='qnorm'
                             )