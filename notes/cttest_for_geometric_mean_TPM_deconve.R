annotation_table <-  read.csv(get_task_outputs("wf-2c26fb32e8",0)["tcga_annotation_data.csv"])

manipulated_deconve_data <- coef(readRDS(get_workflow_outputs('wf-9d21bd3f11')))
head(manipulated_deconve_data[1:5,1:5])

# compare 2 groups
unique(annotation_table$group)
group_a_name <- "bladder cancer_epithelial"
group_b_name <- "bladder cancer_inflamed"
group_a_data <- manipulated_deconve_data[,colnames(manipulated_deconve_data) %in% annotation_table$sample_id[annotation_table$group == group_a_name]]
group_b_data <- manipulated_deconve_data[,colnames(manipulated_deconve_data) %in% annotation_table$sample_id[annotation_table$group == group_b_name]]

# run CT test
new_dataset <- cbind(group_a_data, group_b_data)
group = c(rep(group_a_name, dim(group_a_data)[2]), rep(group_b_name, dim(group_b_data)[2]))
res_manipulated <- cytoreason.deconvolution::service_ct_test(new_dataset, group = group)
plot_df <- res_manipulated$pvalues %>%
  dplyr::mutate(feature_id = reorder(factor(rownames(res_manipulated$pvalues)), estimate),
                Direction = ifelse(estimate > 0, "Up", "Down"))

print(ggplot(plot_df, aes(x = feature_id, y = estimate, fill = Direction)) +
        geom_bar(stat = "identity") +
        xlab("Cell type") + ggtitle(paste0("CT-test for manipulated data ", group_a_name," vs ",group_b_name))+
        scale_fill_brewer(palette = "Set1", direction = -1, aes(colour = Direction)) +
        coord_flip())

###############################

raw_deconve_data <- coef(readRDS(get_workflow_outputs('wf-b0067e79f7')))
head(raw_deconve_data[1:5,1:5])

# compare 2 groups
group_a_data <- raw_deconve_data[,colnames(raw_deconve_data) %in% annotation_table$sample_id[annotation_table$group == group_a_name]]
group_b_data <- raw_deconve_data[,colnames(raw_deconve_data) %in% annotation_table$sample_id[annotation_table$group == group_b_name]]

# run CT test
new_dataset <- cbind(group_a_data, group_b_data)
group = c(rep(group_a_name, dim(group_a_data)[2]), rep(group_b_name, dim(group_b_data)[2]))
res_raw <- cytoreason.deconvolution::service_ct_test(new_dataset, group = group)
plot_df <- res_raw$pvalues %>%
  dplyr::mutate(feature_id = reorder(factor(rownames(res_raw$pvalues)), estimate),
                Direction = ifelse(estimate > 0, "Up", "Down"))
print(ggplot(plot_df, aes(x = feature_id, y = estimate, fill = Direction)) +
        geom_bar(stat = "identity") +
        xlab("Cell type") + ggtitle(paste0("CT-test for raw data ", group_a_name," vs ",group_b_name))+
        scale_fill_brewer(palette = "Set1", direction = -1, aes(colour = Direction)) +
        coord_flip())



