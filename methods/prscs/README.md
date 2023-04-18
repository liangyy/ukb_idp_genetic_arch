Here I am implementing the PRS-cs method on the brain IDP traits

- Starting with the individual level data split it into training and testing
- Run GWAS using the training dataset to obtain the summary stats
- Run the PRS-sc method with the sum stats to obtain the weights
- Using the weights predict the IDP traits
- Evaluate the performance of the model
- Compare the PRS-cs model performance with Ridge and Lasso

```{r}
# one to one comparison
ridge_pred %>% left_join(t1_corr %>% bind_rows(dmri_corr), by = "IDP") %>% 
  mutate(Spearman = sign(Spearman) *(Spearman^2)) %>% 
  mutate(Spearman_prscs = sign(Spearman_prscs) *(Spearman_prscs^2)) %>% 
  ggplot(., aes(x = Spearman_prscs, y = Spearman)) + 
    geom_point() + 
    geom_abline(intercept = 0,slope = 1, color="red") +
    geom_hline(yintercept = 0, color="red", linetype = "dotted") +
    geom_vline(xintercept = 0, color="red", linetype = "dotted") +
    #facet_grid(rows = vars(model_name), cols = vars(method)) +
    facet_wrap(model_name~., ncol = 2, scales = "fixed") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1,
          strip.text = element_text(size = 18, face = "bold")) +
    labs(title = "Comparing BrainXcan e-net, ridge vs PRScs models performance",
         y=expression(italic(R^2)),x=expression(paste("PRS-cs ", R^2)))
```
