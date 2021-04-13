# Work before

I did some explorations on the BrainXcan visualization a while ago. 
The scripts are at `misc_data/prep_vis_data_*.R` and `misc_data/exploring_vis/`.

# Now

Now that we want to finalize the visualization pipeline as part of BrainXcan software.
Since some IDP meta data has been changed, I need to re-do the pipeline again using the
up-to-date meta files.

# What are they?

## T1 visualization

* Cortical (based on FAST (Harvard + Oxford))
* Subcortical total volumes (based on FIRST) 
* Subcortical gray matter (FAST ones): use the same as Subcortical FIRST 
* Cerebellum (based on the Diedrichsen cerebellar atlas)
* Brainstem and Total: **No support**

## dMRI visualization

* TBSS based
* ProbTrack based: **No annotation, so no support**

# Meta data for visualization

Save each visualization category separately.
There are three categories for T1 IDPs and one category for dMRI IDPs.

```
metadata = list(
    df_color = data.frame(
        color_code,
        IDP
    ),
    img = NIfTI-1 image with color_code
)
```

Also, save the background brain image.

```
metadata_bg = NIfTI-1 image with gray scale
```

Also, save the IDP sets with meta information on plotting.

```
meta_plot = list(
    name1 = category1, 
    name2 = category2, 
    ...
)
categoryi = list(
  IDP = list of IDPs,
  vis_rds = filename of metadata for visualization,
  slide_position = c(d1, d2, d3), # the slide position to show for the three dimensions
)
```

