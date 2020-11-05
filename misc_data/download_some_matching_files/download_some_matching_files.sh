wget https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/-/raw/master/bb_IDP/list.txt

tmp='FS SWI T1 T2_FLAIR dMRI rfMRI_AMP rfMRI_CONN tfMRI'
for FF in $tmp
do
  wget https://git.fmrib.ox.ac.uk/falmagro/ukb_unconfound_v2/-/raw/master/data/GROUPS_IDPs_7_groups/$FF.txt
done

cat list.txt | cut -d\" -f2 | sed 's#^ ##g' | sed "s#\\\t#,#g"  | awk -F "float" '{print $2}' | sed 's#^,##g' | sed 's#^ ##g' | sed 's#rigt#right#g' > list_last_col.txt