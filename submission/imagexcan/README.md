Round 2.

```
# T1 IDP
screen -dmS t1_2nd bash -c 'bash run_round_2.screen t1'
# dMRI IDP
screen -dmS dmri_2nd bash -c 'bash run_round_2.screen dmri'
```

Round 3.

```
# T1 IDP
screen -dmS t1_3rd bash -c 'bash run_round_3.screen t1'
# dMRI IDP
screen -dmS dmri_3rd bash -c 'bash run_round_3.screen dmri'
```

Round 4.

The logfile was wrongly labeled as `round_3` ..

```
# T1 IDP
screen -dmS t1_no_pc_4th bash -c 'bash run_round_4.screen t1_no_pc'
screen -dmS t1_w_pc_4th bash -c 'bash run_round_4.screen t1_w_pc'

# dMRI IDP
screen -dmS dmri_no_pc_4th bash -c 'bash run_round_4.screen dmri_no_pc'
screen -dmS dmri_w_pc_4th bash -c 'bash run_round_4.screen dmri_w_pc'
```
