# ARGS: the list of nPCs

for pc in $npcs
do
  screen -dmS ridge_pc_$pc bash run_idp_pcs.screen $pc
done