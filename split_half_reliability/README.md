NETWORKS 
1. Split half dense connectivity matrices were derived from resting-state data from the odd versus even sessions.
2. Those were fed as inputs to the Lynch PFM code (https://github.com/cjl2007/PFM-Depression) to generate networks for each half.
3. Dice overlap for relevant networks was compared using split_half_network_dice.m

TASKS
1. Split half task activation maps were derived using odd versus even task runs (see task_processing folder).
2. Correlation between split halves was compared using split_half_task_corr.m
