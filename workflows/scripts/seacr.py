from snakemake.shell import shell
import os

seacr_sh = os.path.expanduser("~/tools/SEACR/SEACR_1.3.sh")

# Inputs      
treatment_bg = snakemake.input.treatment_bg
control_bg = snakemake.input.get('control_bg', None)

# Outputs
_peaks = snakemake.output._peaks  

# Params
out_prefix = snakemake.params.out_prefix
threshold = snakemake.params.threshold
use_control = snakemake.params.use_control
mode = snakemake.params.mode  

# Log
log = snakemake.log[0]

if use_control and control_bg:
    control_flag = f"{control_bg}"
    norm_flag = "norm" 
else:
    control_flag = f"{threshold}"
    norm_flag = "non"  
 

cmd_ = f"""
bash {seacr_sh} {treatment_bg} {control_flag} \
    {norm_flag} {mode} {out_prefix} &> {log}
"""

shell(cmd_) 