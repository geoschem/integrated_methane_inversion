"""
Intended to configure Slurm by editing slurm.conf with instance-specific info. 
Should output to new_slurm.conf something similar to the following:

NodeName=ip-172-31-6-14 CPUs=36 RealMemory=94432 CoresPerSocket=36 ThreadsPerCore=1 State=UNKNOWN
PartitionName=debug Nodes=ip-172-31-6-14 Default=YES MaxTime=INFINITE State=UP
'ControlMachine=ip-172-31-6-14'
"""

import subprocess

slurm_info = (
    subprocess.run(["slurmd", "-C"], stdout=subprocess.PIPE)
    .stdout.decode("utf-8")
    .split()
)
first_line = "ControlMachine" + slurm_info[0][8:] + "\n"
second_line = (
    " ".join(
        [
            "\n",
            slurm_info[0],
            slurm_info[1],
            slurm_info[6],
            "SocketsPerBoard" + slurm_info[1][4:],
            "ThreadsPerCore=1 State=UNKNOWN",
        ]
    )
    + "\n"
)
third_line = (
    " ".join(
        [
            "PartitionName=debug",
            "Nodes" + slurm_info[0][8:],
            "Default=YES",
            "MaxTime=INFINITE",
            "State=UP",
        ]
    )
    + "\n"
)

with open(
    "/home/ubuntu/install-scripts/new_slurm.conf", "a+"
) as f:
    f.write("\n")
    f.write(first_line)
    f.write(second_line)
    f.write(third_line)
