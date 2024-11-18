# visium_hd_processing
Package installation
-CPU: conda env create -f stardist-cpu-env.yml (tested and working on linux-ubuntu 22.04)
-GPU:
    - Windows: conda env create -f stardist-cpu-env-windows.yml (Working on this section)
    - Linux: conda env create -f stardist-gpu-env-linux.yml     (Working on this section)

Activate environment: conda activate stardist-cpu-env

Install kernel for linux jupyter: python -m ipykernel install --user --name=stardist-cpu-env --display-name "Python (stardist-cpu-env)"



