import pybedtools
#pybedtools.set_bedtools_path('/ref/rmlab/software/spack/opt/spack/linux-rocky8-x86_64/gcc-8.5.0')
a = pybedtools.BedTool('/scratch/rmlab/rmlab_shared/ihreiss/Bmal1_CCs/032924_macs/041624_spikein/output_and_analysis/peaks/whole_brain_acsa2_bound_bound_cccaller_peaks.bed')

import subprocess
subprocess.check_output(['bedtools', '--version'], text=True).strip().split()[1]

c = a.cat(a)
