import os, stat
from pathlib import Path
def _create_hpc_job(self, sSlurm_in=None):
    """create a HPC job for this simulation
    """
    os.chdir(self.sWorkspace_output)
    
    #part 1 python script         
    
    sFilename_pyflowline = os.path.join(str(Path(self.sWorkspace_output)) , "run_pyflowline.py" )
    ofs_pyflowline = open(sFilename_pyflowline, 'w')
       
    sLine = '#!/qfs/people/liao313/.conda/envs/pyflowline/bin/' + 'python3' + '\n' 
    ofs_pyflowline.write(sLine) 
    sLine = 'from pyflowline.pyflowline_read_model_configuration_file import pyflowline_read_model_configuration_file' + '\n'
    ofs_pyflowline.write(sLine)         
    sLine = 'sFilename_configuration_in = ' + '"' + self.sFilename_model_configuration + '"\n'
    ofs_pyflowline.write(sLine)
    sLine = 'oPyflowline = pyflowline_read_model_configuration_file(sFilename_configuration_in,'  \
        + 'iCase_index_in='+ str(self.iCase_index) + ',' \
        + 'dResolution_meter_in=' + "{:0f}".format(self.dResolution_meter)+ ',' \
        +  'sDate_in="'+ str(self.sDate) + '",' \
        +  'sMesh_type_in="'+ str(self.sMesh_type) +'"' \
        + ')'  +   '\n'   
    ofs_pyflowline.write(sLine)
    #if self.iFlag_flowline==1:
    #    sLine = 'oPyflowline.aBasin[0].dLatitude_outlet_degree=' \
    #        +  "{:0f}".format(self.aBasin[0].dLatitude_outlet_degree)+ '\n'   
    #    ofs_pyflowline.write(sLine)
    #    sLine = 'oPyflowline.aBasin[0].dLongitude_outlet_degree=' \
    #        + "{:0f}".format(self.aBasin[0].dLongitude_outlet_degree)+ '\n'   
    #    ofs_pyflowline.write(sLine)        
    sLine = 'oPyflowline.setup()' + '\n'   
    ofs_pyflowline.write(sLine)
    
    if self.iFlag_flowline ==1:                        
        sLine = 'oPyflowline.flowline_simplification()' + '\n'   
        ofs_pyflowline.write(sLine)     
    else:
        pass    
    sLine = 'aCell = oPyflowline.mesh_generation()' + '\n'   
    ofs_pyflowline.write(sLine)      
    sLine = 'oPyflowline.reconstruct_topological_relationship(aCell)' + '\n'   
    ofs_pyflowline.write(sLine)   

    sLine = 'oPyflowline.analyze()' + '\n'   
    ofs_pyflowline.write(sLine)      
    
    #sLine = 'oPyflowline.evaluate()' + '\n'   
    #ofs_pyflowline.write(sLine) 
    
    sLine = 'oPyflowline.export()' + '\n'   
    ofs_pyflowline.write(sLine)
    ofs_pyflowline.close()
    os.chmod(sFilename_pyflowline, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)          
    
    #part 2 bash script    
    sFilename_job = os.path.join(str(Path(self.sWorkspace_output)  ) ,  "submit.job" )
    ofs = open(sFilename_job, 'w')
    sLine = '#!/bin/bash\n'
    ofs.write(sLine)
    sLine = '#SBATCH -A ESMD\n'
    ofs.write(sLine)
    sLine = '#SBATCH --job-name=' + self.sCase + '\n'
    ofs.write(sLine)
    sLine = '#SBATCH -t 1:00:00' + '\n'
    ofs.write(sLine)
    sLine = '#SBATCH --nodes=1' + '\n'
    ofs.write(sLine)
    sLine = '#SBATCH --ntasks-per-node=1' + '\n'
    ofs.write(sLine)
    if sSlurm_in is not None:
        sSlurm = sSlurm_in
    else:
        sSlurm = 'slurm'
    sLine = '#SBATCH --partition='+ sSlurm + '\n'
    ofs.write(sLine)
    sLine = '#SBATCH -o stdout.out\n'
    ofs.write(sLine)
    sLine = '#SBATCH -e stderr.err\n'
    ofs.write(sLine)    
    sLine = 'module purge\n'
    ofs.write(sLine)
    sLine = 'module load gcc/8.1.0' + '\n'
    ofs.write(sLine)
    
    if self.iFlag_dggrid ==1:
        sLine = 'module load gdal/2.3.1' + '\n'
        ofs.write(sLine)

    sLine = 'module load anaconda3/2019.03' + '\n'
    ofs.write(sLine)
    sLine = 'source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh' + '\n'
    ofs.write(sLine)    
    sLine = 'conda activate hexwatershed' + '\n'
    ofs.write(sLine)
    sLine = 'cd $SLURM_SUBMIT_DIR\n'
    ofs.write(sLine)
    sLine = 'JOB_DIRECTORY='+ self.sWorkspace_output +  '\n'
    ofs.write(sLine)
    sLine = 'cd $JOB_DIRECTORY' +  '\n'
    ofs.write(sLine)
    sLine = 'python3 run_pyflowline.py' +  '\n'
    ofs.write(sLine)
    sLine = 'conda deactivate' + '\n'
    ofs.write(sLine)
    ofs.close()
    return