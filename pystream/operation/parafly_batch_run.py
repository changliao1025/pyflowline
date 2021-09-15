from pyearth.toolbox.slurm.slurm_prepare_job_script_python import slurm_prepare_job_script_python

def parafly_batch_run():
    
    slurm_prepare_job_script_python( sDirectory_job, \
        sBasename_job, \
        sBasename_parafly, \
        sJob_name, \
        iWalltime, \
        nNode_in = 1, \
        nTask_in=40, \
        sQueue_in='slurm')
    return