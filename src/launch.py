import os
import subprocess

# Get The current directory 
jobs_dir = os.getcwd()

# Get all job, that have .sh extention 
job_files = [f for f in os.listdir(jobs_dir) if f.endswith('.sh')]

# Launch every job 
for job_file in job_files:
    job_path = os.path.join(jobs_dir, job_file)
    print(f"Lancio il job: {job_file}")

    # Compatibility with Python 3.6
    result = subprocess.run(['sbatch', '--wckey=rop', '--requeue', job_path],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)  

    # For the output
    print(result.stdout)
    print(result.stderr)

print("Tutti i job sono stati lanciati.")