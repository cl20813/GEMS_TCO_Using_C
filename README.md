
### Important!
All the header files (.h), sources files including main.cpp should be located at the same level inside the "GEMS_C" folder.

# Copy folder
scp -r "C:\Users\joonw\source\repos\GEMS_C2" jl2815@amarel.rutgers.edu:/home/jl2815/tco            # r stands for recursive, used to copy a folder.

If face errors, upload files separately. 

scp "C:\Users\joonw\source\repos\GEMS_C\matern_cov.cpp" jl2815@amarel.rutgers.edu:/home/jl2815/tco 
scp "C:\Users\joonw\source\repos\GEMS_C\matern_cov.h" jl2815@amarel.rutgers.edu:/home/jl2815/tco
scp "C:\Users\joonw\source\repos\GEMS_C2\main.cpp" jl2815@amarel.rutgers.edu:/home/jl2815/tco

# Create the Slurm Batch Script
ssh jl2815@amarel.rutgers.edu   # Now log back to amarel
module use /projects/community/modulefiles  #without this, I can't load 2024.06-ts840

nano GEMS_C.sh                  # open a new text editor     
'''
#!/bin/bash          
#SBATCH --job-name=C_test_job        # Job name         
#SBATCH --output=/home/jl2815/tco/output_%j.out            # Standard output file (%j = JobID)        
#SBATCH --error=/home/jl2815/tco/error_%j.err              # Standard error file (%j = JobID)          
#SBATCH --time=01:00:00                   # Maximum time          
#SBATCH --ntasks=1                        # Number of tasks         
#SBATCH --cpus-per-task=8                 # Number of CPU cores per task         
#SBATCH --mem=16G                          # Memory per node          
#SBATCH --partition=main               # Partition to submit to           

### Load the Anaconda module to use srun 
      
module purge                                     # unload every other environment to avoid conflict        
module use /projects/community/modulefiles                  # without this, I can't load 2024.06-ts840          
module load gcc/11.2/openmpi/4.1.6-ez82                       

### Navigate to the directory with the executable
cd /home/jl2815/tco/GEMS_C

### Compile the code (if not done already)
g++ -o GEMS_C_executable main.cpp    #  -o GEMS_C tells the compiler to create an executable file called GEMS_C.      

### Run the executable
./GEMS_C_executable
'''

# Submit the Job using sbatch
sbatch GEMS_C.sh

# Srun job order example
cd /home/jl2815/tco/GEMS_C
g++ -o GEMS_C_executable main.cpp   # check  ls -l /home/jl2815/tco/GEMS_C
srun --ntasks=1 --cpus-per-task=1 --mem=8GB --time=01:00:00 ./GEMS_C_executable

squeue -u jl2815        # Status of my job
sinfo --Node --long     # View all available node
sstat -j 38879656       # Check whether the job has initialized 
scontrol show job 38880272  # Show job details
sprio -j <jobID>            #  job scheduling priority
scancel 38886671

# Check the output
cd ./tco           
cat output_<jobID>.out         #actual job id
cat error_38868490.err

# copy output files to my computer
scp jl2815@amarel.rutgers.edu:/home/jl2815/GEMS/output_38868525.out C:/Users/joonw/Downloads/
