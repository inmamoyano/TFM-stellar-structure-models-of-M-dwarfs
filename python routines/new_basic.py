"""
Grid building routine creating a file tree, inlists and script
This version is a basic single-routine code adapted for M dwarfs

GMM gmirouh@iaa.csic.es     Jun24
"""
import pprint
import subprocess
import itertools
import numpy as np

################################################################################################
# DEFINE PARAMETERS 
# Here par_dic is a dictionary, you should add parameters individually
par_dic = {}

# Define parameters to build the file tree for the runs
# The code will create two folders: one for the useful output ("data_dir") 
# and one for the temporary subfolders in which we run each track ("tmp").
par_dic["src_dir"]    = "/home/gmm/grid_builder/"         #Path to the code
par_dic["data_dir"]   = "/data/scratch/gmm/Mdwarfs/"      #Path to the useful output folder
par_dic["executable"] = '{src_dir}/build/star'.format(**par_dic)  #Path to the executable file "star" from MESA
par_dic["tmp"]        = '{data_dir}/runs/'.format(**par_dic)      #Temporary folder for the caches and runs

# The code creates its own mesa source file ("mesarc") and sources it
# The total number of CPUs used (OMP_NUM_THREADS * number of runs) should be <96 on festelar1/38 
par_dic["mesarc"]     = "{tmp}/mesa_rc.sh".format(**par_dic)   #File to be sourced
par_dic["MESA_DIR"]   = "/data/scratch/gmm/mesa-r23.05.1/"     #Where MESA is installed (don't change)
par_dic["OMP_NUM_THREADS"] =  "4"                              #CPUs per run
par_dic["MESASDK_ROOT"] = "/data/scratch/gmm/mesasdk/"         #Where the SDK is installed (don't change)

# VARIABLE PARAMETERS RANGES
# Define the quantities to be varied and the values they should assume
par_dic['variables']   = ['mass_range','alpha_range','z_range']
par_dic['mass_range']  = [ round(mass, 4) for mass in np.arange(0.1, 0.6+1e-5, 0.05) ]
par_dic['alpha_range'] = [1.82, 1.6]
par_dic['z_range']     = [0.0045,0.0142,0.028]

#cmd = ['lscpu'] #, '|', 'grep', 'CPU(s):']#, "| head -1 | awk '{print $2}'"]
#result = subprocess.check_output(cmd)
#for item in str(result).split('\\n'):
#    if item.split()[0] == "CPU(s):":
#        cpus = item.split()[1]
#par_dic["JOB_BATCH_SIZE"] = int(int(cpus)/int(par_dic["OMP_NUM_THREADS"]))
################################################################################################
# INITIALIZE RUN

# Make the temp subfolder where we'll run MESA instances
setup = subprocess.Popen("mkdir -p {tmp}".format(**par_dic), shell=True, executable="/bin/bash")
setup.wait()

# Create the appropriate MESA source file from the parameters above
f = open(par_dic["mesarc"], "w")
f.write("""
export OMP_NUM_THREADS={OMP_NUM_THREADS}
export MESA_DIR={MESA_DIR}
export MESASDK_ROOT={MESASDK_ROOT}
source $MESASDK_ROOT/bin/mesasdk_init.sh""".format(**par_dic)
)
f.close()

# Source the MESA environment variable file and compile the executable 
# Create the command
commands = ["source {}".format(par_dic["mesarc"])]
commands.extend(["cd {}".format(par_dic["executable"].replace("star", "")),
    "./mk",
    "cd -",
])

# Run the command 
setup = subprocess.Popen(";".join(commands), shell=True, executable="/bin/bash")
setup.wait()

################################################################################################
#SETUP RUN

# Initialize job counter
cnt_job = 0
# Create a script that goes in every folder and runs MESA
script = "#!/bin/bash"

# Create the sets of parameters for each run
sets   = list( itertools.product(par_dic['mass_range'], par_dic['alpha_range'], par_dic['z_range']) )

run_dic = {}
for element in sets:
    # For each set of variables, create an inlist and add a paragraph in the script file

    # Define run-specific parameters
    # This is not exactly automated
    run_dic["mass"]   = element[0]
    run_dic["alpha"]  = element[1]
    run_dic["z"]      = element[2]

    # Prepare all the subfolders
    run_dic["tmp"]    = par_dic["tmp"]
    run_dic["mesarc"] = par_dic["mesarc"]
    run_dic["tmp"]    = "{tmp}/tmp_mass{mass}_alpha{alpha}_z{z}".format(**run_dic)
    run_dic["output"] = par_dic["data_dir"]
    run_dic["output"] = "{output}/results/mass{mass}_alpha{alpha}_z{z}".format(**run_dic)
    run_dic["src_dir"]= par_dic["src_dir"]

    # Make tmp and output folders
    setup = subprocess.Popen("mkdir -p {tmp}; mkdir -p {output}".format(**run_dic), shell=True, executable="/bin/bash")
    setup.wait()

    # Adjust composition
    run_dic["he3"] = 2.0e-5               # Jim's recommendation
    run_dic["h2"]  = 2.3e-5               # Jim's recommendation
    run_dic["y"]   = 0.225 + (0.2703-0.255)*run_dic["z"]/0.0142            # Formula for Y
    run_dic["he4"] = run_dic["y"] - run_dic["he3"]
    run_dic["h1"]  = 1. - run_dic["h2"] - run_dic["he3"] - run_dic["he4"]  # X is what's left

    # Prepare the inlist file
    run_dic["inlist"]="""
! Auto-generated M dwarf inlist for Inma Moyano's JAE project

&star_job
! FILE MGMT
    history_columns_file = '{src_dir}/src/hist.list'
    profile_columns_file = '{src_dir}/src/prof.list'
    eosDT_cache_dir = '{tmp}/caches/eosDT_cache'
    kap_cache_dir   = '{tmp}/caches/kap_cache'
    rates_cache_dir = '{tmp}/caches/rates_cache'

! PMS
    create_pre_main_sequence_model = .true.
    save_model_when_terminate = .true.
    save_model_filename = '{output}/final.mod'
    show_log_description_at_start = .false.
    create_pre_main_sequence_model = .true.
    pre_ms_T_c = 5d5  ! default is 3d5, but MIST uses 5d5
    pre_ms_relax_to_start_radiative_core = .false.
    pre_ms_relax_num_steps = 300        !to be tested more

! NETWORK, RATE
    change_net = .true. ! switch nuclear reaction network
    new_net_name = 'pp_extras.net' ! 
    show_net_reactions_info = .false.
    show_net_species_info = .false.
      
! CHANGE SOLAR ABUNDANCES      
    set_uniform_initial_composition = .true.
    initial_h1  = {h1}
    initial_h2  = {h2}
    initial_he3 = {he3}
    initial_he4 = {he4}
    initial_zfracs = 6  ! AGSS09_zfracs
        
    pgstar_flag = .false.
/ !end of star_job

&eos
/ ! end of eos namelist

&kap
! OPACITY
    Zbase = {z}
    kap_file_prefix = 'a09'    ! default was 'gs98'
    kap_CO_prefix = 'a09_co'   ! default was 'gs98_co'
    kap_lowT_prefix = 'lowT_fa05_a09p'      ! Ferguson2005 for a09
/ ! end of kap namelist

&controls
! OUTPUT OPTIONS
    log_directory = '{output}/'
    photo_directory = '{tmp}/photos'
    extra_terminal_output_file = '{output}/term.out'
    photo_interval = 10000000
    profile_interval = 50 
    max_num_profile_models = -1
    terminal_interval = 20
    star_history_name = 'history.data'
    history_interval = 1   

! STARTING SPECIFICATIONS
    initial_mass = {mass} ! in Msun units
    initial_Z = {z}
    initial_Y = {y}

! MESH AND TIMESTEP PARAMETERS
    ! composition controls
    mesh_delta_coeff = 1.0    ! 0.5 to double the number of mesh points. default is 1.0
    max_dq = 1d-3             ! recommendation from mesa mailist
    time_delta_coeff = 0.10   ! 0.01 to increase timestep dramatically
      
    ! timestep limit
    min_timestep_limit = 1d-6        
      
    ! better resolution of the Henyey hook
    delta_lg_XH_cntr_max = -1        ! default was 0
      
! WHEN TO STOP
    max_age = 1.3787d10  !Age of Universe (Planck 2018) 
      
! PHYSICS      
    use_Ledoux_criterion = .true.
    alpha_semiconvection = 0.1
    thermohaline_coeff = 666.0
      
    ! mixing parameters
    mixing_length_alpha = {alpha}    ! sun has 1.82, Baraffe uses 1.6 
    MLT_option = 'Henyey'
      
    ! winds
    mass_change = 0d0     
   
    ! atmosphere options
    atm_option = 'table'    
    atm_table = 'photosphere'
/

&pgstar
/ ! end of pgstar namelist
    """.format(**run_dic)

    # Write the inlist file in the right subfolder
    f_i = open("{output}/inlist".format(**run_dic), "w") 
    f_i.write(run_dic["inlist"])
    f_i.close()

    # Define niceness
#    run_dic["nice"] = int(cnt_job/par_dic["JOB_BATCH_SIZE"])
    
    # Add appropriate lines to script text
    script +="""
job{}() {{""".format(cnt_job)
    script +="""
    source {mesarc}
    ulimit -v 5191680
    cd {tmp}
    export MESA_INLIST={output}/inlist
    /home/gmm/grid_builder//build/star > out 2>err &
}}
""".format(**run_dic)

    # Increase the job counter
    cnt_job +=1

#Add list of all job functions defined above
script +="""jobs=(\n"""
for i in range(cnt_job):
    script+="""job{}\n""".format(i)
script+=""")\n"""
script+="""
# Function to get the number of running jobs
running_jobs() {
    jobs=$(ps -e -o state=,pid=,comm= | grep "star" | wc -l)
    echo $jobs
}

# Maximum number of parallel jobs allowed
max_parallel_jobs=24

# Loop through the jobs array and run each job function with control
for job in "${jobs[@]}"; do
    while (( $(running_jobs) >= max_parallel_jobs )); do
        echo "Too many jobs running, waiting..."
        sleep 60
    done
    ($job) &
    sleep 3
    echo "Started $job with memory limit"
done

# Wait for all background jobs to finish
wait
echo "All jobs have been submitted."
"""

# Write script file in temp folder
# This script file normally can be run from anywhere
f_s = open("{tmp}/script".format(**par_dic), "w")
f_s.write(script)
f_s.close()

# Change permissions to make script file executable
setup = subprocess.Popen("chmod 777 {tmp}/script".format(**par_dic), shell=True, executable="/bin/bash")
setup.wait()
