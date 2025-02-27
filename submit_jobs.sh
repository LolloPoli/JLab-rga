#!/bin/bash

# Directory delle liste create con make_lists.sh
#config_dir="RGB_fall2019_torus+1_lists"
config_dir="RGA_fall2018_torus+1_lists"
# Parametri SLURM
time_limit="03:00:00"
memory="2000"

output_dir="jobs_list"

echo "Lanciando job per i gruppi in $config_dir"

# Per ogni file di gruppo nella cartella
for group_file in "$config_dir"/group_*.txt; do
    # Estrai il numero del gruppo
    group_number=$(basename "$group_file" | grep -oP '\d+')

    # Crea uno script di sottomissione SLURM
    job_script="$output_dir/job_group_${group_number}.sh"
    log_file="logs/log_group_${group_number}.txt"

    cat << EOF > "$job_script"
#!/bin/bash
#SBATCH --job-name=job_group_$group_number
#SBATCH -t $time_limit
#SBATCH --mem=$memory
#SBATCH --output=$log_file

module load clas12root  # Carica l'ambiente (se necessario)

# Passa il file di gruppo al tuo programma
clas12root '4processing.cpp("'$group_file'")' -l -b -q 
EOF

    # Sottometti il job a SLURM
    sbatch "$job_script"
done

echo "Tutti i job sono stati sottomessi."
