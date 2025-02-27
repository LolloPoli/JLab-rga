#!/bin/bash

# Percorsi
input_file="data_lists/rga_fall2018_torus+1_pass2.txt"  # File con la lista di tutti i .hipo
output_dir="RGA_fall2018_torus+1_lists"  # Directory dove salvare i gruppi di liste
files_per_group=5  # Numero di file per ogni job

# Crea la cartella di output se non esiste
mkdir -p "$output_dir"

# Leggi i percorsi dei file dal file di input
mapfile -t file_list < <(tail -n +2 "$input_file")  # Legge la lista ignorando la prima riga

# Conta il numero totale di file
total_files=${#file_list[@]}
group_index=0

# Suddivide i file in gruppi da 10
for ((i=0; i<total_files; i+=files_per_group)); do
    group_file="$output_dir/group_${group_index}.txt"
    
    # Definisce il nome del file di output (puoi cambiarlo se vuoi)
    output_file="output_group_${group_index}.root"

    # Scrive il nome del file di output nella prima riga
    echo "/volatile/clas12/lpolizzi/sidis/rga/$output_file" > "$group_file"

    # Aggiunge i file .hipo
    for ((j=i; j<i+files_per_group && j<total_files; j++)); do
        echo "${file_list[j]}" >> "$group_file"
    done

    ((group_index++))
done

echo "Generati $(ls $output_dir | wc -l) gruppi di file nella cartella '$output_dir'."
