#!/bin/bash
 
#setting the number of jobs to be executed
#$ -t 1-2

#infile=`sed -n "$SGE_TASK_ID p" list.txt`
 
#echo $infile

#cp -r ./${infile}/clusterblast/*.txt ./textfiles/${infile}

# Loop through each line in list.txt
while IFS= read -r infile; do
    echo $infile
    echo "Copying files for directory: ${infile}"

    # Check if source directory exists
    if [ -d "./${infile}/clusterblast/" ]; then
        # Create destination directory if it doesn't exist
        mkdir -p "./textfiles/${infile}"

        # Copy *.txt files from source to destination
        cp -r "./${infile}/clusterblast/"*.txt "./textfiles/${infile}/"
        
        echo "Files copied successfully."
    else
        echo "Error: Source directory ./${infile}/clusterblast/ not found."
    fi

done < list.txt
 
