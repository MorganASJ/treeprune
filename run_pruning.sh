#!/bin/bash

# Iterate over the directories starting with "C_C"
for directory in C_C*; do
    # Check if the current item is a directory
    if [ -d "$directory" ]; then
        # Change to the directory
        cd "$directory"

        # Execute the Python script
        echo "-------------------- Found $directory: Running python script --------------------"
        python3 ../prune_trees_script.py

        echo ""
        echo "-------------------- Python script has been run --------------------"
        # Return to the original directory
        cd ..
    fi
done


# Iterate over the directories starting with "C_C"
for directory in C_C*; do
    # Check if the current item is a directory
    if [ -d "$directory" ]; then
        # Change to the directory
        cd "$directory"


    for file in *.treels; do
    if [ -f "$file" ]; then
        # Run the sed command on the matching file
        sed -i "s/'//g" $file
    fi
    done

    for file in *.pruned.tre; do
    if [ -f "$file" ]; then
        # Run the sed command on the matching file
        sed -i "s/'//g" $file
    fi
    done

        # Execute the bash script
        echo "-------------------- Found $directory: Running AU test --------------------"
        bash ../au_test.sh

        echo ""
        echo "-------------------- AU test completed --------------------"
        # Return to the original directory
        cd ..
    fi
done