#!/bin/bash

# Loop through all directories that match the subsample_run? string
for dir in ./subsample_run?
do
  # Get the directory name
  directory=${dir}
  echo "$directory"
 
  # Copy the params file to the directory
  cp ./params.yml $directory

  # Run sisana
  path_to_params=$directory/params.yml
  echo "$path_to_params"
  # sisana generate $path_to_params
done
