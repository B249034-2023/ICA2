#!/usr/bin/python3

#import necessary modules first
import subprocess, sys

#specify user input information about query
print("Please enter the necessary information about the data you want to process")
#be careful about the input order here, need an error trap later
Prot_name = input("Please enter the name of PROTEIN you want to process\n")

print("OK, thank you")
Species_name = input("Please enter the name of SPECIES you want to process\n")

print("OK, thank you")

#use the information defined by user to download relevant data
species = str(Species_name.upper()) + "[organism]"
protein = str(Prot_name.upper())
user_query = f'{species} AND {protein}'
print("Now, downloading the data of " + species + protein)
output = f'{Species_name}.fasta'
subprocess.call(f'esearch -db protein -query "{user_query}" | efetch -format fasta > {output}',shell=True)
print("Downloading Completed")
print("Downlaod file has been saved as " + output)
