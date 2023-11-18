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
print("Downlaod file has been saved as " + f'{output}')

#count the downloaded sequences quantity
with open(f'{output}') as count_open:
	seq_count = count_open.read().count('>')
print('There are ' + f'{seq_count}' + ' sequences in the downloaded file')

#use the downloaded user defined data to build taxonomic group
#before build taxonomic group, we need to analyse the data first by clustalo
print("Now, starting analysing your data by clustalo. And your input file is " + f'{output}')
output_clustalo = f'{Species_name}_clustalo.fasta'
subprocess.call(f'clustalo -i {output} --full -o {output_clustalo}',shell=True)
print("Clustalo analysing finished. And the output file from clustalo has been saved as " + f'{output_clustalo}')

#after clustalo analyse, use the output file of clustalo to build taxonomic group. and here I used plotcon from EMBOSS
print("Now, starting taxonomic group buidling.\nNow, using plotcon to build taxonomic group")
print("Please choose 1-4 for the number of columns to average alignment quality\nThe larger this input number is, the smoother the plot will be")
subprocess.call(f'plotcon -sequences {output_clustalo} -graph png',shell=True)
print("OK, taxonomic group building completed.")
