#!/usr/bin/python3

#import necessary modules first
import subprocess, sys, os

#specify user input information about query
print("\nPlease enter the necessary information about the data you want to process\n")

#be careful about the input order here, need an error trap later
Prot_name = input("\nPlease enter the name of PROTEIN you want to process\n")
print("\nOK, thank you\n")
Species_name = input("\nPlease enter the name of SPECIES you want to process\n")
print("\nOK, thank you\n")

#use the information defined by user to download relevant data
species = str(Species_name.upper()) + "[Organism]"
protein = str(Prot_name.upper()) + "[Protein Name]"
user_query = f'{species} AND {protein}'
print("\nNow, downloading the data of " + species + protein)
output = f'{Species_name}.fasta'
subprocess.call(f'esearch -db protein -query "{user_query}" | efetch -format fasta > {output}',shell=True)
print("\nDownloading Completed\n")
print("\nDownlaod file has been saved as " + f'{output}')

#count the downloaded sequences quantity
with open(f'{output}') as count_open:
	seq_count = count_open.read().count('>')
print('\nThere are ' + f'{seq_count}' + ' sequences in the downloaded file\n')

#use the downloaded user defined data to build taxonomic group
#before build taxonomic group, we need to analyse the data first by clustalo
print("\nNow, starting analysing your data by clustalo\nAnd your input file is " + f'{output}')
output_clustalo = f'{Species_name}_clustalo.fasta'
subprocess.call(f'clustalo -i {output} --full -o {output_clustalo} --threads=64',shell=True)
print("\nClustalo analysing finished\nAnd the output file from clustalo has been saved as " + f'{output_clustalo}')

#after clustalo analyse, use the output file of clustalo to build taxonomic group. and here I used plotcon from EMBOSS
print("\nNow, starting taxonomic group buidling\nNow, using plotcon to build taxonomic group\n")
print("\nPlease choose 1-4 for the number of columns to average alignment quality\nThe larger this input number is, the smoother the plot will be\n")
subprocess.call(f'plotcon -sequences {output_clustalo} -graph png',shell=True)
print("\nOK, taxonomic group building completed\n")

#then scan the data for motif information
print("\nNow, starting MOTIF scanning\n")
motif_input = {}
motif_key = None
motif_content = ''

print("\nstarting variables created\n")

with open(f'{output}') as motif_open:
	motif_read = motif_open.readlines()
	for line in motif_read:
		if line.startswith('>'):
			motif_key = line
			motif_input[motif_key] = ''
		else:
			motif_content = line
			if motif_key != None:
				motif_input[motif_key] += motif_content

print("\nfasta dictionary created\n")

motif_output_filename = Species_name + '_motif_ana' + '.txt'
motif_output = open(motif_output_filename,'w')
cwd = os.getcwd()
file_count = 0 #set count number for break the loop
for motif_key_input in motif_input.keys(): #the downloaded data have been split into a dictionary
	motif_output.write(f'{motif_key_input}') #creat a file for store the output from patmatmotifs
	file_count += 1 #count the processing procedure
	motif_prepare = open(f'{file_count}_process.fasta','w') #creat the temporary file as an input file
	motif_prepare.write(f'{motif_key_input}{motif_input[motif_key_input]}') #write the '>' line and main body into input file
	motif_prepare.close() #close the input file and wait for scanning
	motif_process = f'{file_count}_process.fasta' #make input file readable for patmatmotif
	motif_out_prepare = open(f'{file_count}_out_tem.fasta','w')
	motif_out_prepare.close()
	motif_out_tem = f'{file_count}_out_tem.fasta'
	subprocess.call(f'patmatmotifs -sequence {motif_process} -full -outfile {motif_out_tem}',shell=True) #process the input file with patmatmotif
	motif_out_tem2 = open(motif_out_tem,'r')
	motif_out_tem3 = motif_out_tem2.read()
	motif_output = open(motif_output_filename,'a')
	motif_output.write('\n' + motif_out_tem3)
	motif_out_tem2.close()
	cwd_files = os.listdir(cwd) #set up a variable for further usage
	for file_name in cwd_files: #because of the high quantity of input files, I plan to delete every input file after motif scanning
		if file_name.endswith("_process.fasta"):
			file_remove = os.path.join(cwd, file_name)
			os.remove(file_remove)
		if file_name.endswith("_out_tem.fasta"):
			file_tem_remove = os.path.join(cwd, file_name)
			os.remove(file_tem_remove)
	if file_count == len(motif_input): #show the processing and break the loop
		print("\nMotif scaning finished\n")
		break
	else:
		print(f'Now is analysing the No.{file_count} data')
motif_output.close()

print("\nMotif analysing has completed, the output has been saved as " + motif_output_filename)

