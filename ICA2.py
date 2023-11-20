#!/usr/bin/python3

#import necessary modules first
import subprocess, sys, os
import numpy as np
import pandas as pd

#specify user input information about query
print("\nPlease enter the necessary information about the data you want to process\n")
Prot_name = input("\nPlease enter the name of PROTEIN you want to process\n") #store the user-input protein name as a variable
print("\nOK, thank you\n")
Species_name = input("\nPlease enter the name of SPECIES you want to process\n") #store the user-input species name as another variable
print("\nOK, thank you\n")

#use the information defined by user to download relevant data
species = str(Species_name.upper()) + "[Organism]"
protein = str(Prot_name.upper()) + "[Protein Name]"
user_query = f'{species} AND {protein}' #make query sentence
print("\nNow, downloading the data of " + species + protein)
output = f'{Species_name}.fasta' #set up the name for download file
subprocess.call(f'esearch -db protein -query "{user_query}" | efetch -format fasta > {output}',shell=True) #download file from NCBI database
print("\nDownloading Completed\n")
print("\nDownlaod file has been saved as " + f'{output}')

#count the downloaded sequences quantity
with open(f'{output}') as count_open:
	seq_count = count_open.read().count('>')
print('\nThere are ' + f'{seq_count}' + ' sequences in the downloaded file\n')

#make sure the input sequences are less than 1000
if int(seq_count) > 1000:
	print("\nSorry\nFor saving the time, this program can only process sequences less than 1000\n")
	exit()

#use the downloaded user defined data to build taxonomic group
#before build taxonomic group, we need to analyse the data first by clustalo
print("\nNow, starting analysing your data by clustalo\nAnd your input file is " + f'{output}')
output_clustalo = f'{Species_name}_clustalo.fasta' #set up the output file name of clustalo
subprocess.call(f'clustalo -i {output} --full -o {output_clustalo} --threads=64',shell=True) #using clustalo for analysing
print("\nClustalo analysing finished\nAnd the output file from clustalo has been saved as " + f'{output_clustalo}')

#blastp
blastp_input = input("\nDo you want to use BLSATP for extra alignment?\n[Y/N]") #ask the user whether need extra BLASTP analysing, and depends that to decide whether use BLASTP or not
if blastp_input.upper() in ('Y', 'YES'):
	print("\nOK, now using BLASTP for sequence alignment\n")
	subprocess.call(f'makeblastdb -in {output} -dbtype prot -out {Species_name}_local_database',shell=True)
	subprocess.call(f'blastp -query {output} -db {Species_name}_local_database -out {Species_name}_blastp_output.txt',shell=True)
	print(f"\nBLASTP analysing finished\nYour output file has been saved as {Species_name}_blastp_output.txt")
elif blastp_input.upper() in ('N', 'NO'):
	print("\nOK, give up BLASTP analysing\n")
else:
	print("\nBecause of the invalid option input, BLASTP alignment abandoned automatically\n")

#guide tree creation
guide_tree = input("\nDo you want build a guide tree for your data processed by clustalo?\nThe multiple substitution correction method used here is Kimura\nAnd Neighbor-joining is used for guide tree creation\n[Y/N]") #ask the user whether build a guide tree or not
if guide_tree.upper() in ('Y', 'YES'):
	tree_matrix = f'{Species_name}_matrix.ph'
	subprocess.call(f'distmat -sequence {output_clustalo} -protmethod 2 -outfile {tree_matrix}',shell=True)
	tree_output = f'{Species_name}_tree.nj'
	subprocess.call(f'neighbor < {tree_matrix} > {tree_output}',shell=True)
	print("\nGuide tree creation completed\nThe results have been saved as " + tree_matrix + " and " + tree_output)
	print("\nYour result may need FigTree for visible reading\nAnd this program is NOT provided here")
elif guide_tree.upper() in ('N', 'NO'):
	print("\nOK, give up guide tree creating\n")
else:
	print("\nBecause of the invalid option input, guide tree creation abandoned automatically\n")

#after clustalo analyse, use the output file of clustalo to build taxonomic group. and here I used plotcon from EMBOSS
print("\nNow, starting taxonomic group buidling\nNow, using plotcon to build taxonomic group\n")
print("\nPlease choose 1-4 for the number of columns to average alignment quality\nThe larger this input number is, the smoother the plot will be\n")
subprocess.call(f'plotcon -sequences {output_clustalo} -graph png',shell=True)
cwd = os.getcwd()
plot_graph = os.listdir(cwd)
print("\nOK, taxonomic group building completed\n")
for png_graph in plot_graph:
	if png_graph.endswith(".png"):
		os.rename(png_graph, f'{Species_name}.png') #rename the file to make the result recognizable
print("\nThis is your PNG photo of plotcon result\nPlease click CLOSE after reading the photo for further motif scanning\n")
subprocess.call(f'display {Species_name}.png',shell=True) #display the photo on the screen

#then scan the data for motif information
print("\nNow, starting MOTIF scanning\n")
motif_input = {}
motif_key = None
motif_content = ''

print("\nstarting variables created\n")
#store the header and mainbody of fasta file separately in a dictionary
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
#set up some variables for later use
motif_output_filename = Species_name + '_motif_ana' + '.txt'
motif_output = open(motif_output_filename,'w')
motifs_count_filename = f'{Species_name}_motif_scan_result.txt'
motifs_count_file = open(f'{Species_name}_motif_scan_result.txt','w')
motifs_count_read = motifs_count_file.write('Motifs')
file_count = 0 #set count number for break the loop

#analyse the fasta file with patmatmotif
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
	#store the splited fasta sequences in the temporary files
	subprocess.call(f'patmatmotifs -sequence {motif_process} -full -outfile {motif_out_tem}',shell=True) #process the input file with patmatmotif
	motif_out_tem2 = open(motif_out_tem,'r')
	motif_out_tem3 = motif_out_tem2.read()
	motif_output = open(motif_output_filename,'a')
	motif_output.write('\n' + motif_out_tem3) #write the output result within one file
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

#analyse and calculate the motifs in the result and display them
motif_output2 = open(motif_output_filename,'r')
motif_output_read = motif_output2.readlines()
for motif_out_line in motif_output_read:
	if motif_out_line.startswith('Motif = '):
		motifs_count_read = motifs_count_file.write('\n' + motif_out_line)
motif_output2.close()
motifs_count_file.close()

print("\nMotif analysing has completed, the output has been saved as " + motif_output_filename)
print("\nHere is the motifs occurred in the data and their counts:\n")

df = pd.read_csv(motifs_count_filename, sep="\t")
motif_result_screen = df['Motifs'].value_counts()
print(motif_result_screen)


