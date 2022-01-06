#!/usr/bin/python3

import os
import sys
import json
from pprint import pprint
from warnings import filterwarnings
import multiprocessing as mp
import itertools
import time
from datetime import datetime
import shutil

filterwarnings("ignore")

# try to import libraries, but install them if they
# are not installed already
try:
	from tqdm import tqdm
except ModuleNotFoundError:
	os.system('python3 -m pip install tqdm')
	from tqdm import tqdm

try:
	import requests
except ModuleNotFoundError:
	os.system('python3 -m pip install requests')
	import requests

try:
	from rdkit import Chem, DataStructs
except ModuleNotFoundError:
	os.system('python3 -m pip install rdkit-pypi')
	from rdkit import Chem, DataStructs

try:
	from sklearn.cluster import DBSCAN
except ModuleNotFoundError:
	os.system('python3 -m pip install sklearn')
	from sklearn.cluster import DBSCAN

try:
    import pandas as pd
except ModuleNotFoundError:
	os.system('python3 -m pip install pandas')
	import pandas as pd

try:
	import numpy as np
except ModuleNotFoundError:
	os.system('python3 -m pip install numpy')
	import numpy as np


# list for yes as user input to query
yes = ['y','Y']

# amino acids to distinguish hetatms/ligands in pdb structures
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'GLX', 'GLI', 'NLE', 'CYC']

# set default thresholds for program to operate
sequence_count_threshold = 1000

ligand_clustering_strictness = 0.8

evalue_cutoff = 0.01

identity_cutoff = 0.9

def parse_args(args): # function to parse command line user inputs

	args_dict = dict()

	# threads to use, defaults to 1 if not supplied
	args_dict["threads"] = [args[args.index('-threads') + 1] if '-threads' in args else '1'][0]

	# organism search term with spaces replaced to create folder and filenames
	args_dict["organism"] = args[args.index('-organism') + 1].replace(' ','_')

	# protein search term with spaces replaced to create folder and filenames
	args_dict["protein"] = args[args.index('-protein') + 1].replace(' ',"_")

	# boolean values for excluding poor quality sequences
	# these default to false and are true if supplied as arguments
	args_dict["partial"] = [True if "-partial" in args else False]

	args_dict["subdomain"] = [True if "-subdomains" in args else False]

	args_dict["isoform"] = [True if "-isoforms" in args else False]

	# working folder name and file name are based on search terms joined by underscore
	args_dict["foldername"] = f'{args_dict.get("organism")}_{args_dict.get("protein")}'
	args_dict["filename"] = f'{args_dict.get("foldername")}/{args_dict.get("organism")}_{args_dict.get("protein")}'
	args_dict["foldername"] = args_dict["foldername"].replace("*","glob")
	args_dict["filename"] = args_dict["filename"].replace("*","glob")

	return args_dict

def define_fetch_method(searchkey, args_dict): # function to determine correct esearch and efetch commands based on number of
											   # sequences and user specified input

	# get initial esearch result and number of sequences
	result = os.popen(searchkey).read()
	num_sequences = int(result.split('</Count>')[0].split('>')[-1])

	print(f'Found {num_sequences} relevant sequences...')

	# if there are more sequences than specfied threshold, get user input on solution
	if num_sequences > sequence_count_threshold:

		# define possible sorting methods for the esearch results
		sort_method_dict = {
							'1':'accession',
							'2':'date modified',
							'3':'date released',
							'4':'default order',
							'5':'organism name',
							'6':'taxonomy id'
							}

		# notify user and ask for input on solution
		print('**********************')
		print(f'WARNING: Large number of sequences detected! ({num_sequences})')
		print(f'It is recommended not to exceed {sequence_count_threshold} sequences in a single analysis')
		print(f'Available filtering methods: \n')

		# present the possible filtering methods to reduce the number of sequences
		print(f'***COMPLEX METHODS***')
		print(f'0 - Equal number of sequences from each species will be fetched (recommended)')
		print('\n')
		print(f'***SIMPLE METHODS***')
		print(f'Only the first {sequence_count_threshold} sequences will be fetched - please pick the sort method:')
		print('1 - Accession')
		print('2 - Date Modified')
		print('3 - Date Released')
		print('4 - Default Order')
		print('5 - Organism Name')
		print('6 - Taxonomy ID')
		print('\n')
		print('Or type EXIT to end the analysis')

		# get the users chosen method
		sort_method = input('Method choice: ')

		# exit if input is EXIT or exit
		if sort_method.upper() == 'EXIT':
			sys.exit()

		# otherwise implement method
		else:

			# exit program if chosen sort method is invalid or not in dictionary
			if sort_method is None:
				print('Invalid sort method - search aborted')
				sys.exit()

			# if sort method is simple then fetch the first x results sorted by
			# method of users choosing
			elif int(sort_method) > 0:

				# get the search method as a string to pass to eseatrch
				sort_method = sort_method_dict.get(sort_method)

				# add the sort argument onto the existing default esearch searchkey
				print(f'Fetching first {sequence_count_threshold} sequences sorted by {sort_method}...')
				sorted_searchkey = searchkey + f' -sort "{sort_method}"'

				# define the full command to send to the shell to fetch the first x sorted sequences
				full_query_command = f'echo "{result}" | efetch -format fasta -stop {sequence_count_threshold} > {args_dict.get("filename")}.prot.fa'

			# if sort method is complex then get as many sequences as possible from each
			# unique species returned in the search
			else:

				# get the taxonomy ids and accession numbers of all the returned sequences
				# and write them to a tsv file
				print(f'Using complex method 0 to fetch sequences...')
				varied_searchkey = searchkey + f"| efetch -format docsum | xtract -pattern DocumentSummary -element TaxId AccessionVersion > {args_dict.get('filename')}_varied_results.tsv"
				os.system(varied_searchkey)

				# read the produced tsv as a dataframe
				accession_df = pd.read_csv(f'{args_dict.get("filename")}_varied_results.tsv' ,sep="\t", names=['TaxID','AccID'])

				# count the unique species
				unique_species = len(accession_df['TaxID'].unique())

				# determine how many sequences per species can be used while staying below total sequence threshold
				sequences_per_species = int(sequence_count_threshold/unique_species)

				# if there are more species than the sequence threshold, then just truncate the list of accession ids
				if sequences_per_species == 0:
					print(f'Number of unique species is greater than {sequence_count_threshold}')
					print(f'Only fetching first {sequence_count_threshold} species in default order...')

					# drop duplicate species so one accession number/sequence per species remains
					query_df = accession_df.drop_duplicates(subset='TaxID')

					# define list of accession ids to fetch
					query_accession_ids = list(query_df['AccID'])[:sequence_count_threshold]

				# if there are less species than the sequence threshold, get accession ids of sequences to keep
				# based on the sequences_per_species variable
				else:
					print(f'Fetching {sequences_per_species} sequences per species...')

					# add count column for sequence per species
					accession_df['num_sequences'] = accession_df.groupby('TaxID').cumcount() + 1

					# subset the dataframe of accession ids to entries where count column
					# is less than sequences per species variable
					query_df = accession_df.loc[accession_df['num_sequences'] <= sequences_per_species]

					# define list of accession ids to fetch
					query_accession_ids = list(query_df['AccID'])

				# remove temporary esearch file
				os.system(f'rm -f {args_dict.get("filename")}_varied_results.tsv')

				# define the full command to send to the shell to fetch x sequences per species
				full_query_command = f'efetch  -db protein -format fasta -id {" ".join(query_accession_ids)} > {args_dict.get("filename")}.prot.fa'

	else:

		# if there are less sequences than the threshold, fetch them all
		full_query_command = f'echo "{result}" | efetch -format fasta > {args_dict.get("filename")}.prot.fa'

	return full_query_command





def count_motifs(filepath): # parse the motifs in a patmatmotifs output file into a list

	# open the patmatmotifs results file
    output = open(filepath,"r")

	# empty list to populate with motifs
    motifs = list()

	# if line contains motif count, then add name of motif to list
    for line in output:
            if 'Motif = ' in line:
                    motifs.append(line.split('Motif = ')[-1])

	# clean up and standardise list
    motifs = [m.strip() for m in motifs]

	# return the list of motifs if any were found
    if len(motifs) == 0:
            return None
    else:
            return motifs

def get_ligand_information(identifier): # gets the SMILE for a ligand based on its three letter PDB identifier

	# use requests to get html of ligand information page and parse
	# the html as a string
	ligand_info = str(requests.get(f'http://ligand-expo.rcsb.org/reports/{identifier[0]}/{identifier}/index.html').content)

	# if page exists, extract ligand SMILE string and return it
	try:
		ligand_smile = ''.join(ligand_info.split('<b>SMILES (CACTVS)</b></td> \\n<td valign=top align="left">')[1].split('</td>')[0].split('<wbr />'))
		return ligand_smile

	# if page does not exist, return none
	except:
		print(f'Error - could not get ligand information for {identifier}')
		return None


def get_pdb_ligands(raw_sequence_tuple): # query the PDB for structures with similar sequence to sequence argument
										 # then identify non amino acid residues in results and get their formula
										 # in SMILE format

    # unpack tuple supplied by the multiprocessing wrapper function
	sequence = raw_sequence_tuple[0]
	evalue_cutoff = raw_sequence_tuple[1]
	identity_cutoff = raw_sequence_tuple[2]

	# define pdb query parameters
	parameters = {
	"query": {
	        "type": "terminal",
	        "service": "sequence",
	        "parameters": {
	        "evalue_cutoff": evalue_cutoff,
	        "identity_cutoff": identity_cutoff,
	        "target": "pdb_protein_sequence",
	        "value":sequence
	        }
	},
	"return_type": 'mol_definition',
	"request_options": {
	"pager": {
	        "start": 0,
	        "rows": 100
	        },
	"scoring_strategy": "combined",
	"sort": [
	            {
	            "sort_by": "score",
	            "direction": "desc"
	            }
	            ]
	            }
	    }

	# submit post request to PDB query api
	response = requests.post('https://search.rcsb.org/rcsbsearch/v1/query', json=parameters).content

	# make empty dictionary to populate with and ligands found
	found_ligands = dict()

	# if no PDB results from sequence query, return empty dictionary
	if len(str(response)) == 3:
		return dict()

	# otherwise parse the response into a dictionary
	else:
		response = json.loads(response)

	# check if there is a problem with the query, and if so skip it by returning an empty dict
	if response.get('message') is not None and 'Invalid' in response.get('message'):
		return dict()

	# otherwise get chemical information for all non-amino acid residues
	else:
		# for all residues in all the returned structures
		for result in response['result_set']:

			# skip if residue is unknown or an amino acid
			if result['identifier'] in amino_acids or result['identifier'] == 'UNX':
				pass

			# otherwise use above function to get ligand chemical information in SMILE format
			else:
				ligand = get_ligand_information(result['identifier'])

				# if ligand information could not be obtained, then pass
				if ligand is None:
					pass

				# otherwise add ligand information and name to dictionary
				else:
					found_ligands[ligand] = result['identifier']

	# return dictionary of found ligands
	return found_ligands

def getTanimoto(mol1, mol2): # function to calculate similarity score between two ligands using Tanimoto Similarity and RDKit

	# turn ligands into RDKit molecule objects
    mol1 = Chem.MolFromSmiles(mol1)
    mol2 = Chem.MolFromSmiles(mol2)

	# calculate fingerprints from molecules
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)

	# return similarity score from fingerprint similarity
    return (DataStructs.TanimotoSimilarity(fp1,fp2))

def calculate_similarity_matrix(ligands): # calculate 2D similarity matrix of Tanimoto Similarity from list of unique ligands

	# turn list into numpy array
	lig_array = np.asarray(ligands)

	# create 2D pairwise similarity matrix for clustering algorithm
	similarity = np.array([[getTanimoto(l1,l2) for l1 in lig_array] for l2 in lig_array])

	return similarity

def cluster_ligands_and_cofactors(ligands, strictness): # clusters list of ligands into groups based on fingerprint similarity

	# calculate 2D similarity matrix of Tanimoto Similarity from list of unique ligands
	similarity_matrix = calculate_similarity_matrix(ligands)

	# define clustering algorithm
	affprop = DBSCAN(min_samples=1, eps=strictness)

	# fit the clustering algoritm to the similarity matrix
	affprop.fit(similarity_matrix)

	# set up list to create nested list of clusters
	cluster_members = list()

	# loop through clustering results and make list of clusters of ligands

	# for each cluster id
	for cluster_id in np.unique(affprop.labels_):

		# get list indexes of ligands that belong to this cluster
		cluster_indexes = np.nonzero(affprop.labels_==cluster_id)[0]

		# add the identified ligands to a list
		cluster = [ligand for index, ligand in enumerate(ligands) if index in cluster_indexes]

		# add this list to nested list above
		cluster_members.append(cluster)

	# return the nested list of ligands in clusters
	return cluster_members

def multiprocess_wrapper(function, items, threads): # generic function to multiprocess other functions

	# ensure processes cannot exceed number of cpus
	processes = min(threads, mp.cpu_count())

	# pass each item from items to the defined function in a multiprocessing pool
	# track progess with tqdm
	# fetch results with list wrapper of p.imap
	with mp.Pool(processes) as p:
		result = list(tqdm(p.imap(function, items), total=len(items)))
		p.close()
		p.join()

	# return list of results from defined function
	return result

def blastp_against_bindingdb(sequence_tuple): # function to blast a given protein sequence against a blast database of BindingDB protein targets

	# unpack tuple supplied as argument to individual variables
	seq = sequence_tuple[0]
	args_dict = sequence_tuple[1]
	evalue_cutoff = sequence_tuple[2]

	# define output filepaths and sequence filepath
	blast_dfs = f'{args_dict.get("foldername")}/blast_results/'
	seq_file = f'{args_dict.get("foldername")}/split_seqs/{seq}'

	# run blastp of the sequence file and write the results as a csv to the output filepath
	os.system(f'blastp -query {args_dict.get("foldername")}/split_seqs/{seq} -db {args_dict.get("foldername")}/bindingdb -evalue 0.01 -num_threads 1 -outfmt 10 > {blast_dfs}{seq}.csv')

def check_for_bindingdb_update(): # check if BindingDB reference files need updating

	print('Checking binding database is up to date...')

    # check program files folder exists
	if not os.path.isdir('program_files'):
		os.mkdir('program_files')

	# get current month and year
	currentMonth = int(datetime.now().month)
	currentYear = int(datetime.now().year)

	# define the url for BindingDB tsv file of all targets and ligands
	while True:
		update_url = f'https://www.bindingdb.org/bind/downloads/BindingDB_All_{currentYear}m{currentMonth}.tsv.zip'
		response = requests.head(update_url)
		if int(response.status_code) == 200:
			break
		else:
			if currentMonth == 1:
				currentYear = currentYear - 1
				currentMonth = 12
			else:
				currentMonth -= 1

	# check the checkpoint file created by this function to see if
	# file is already up to date
	# if it is then exit the function
	if os.path.isfile(f'program_files/{currentYear}m{currentMonth}.checkpoint'):
		print('BindingDB database already latest version!')
		return None

	# do not update unless a reason is found
	update = False

	# if there is no BindingDB tsv file of targets and ligands, fetch it
	if not os.path.isfile('program_files/BindingDB_All.tsv'):
		print('No binding database file is present - obtaining file...')

		# fetch and unzip the file
		os.system(f'wget {update_url}')
		os.system(f'unzip BindingDB_All_{currentYear}m{currentMonth}.tsv.zip')

		# only read in columns that need to be kept as file is very large
		df = pd.read_csv('BindingDB_All.tsv', sep="\t", usecols=['BindingDB MonomerID','Ligand SMILES','Ki (nM)','ZINC ID of Ligand','PubChem CID','PDB ID(s) of Target Chain'])

		# save BindingDB tsv to program files
		# keeping only subset of columns to save space
		df.to_csv('program_files/BindingDB_All.tsv', index=False)

		# remove redundant copies of file now that tsv has been processed
		os.system('rm -f BindingDB_All.tsv')
		os.system('rm -f BindingDB_All_{currentYear}m{currentMonth}.tsv.zip')

		# add a checkpoint file to log the update
		with open(f'program_files/{currentYear}m{currentMonth}.checkpoint','a+') as checkpoint:
			checkpoint.write('*******')
		checkpoint.close()

		# change update variable to true to update other files
		update = True

	# if there is a new update file available then ask user if they want to update
	elif requests.head(update_url).status_code == 200:
		answer = input('Update available! would you like to update bindingdb file? (y/N)')

		# if yes, update the file
		if answer == 'y' or answer == 'Y':

			# fetch and unzip the file
			os.system(f'wget {update_url}')
			os.system(f'unzip BindingDB_All_{currentYear}m{currentMonth}.tsv.zip')

			# only read in columns that need to be kept as file is very large
			df = pd.read_csv('BindingDB_All.tsv', sep="\t", usecols=['BindingDB MonomerID','Ligand SMILES','Ki (nM)','ZINC ID of Ligand','PubChem CID','PDB ID(s) of Target Chain'])

			# save BindingDB tsv to program files
			# keeping only subset of columns to save space
			df.to_csv('program_files/BindingDB_All.tsv', index=False)

			# remove redundant copies of file now that tsv has been processed
			os.system('rm -f BindingDB_All.tsv')
			os.system('rm -f BindingDB_All_{currentYear}m{currentMonth}.tsv.zip')

			# add a checkpoint file to log the update
			with open(f'program_files/{currentYear}m{currentMonth}.checkpoint','a+') as checkpoint:
				checkpoint.write('*******')
			checkpoint.close()

			# change update variable to true to update other files
			update = True

	# if there is no updated file available then user already has newest version
	else:
		print(requests.head(update_url).status_code)
		print('No bindingdb update available')

	# if there is no BindingDB file of protein target sequences then fetch it
	if not os.path.isfile(f'program_files/BindingDBTargetSequences.fasta') or update:
		os.system(f'wget https://www.bindingdb.org/bind/BindingDBTargetSequences.fasta')

		# move fetched file to program files
		shutil.move(f'BindingDBTargetSequences.fasta',f'program_files/BindingDBTargetSequences.fasta')

	# if BindingDB file of protein target sequences has not been formatted
	# or there is an updated version available, then update and format the
	# BindingDB file of protein target sequences
	if not os.path.isfile(f'program_files/formatted_bindingdb.fasta') or update:

		# open the unformatted fasta
		fasta = open(f'program_files/BindingDBTargetSequences.fasta','r')

		# remove the current formatted fasta if it exists
		if update:
			if os.path.isfile(f'program_files/formatted_bindingdb.fasta'):
			    os.system('rm -f program_files/formatted_bindingdb.fasta')

		# create empty file to write formatted fasta sequences to
		with open('program_files/formatted_bindingdb.fasta','a+') as newfasta:

			# for each line in unformatted fasta
		    for line in fasta:

					# remove non-ascii characters to prevent blast errors
		            line = line.encode("ascii","ignore")
		            line = line.decode()

					# write startlines of sequences with no changes
		            if '>' in line:
		                    newfasta.write(line)

					# if sequence lines are too long, split onto newlines of max 70 characters
					# and write to formatted file
		            elif len(line) > 70:
		                    total = len(line)
		                    newlines = int(total/70)+1
		                    for i in range(newlines):
		                            start = i*70
		                            end = start + 70
		                            print(newlines,start,end)
		                            if i == (newlines - 1):
		                                    newline = line[start:]
		                                    print(newline)
		                                    newfasta.write(f'{newline}\n')
		                            else:
		                                    newline = line[start:end]
		                                    print(newline)
		                                    newfasta.write(f'{newline}\n')
		# close new formatted fasta file
		newfasta.close()

	if os.path.isfile(f'BindingDB_All_{currentYear}m{currentMonth}.tsv.zip'):
		os.remove(f'BindingDB_All_{currentYear}m{currentMonth}.tsv.zip')

def main(): # run the program

	# parse user arguments
	args_dict = parse_args(sys.argv)

	# create directory for results of the specified search
	# based on user search query input
	if not os.path.exists(f'{args_dict.get("foldername")}'):
			os.mkdir(f'{args_dict.get("foldername")}')
			os.mkdir(f'{args_dict.get("foldername")}/split_seqs')
			os.mkdir(f'{args_dict.get("foldername")}/summary_files')

	# if identical search results folder exists, warn user
	# can either exit and backup the search results in different
	# location, or delete them and repeat search
	else:
			print('WARNING - data from an identical search exists already!')
			answer = input('Continuing will delete the existing data - continue? (y/N)')

			# delete existing search results folder and make new empty folers
			if answer in yes:
				os.system(f'rm -r -f {args_dict.get("foldername")}')
				os.mkdir(f'{args_dict.get("foldername")}')
				os.mkdir(f'{args_dict.get("foldername")}/split_seqs')
				os.mkdir(f'{args_dict.get("foldername")}/summary_files')

			# or exit search and let user backup data
			else:
				print('Please move or backup the existing data in a different location')
				sys.exit()

	# make program files folder if it does not exist already
	if not os.path.exists('program_files'):
	        os.mkdir('program_files')

	# build the esearch query based on user inputs
	searchkey = f'esearch -db protein -spell -query "{args_dict.get("organism").replace("_"," ")} [ORGN] AND {args_dict.get("protein").replace("_"," ")} [PROT]'

	# add exclusion of partial sequences, isoforms and subdomains unless user has asked to include them
	for key, value in args_dict.items():
		if value[0] is False:
			searchkey = searchkey + f' NOT {key.upper()}'
	searchkey = searchkey + '"'

	print(f'Searching for {args_dict.get("organism")} {args_dict.get("protein")} protein sequences...')

	# get correct search and fetch command based on total number of sequences initially returned
	# in case number exceeds sequence search threshold number
	fetch_command = define_fetch_method(searchkey, args_dict)

	# fetch the sequences
	os.system(fetch_command)

	print('Sequences fetched!')

	print('Performing MSA...')

	# use clustalo to do mutliple sequence alignment of sequences
	os.system(f'clustalo -i {args_dict.get("filename")}.prot.fa -o {args_dict.get("filename")}.aln --outfmt=clustal -v --force --threads {args_dict.get("threads")}')

	print('Sequences aligned!')
	print('\n')
	print('Checking conservation:')

	# use infoalign to analyse results of multiple sequence alignment
	os.system(f'infoalign {args_dict["filename"]}.aln -outfile {args_dict["filename"]}_alignment.info -nousa -noname -noweight -nodescription 2> /dev/null')

	# read in infoalign output
	alignment_info = pd.read_csv(f'{args_dict["filename"]}_alignment.info', sep="\t")

	# create row of mean counts as last row of dataframe
	alignment_info.loc['mean'] = alignment_info.mean().round(3)

	# display means of alignment statistics to user
	for col in alignment_info:

		# remove unwanted columns
		if 'Unnamed' in col:
			del alignment_info[col]

		# else display column name and mean
		else:
			print(f'Mean {col} : {alignment_info[col].mean().round(3)}')

	# write summary of alignment information as csv to summary_files folder
	alignment_info.to_csv(f'{args_dict["foldername"]}/summary_files/{args_dict["foldername"]}_alignment_info.csv')

	# create an svg of the conservation across residue numbers and save it in summary files folder
	os.system(f'plotcon {args_dict["filename"]}.aln -winsize 6 -graph svg -goutfile {args_dict["foldername"]}_conservation_plot -gdirectory {args_dict["foldername"]}/summary_files/ 2> /dev/null > /dev/null')

	# split fasta from efetch into separate sequences
	os.system(f'seqretsplit -sequence {args_dict.get("filename")}.prot.fa -outseq {args_dict.get("foldername")}/split_seqs/ -osdirectory2 {args_dict.get("foldername")}/split_seqs/ > /dev/null 2> /dev/null')

	# define filepaths for patmatmotifs input and output
	sequences = f'{args_dict.get("foldername")}/split_seqs/'
	motif_counts = f'{args_dict.get("foldername")}/motif_outputs/'

	print('\n')
	print('Counting PROSITE db motifs...')

	# make output directory for patmatmotifs output
	os.mkdir(motif_counts)

	# for each individual sequence file, count the motifs with patmatmotifs and save output to output directory
	for file in tqdm(os.listdir(sequences)):
	  os.system(f'patmatmotifs -full -sequence {args_dict.get("foldername")}/split_seqs/{file} -outfile {file.replace(".fasta",".motifs")} -rdirectory {motif_counts} 2> /dev/null')

    # count how many sequences there are in total from the search
	total_sequences = len(os.listdir(f'{args_dict.get("foldername")}/motif_outputs/'))

	# make empty lists and dictionaries for populating with motif information from loop below
	# total_motifs list will be a list of names of motifs, with one name for each individual occurence
	total_motifs = list()

	# motif_sources will be a dictionary of motif names as keys, and lists of sequence
	# files where that motif is present as values
	motif_sources = dict()

	# for each of the patmatmotif output files
	for file in os.listdir(f'{args_dict.get("foldername")}/motif_outputs/'):

		# get a list of the motifs that are present with the count_motifs() function
		motifs = count_motifs(f'{args_dict.get("foldername")}/motif_outputs/{file}')

		# add these to the overall list of motifs (total_motifs)
		if motifs is not None:
			total_motifs.extend(motifs)

			# then for each motif
			for motif in set(motifs):

				# if the motif is not in the motif_sources dictionary
				# then add it, and add its value as a list with the sequence
				# filename as a member
				if motif_sources.get(motif) is None:
					motif_sources[motif] = [file.replace('.motif','')]

				# if motif is already in the motif_sources dictionary
				# then add the sequence filename to the list stored as its value
				else:
					motif_source_list = motif_sources.get(motif)
					if file.replace('.motif','') not in motif_source_list:
						motif_source_list.append(file.replace('.motif',''))
						motif_sources[motif] = motif_source_list


	# count the overall occurences of each motif and store them in
	# motif_counts dictionary
	motif_counts = dict()
	for motif in set(total_motifs):
	     motif_counts[motif] = total_motifs.count(motif)

	# print the motif counts to screen
	print('\n')
	print('Total motif counts:')
	for motif, count in motif_counts.items():
		print(f'Motif {motif} found {count} times')

	print('\n')
	print('Relative motif frequencies: ')
	if len(motif_counts) != 0:

		# calculate relative frequencies of motifs based on total sequences
		# and raw counts and print these values to screen
		for motif, count in motif_counts.items():
			print(f'Motif: {motif} present {str(count/total_sequences)[:8]} times per sequence')

		# construct summary dataframe of motif count information
		motif_df = pd.DataFrame({'Motif':motif_counts.keys(),'Count':motif_counts.values()})

		# add count per sequence and source sequences
		motif_df['Count Per Sequence'] = motif_df['Count']/total_sequences
		motif_df['Found in Sequences'] = motif_df['Motif'].map(motif_sources)

		# add sequence search information columns
		motif_df['Sequences Queried'] = total_sequences
		motif_df['Taxonomy Search Key'] = args_dict.get('organism')
		motif_df['Protein Search Key'] = args_dict.get('protein')

		# save to summary files
		motif_df.to_csv(f'{args_dict.get("foldername")}/summary_files/{args_dict.get("foldername")}_PROSITE_motif_summary.csv', index=False)

	# if no motifs found then notify user
	else:
	      print(f'No known motifs found in requested sequences')

    ######### WILDCARD ANALYSIS #########

	print('\nQuerying BindingDB for known ligands...')

	# make a blast searchable database from the formatted BindingDB fasta file of protein targets
	os.system(f'makeblastdb -in program_files/formatted_bindingdb.fasta -dbtype prot -out {args_dict.get("foldername")}/bindingdb > /dev/null')

	# define folder names and make output directory for blast results
	sequences = f'{args_dict.get("foldername")}/split_seqs/'
	blast_dfs = f'{args_dict.get("foldername")}/blast_results/'
	os.mkdir(blast_dfs)

	# get filenames of all the individual search query sequence files
	allseqs = os.listdir(sequences)

	# build the sequence filename, arguments dictionary and evalue cutoff
	# into a tuple, and then pass this list of tuples to the multiprocess_wrapper
	# function to multiprocess blastp search of BindingDB blastdb
	allseqs_tuple = [(seq, args_dict, evalue_cutoff) for seq in allseqs]
	multiprocess_wrapper(blastp_against_bindingdb, allseqs_tuple, int(args_dict.get("threads")))

	print('Querying the Protein Data Bank for known ligands...')

	# make empty list to populate with strings of sequences
	raw_seqs = list()

	# for each sequence file, open it and extract the protein sequence as a string
	for seq in allseqs:

		# define filename
		seq_file = f'{args_dict.get("foldername")}/split_seqs/{seq}'

		# open the file and read it skipping the header line
		raw_seq_lines = open(seq_file,'r').readlines()[1:]

		# turn lines into one string without newline characters
		raw_seq = "".join([line.strip() for line in raw_seq_lines])

		# add the raw sequence string to list above
		raw_seqs.append(raw_seq)

	# build the sequence filename, evalue cutoff and identity cutoff
	# into a tuple, and then pass this list of tuples to the multiprocess_wrapper
	# function to multiprocess querying of protein data bank
	raw_seqs_tuple = [(raw_seq, evalue_cutoff, identity_cutoff) for raw_seq in raw_seqs]

	# multiprocess querying of pdb with sequences with threads as 2 to not exceed the query API limit
	pdb_ligands = multiprocess_wrapper(get_pdb_ligands, raw_seqs_tuple, 2)

	# convert pdb_ligands variable from a list of dictionaries to a single dictionary
	# with ligand smiles as keys and PDB identifier codes as values
	dict_pdb_ligands = dict()
	for subdict in pdb_ligands:
		dict_pdb_ligands.update(subdict)

	# define headers for blastp output results to be read in as a dataframe
	blast_headers = ['query acc.ver', 'subject acc.ver', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']

	# define empty list to add dataframes to for concatenation
	blast_result_dataframes = list()

	# open each blastp output csv file as a dataframe and add it to the blast_results_dataframes list
	for blast_result in tqdm(os.listdir(blast_dfs)):
	        blast_file = f'{blast_dfs}{blast_result}'
	        df = pd.read_csv(blast_file, names=blast_headers)
	        blast_result_dataframes.append(df)

	# concatenate blast results into one single dataframe
	blast_results_df = pd.concat(blast_result_dataframes)

	# keep only the hits which are above the identity cutoff value
	good_hits = blast_results_df.loc[blast_results_df['% identity'] > (identity_cutoff*100)]

    # remove leading 'p' from the BindingDB protein target IDs and convert them to integers so they can
	# be used to query the BindingDB ligands database
	good_hits['subject acc.ver'] = good_hits['subject acc.ver'].apply(lambda x: x.replace('p','')).astype(int)

	# load the BindingDB ligand database as a csv file
	ligand_db = pd.read_csv(f'program_files/BindingDB_All.tsv')

	# query the BindingDB ligand database file for only ligands that are known to bind targets that are present
	# in our blastp search results
	bound_ligand_db = ligand_db.loc[ligand_db['BindingDB MonomerID'].isin(list(good_hits['subject acc.ver']))]

	# get all the unique ligands identified from the BindingDB search
	binding_db_ligands = list(bound_ligand_db['Ligand SMILES'].unique())

	# add the identified PDB ligands to the list of BindingDB ligands to get list of all found ligands
	binding_db_ligands.extend(list(dict_pdb_ligands.keys()))

	# remove duplicate ligands from list of all found ligands
	unique_ligands = list()
	[unique_ligands.append(ligand.strip().lstrip()) for ligand in binding_db_ligands if ligand.strip().lstrip() not in unique_ligands]

	# make empty dictionary to populate with clustering results
	cluster_key = dict()

	# if there are any ligands found
	if len(unique_ligands) != 0:

		# cluster the ligands into groups based on Tanimoto Similarity using the DBSCAN algorithm
		cluster_members = cluster_ligands_and_cofactors(unique_ligands, ligand_clustering_strictness)

		# for each cluster identified in the nested list variable cluster_members
		for index, cluster in enumerate(cluster_members):

			# tell the user the cluster number
			print(f'Ligand Cluster {index+1}: ({len(cluster_members[index])} members)')

			# print the members of the cluster
			for member in cluster:
				cluster_key[member] = index+1
				print(f'{member}')
			print('*************')

		# make a summary dataframe of found ligands from the cluster_key dictionary
		ligand_summary_df = pd.DataFrame({'Ligand':cluster_key.keys(),'Cluster Group':cluster_key.values()})

		# add ZINC database ligand IDs for all ligands where they are present in the BindingDB ligand database file
		ligand_summary_df['ZINC Database ID'] = ligand_summary_df['Ligand'].map(dict(zip(ligand_db['Ligand SMILES'],ligand_db['ZINC ID of Ligand'])))

		# add PubChem ligand IDs for all ligands where they are present in the BindingDB ligand database file
		ligand_summary_df['PubChem ID'] = ligand_summary_df['Ligand'].map(dict(zip(ligand_db['Ligand SMILES'], ligand_db['PubChem CID'])))

		# add PDB residue identifier codes for all ligands that were sourced from the PDB
		ligand_summary_df['PDB Ligand ID'] = ligand_summary_df['Ligand'].map(dict_pdb_ligands)

		# add inhibition constant values for all ligands where the information is available from BindingDB ligand database file
		ligand_summary_df['Ki (nM)'] = ligand_summary_df['Ligand'].map(dict(zip(bound_ligand_db['Ligand SMILES'], bound_ligand_db['Ki (nM)'])))

		# add details of PDB structures of identified high identity BindingDB protein targets
		ligand_summary_df['Homologous PDB Structures'] = ligand_summary_df['Ligand'].map(dict(zip(bound_ligand_db['Ligand SMILES'], bound_ligand_db['PDB ID(s) of Target Chain'])))

		# save the ligand_summary_df to the summary_files folder
		ligand_summary_df.to_csv(f'{args_dict["foldername"]}/summary_files/{args_dict["foldername"]}_known_ligand_structural_summary.csv', index=False)

	# otherwise, if no ligands have been found in the PDB or BindingDB, notify the user
	else:
		print('No structures or ligands found!')

	# clean up the working files from the created search directory
	os.mkdir(f'{args_dict.get("foldername")}/working_files')
	for file in os.listdir(f'{args_dict.get("foldername")}'):
		if not os.path.isdir(f'{args_dict.get("foldername")}/{file}'):
			shutil.move(f'{args_dict.get("foldername")}/{file}',f'{args_dict.get("foldername")}/working_files/{file}')

	print(f'Analysis results have been saved to {args_dict["foldername"]}/summary_files/')

# run script when called from the command line
if __name__ == "__main__":
	check_for_bindingdb_update()
	main()
