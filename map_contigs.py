import argparse
import csv
from Bio import SeqIO

parser = argparse.ArgumentParser(description = "reads two .fasta files")
parser.add_argument("--fasta_1", help = ".fasta file 1", required = False)
parser.add_argument("--fasta_2", help = ".fasta file 2", required = False)
parser.add_argument("--csv_file",help = ".csv file with contig ids but without sequences", required = False)
args = parser.parse_args()

fasta_1 = args.fasta_1
fasta_2 = args.fasta_2
csv_file= args.csv_file
#contig dictionary with sequences as keys and contig names as values
def get_contig_dict(fasta_file):
        contig_dict = {}
        with open(fasta_file, 'rU') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                        contig_dict[str(record.seq)] = record.id
        return contig_dict

#contig dictionary with contig names as keys and sequences as values
def get_contig_dict_2(fasta_file):
	contig_dict = {}
	with open(fasta_file, 'rU') as handle:
		for record in SeqIO.parse(handle, 'fasta'):
			contig_dict[record.id] = str(record.seq)
	return contig_dict

def map_ids(dict1, dict2):
        mapped_ids = []
        differ_ids = 0
        for key in dict1.keys():
                if key in dict2.keys():
                        mapped_ids.append( [dict1[key], dict2[key] ] )
                        print mapped_ids
                        if dict1[key] != dict2[key]:
                                differ_ids += 1
                else:
                        mapped_ids.append( [dict1[key] , "NA" ] )
        if differ_ids  == 0:
                mapped_ids.insert( 0 , ['Different_Ids' , differ_ids] )
        return mapped_ids

def read_csv(csv_file):
	with open(csv_file ,'r') as f:
		reader = csv.reader(f)
		table = [ [ str(j) for j in i] for i in reader]
	return table

def get_new_name(old_name, new_string):
	new_name = str(old_name)
	new_name = new_name.strip().split('.')
	new_name = new_name[:6] + [new_string] + new_name[6:] 
	new_name = '.'.join(new_name)
	print new_name
	return new_name

def write_csv(new_name , table):
	with open(new_name , 'w+') as output_f:
		writer = csv.writer(output_f, quoting = csv.QUOTE_NONE, escapechar = ' ')
		for row in table:
			writer.writerow(row)

def get_tr_name(csv_file , fasta_file):
	contig_dict = get_contig_dict(fasta_file) #keys: sequence, values: transcript ID
	table = read_csv(csv_file)
	for row in table[1:]: #skip header
		if row[-1] in contig_dict:
			row.insert( 7, contig_dict[row[-1]] )
		else:
			row.insert(7, 'NA')
	table[0].insert(7, 'tr_id')
	new_name = get_new_name(csv_file , 'tr')
	write_csv(new_name , table)
#Appends sequence to a .csv file given contig name
def get_sequence(csv_file, fasta_file):

	contig_dict = get_contig_dict_2(fasta_file) 
	table = read_csv(csv_file)
	#write a table appending the sequences that have the corresponding contig name
	for row in table[1:]: #skip header
		if row[8] in contig_dict:
			row.append(contig_dict[row[8]])
		else:
			row.append('MISSING')
	table[0].append('sequence')
	new_name = get_new_name(str(csv_file) , 'seq')
	#write new csv file with sequences
	write_csv(new_name , table)


def main():
	#get_sequence(csv_file , fasta_1)
	get_tr_name(csv_file , fasta_1)
main()
