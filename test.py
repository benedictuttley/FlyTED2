import os
import json
import sys
from Bio.Seq import Seq
import csv
from termcolor import colored
from intermine.webservice import Service
import textwrap
import mysql.connector

final = []  # Will contain the final data sets --> [PROBE, TRANSCRIPT ID, 3', 5', TARGET SEQUENCE]
row_count = 0

# [1] TEST WITH DUMMY TRANSCRIPT SEQUENCE --> STATUS: PASSED
# [2] FETCH TRANSCRIPT FROM FLYMINE USING API --> STATUS: PASSED
# [3] PIPE DATA INTO Probe_BLAST_for_target_sequence.py --> STATUS: PASSED
# [4] CREATE Probe_BLAST_for_target_sequence.py FILE FOR TRANSCRIPT SEQUENCE --> STATUS: PASSED
# [5] PERFORM BLAST and output target sequence to text file --> STATUS: PASSED

# Find number of entries in the csv, two lines per gene, one for 3' and one for 5'

with open('/var/www/html/FlyTED2/'+ sys.argv[1]) as csvfile:

    csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    next(csvreader)
    rows = csvreader
    row_count = sum(1 for row in rows)

# Write eacg dataset for each gene on one line
with open('preprocessed.csv', 'w') as csvfile:
    fieldnames = ['PROBE', '3-END', '5-END']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    with open('/var/www/html/FlyTED2/' + sys.argv[1]) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(csvreader)
        rows = csvreader
        counter = 0
        for row in rows:
            probe_id = row[0].replace('-T3', '').replace('3_','')
            three_end = row[1]

            if (counter < (row_count / 2)):
                five_end = next(rows)[1]

            counter = counter + 1


            # CONSTANTS FOR TWO.csv
            # 3 - TAACCCTCACTAAAGGG
            # 5 - ATTTAGGTGACACTATAGAA

            # REMOVE CONSTANT SEQUENCE FROM THE 3' END:
            three_constant_type_A = 'TAACCCTCACTAAAGGG'
            three_constant_type_B = 'GTAATTAACCCTCACTAAAGGG'
            three_constant_type_C = 'CTAATTAACCCTCACTAAAGGG'
            three_constant_type_D = 'TAACCCTCACTAAAGGG'
            five_constant_type_A = 'ATTTAGGTGACACTATAGAA'
            three_end = three_end.replace(three_constant_type_A, '')
            five_end = five_end.replace(five_constant_type_A, '')
            writer.writerow({'PROBE': probe_id, '3-END': three_end, '5-END': five_end})


# READ PROBE SETS LINE BY LINE AND PERFORM BLASTN ON EACH PROBE AND THE TRANSCRIPT(S) FASTA FILES FOR A GIVEN GENE:
with open('preprocessed.csv') as csvfile:
    gene_probe_sets = csv.reader(csvfile)
    next(gene_probe_sets)
    for data_set in gene_probe_sets:

        # Write 3' probe to three.seq
        with open("three.seq", "w") as seq_file_three:
            three_end = data_set[1]
            seq = Seq(three_end)
            data_set[1] = str(seq.reverse_complement())
            seq_file_three.write(data_set[1])


        # Write 5' probe to five.seq
        with open("five.seq", "w") as seq_file_five:
            seq_file_five.write(data_set[2])


        # Need to fetch the transcript sequences for the respective gene from FlyMine:
        service = Service("http://www.flymine.org/flymine/service")
        query = service.new_query("Gene")
        query.add_view("transcripts.primaryIdentifier", "transcripts.length")
        query.add_sort_order("transcripts.length", "DESC")
        query.add_constraint("Gene", "LOOKUP", data_set[0], code="A")
        transcript_primary_identifiers = []

        for row in query.rows():

            transcript_primary_identifiers.append(row["transcripts.primaryIdentifier"])



        # Need to fetch the associated sequence for each of the gene transcript ID's that were returned from the above query:
        for transcript in transcript_primary_identifiers:


            service = Service("http://www.flymine.org/flymine/service")
            query = service.new_query("Transcript")
            query.add_view("Transcript.sequence.residues")
            query.add_constraint("Transcript", "LOOKUP", transcript, "D. melanogaster", code="A")


            # Write each of the returned transcript sequences to a new FASTA file:
            for row in query.rows():

                with open("/var/www/html/FlyTED2/transcription.FASTA", "w") as trans:
                    trans.write(">test")
                    trans.write("\n")
                    trans.write(row["Transcript.sequence.residues"])
                os.system("makeblastdb -in /var/www/html/FlyTED2/transcription.FASTA -dbtype nucl>/dev/null")

            bounds = []
            end = []
            res = []
            restwo = []

            # Three prime probe blast:
            three_query = "blastn -query three.seq -db /var/www/html/FlyTED2/transcription.FASTA -task blastn-short -strand plus -out three.txt -outfmt 10"
            os.system(three_query)



            # Open the result of the 3' Blast on the current transcript and check that their was indeed a result by checking if the file was empty.
            with open('three.txt') as blast_results_file:

                if not (os.stat("three.txt").st_size == 0):
                    results = blast_results_file.readline()
                    formatted_match = results.split(",")

                    bounds.append(formatted_match[8])  # Log the base position of start of BLAST match.
                    end.append(formatted_match[9])  # Log the base position of the end of the BLAST match.

                    low = int(formatted_match[8])  # IS THIS NEEDED?
                    high = int(formatted_match[9])
                    res.append(low)
                    res.append(high)

            # Five prime probe blast:
            five_query = "blastn -query five.seq -db /var/www/html/FlyTED2/transcription.FASTA -task blastn-short -strand plus -out five.txt -outfmt 10"
            os.system(five_query)


            # Open the result of the 5' Blast on the current transcript and check that their was indeed a result by checking if the file was empty.
            with open('five.txt') as blast_results_file:
                if not (os.stat("five.txt").st_size == 0):
                    results = blast_results_file.readline()
                    formatted_match = results.split(",")


                    bounds.append(formatted_match[9])  # Log the end of the sequence Helen and the team studied.
                    end.append(formatted_match[8])

                    low = int(formatted_match[8])  # DO WE NEED THIS?
                    high = int(formatted_match[9])
                    restwo.append(low)
                    restwo.append(high)

            # Read original transcript sequence from FASTA file into string:
            gene_sequence = ""
            line1 = True

            with open('/var/www/html/FlyTED2/transcription.FASTA', 'r') as myfile:

                for line in myfile:
                    if not line1:
                        gene_sequence += line
                    if line1:
                        line1 = False



            if (len(restwo) > 1 and len(res) > 1 and (res[1] - restwo[0] + 1) > 0):
                # Calculate subsequence that was the focus sequence studied by Helen and the team:

                target_sequence = colored(textwrap.fill(gene_sequence[restwo[0] - 1:res[1]], 50), 'magenta', 'on_grey',
                                          attrs=['bold'])


                temp = [data_set[0], gene_sequence[res[0] - 1:res[1]].replace('\n', '').replace('\r', ''),
                        gene_sequence[restwo[0] - 1:restwo[1]].replace('\n', '').replace('\r', ''),
                        gene_sequence[restwo[0] - 1:res[1]], transcript, len(gene_sequence[restwo[0] - 1:res[1]]), gene_sequence, len(gene_sequence)]
                final.append(temp)






# WRITE RESULTS TO FINAL FILE:
r = []
with open('probes.csv', 'w') as csvfile:
    fieldnames = ['PROBE', 'TRANSCRIPT', '3-PRIME-SEQUENCE', '5-PRIME-SEQUENCE', 'TARGET-SEQUENCE', 'TARGET-SEQUENCE-LENGTH', 'TRANSCRIPT-SEQUENCE', 'TRANSCRIPT-SEQUENCE-LENGTH' ]
    tester = csv.DictWriter(csvfile, fieldnames=fieldnames)
    for row in final:
        r.append({'PROBE': row[0], 'TRANSCRIPT': row[4], '3-PRIME-SEQUENCE': row[1], '5-PRIME-SEQUENCE': row[2],
                         'TARGET-SEQUENCE-LENGTH': row[5], 'TARGET-SEQUENCE': row[3], 'TRANSCRIPT-SEQUENCE': row[6], 'TRANSCRIPT-SEQUENCE-LENGTH':row[7]})
        tester.writerow({'PROBE': row[0], 'TRANSCRIPT': row[4], '3-PRIME-SEQUENCE': row[1], '5-PRIME-SEQUENCE': row[2],
                         'TARGET-SEQUENCE': row[3], 'TARGET-SEQUENCE-LENGTH': row[5], 'TRANSCRIPT-SEQUENCE': row[6], 'TRANSCRIPT-SEQUENCE-LENGTH':row[7]})



#######################################################################################
class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)
#print(r)
python_array = [1,2,3,4,5,6]
json_string = json.dumps(r)
# Now you have a string containing the python array, you can print it if you like
print(json_string)

#now we convert the json string back into a python array
new_var = json.loads(json_string)
#print(new_var[0]['PROBE'])
#print(type(new_var))

# filename = '/home/benedict/Desktop/FlyTED2/'+ sys.argv[1]
# # opening the file with w+ mode truncates the file
# f = open(filename, "w+")
# f.close()
