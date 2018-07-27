import os
from Bio.Seq import Seq
import csv
from termcolor import colored
from intermine.webservice import Service
import textwrap

final = []  # Will contain the final data sets --> [PROBE, TRANSCRIPT ID, 3', 5', TARGET SEQUENCE]
row_count = 0

# [1] TEST WITH DUMMY TRANSCRIPT SEQUENCE --> STATUS: PASSED
# [2] FETCH TRANSCRIPT FROM FLYMINE USING API --> STATUS: PASSED
# [3] PIPE DATA INTO Probe_BLAST_for_target_sequence.py --> STATUS: PASSED
# [4] CREATE Probe_BLAST_for_target_sequence.py FILE FOR TRANSCRIPT SEQUENCE --> STATUS: PASSED
# [5] PERFORM BLAST and output target sequence to text file --> STATUS: PASSED

# Find number of entries in the csv, two lines per gene, one for 3' and one for 5'
with open('probes_disc_1.csv') as csvfile:
    rows = csv.reader(csvfile, delimiter=',', quotechar='|')
    row_count = sum(1 for row in rows)

# Write eacg dataset for each gene on one line
with open('preprocessed.csv', 'w') as csvfile:
    fieldnames = ['PROBE', '3-END', '5-END']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    with open('probes_disc_1.csv') as csvfile:
        rows = csv.reader(csvfile, delimiter=',', quotechar='|')
        counter = 0
        for row in rows:
            probe_id = row[0].replace('-3T3', '').replace('-5', '')
            three_end = row[1]

            if (counter < (row_count / 2)):
                five_end = next(rows)[1]

            counter = counter + 1

            # REMOVE CONSTANT SEQUENCE FROM THE 3' END:
            three_constant_type_A = 'GTAATTAACCCTCACTAAAGGG'
            three_constant_type_B = 'CTAATTAACCCTCACTAAAGGG'
            three_end = three_end.replace(three_constant_type_A, '').replace(three_constant_type_B, '')
            writer.writerow({'PROBE': probe_id, '3-END': three_end, '5-END': five_end})

        print("FINISHED THE WRITE.")

print("")
print("")

# READ PROBE SETS LINE BY LINE AND PERFORM BLASTN ON EACH PROBE AND THE TRANSCRIPT(S) FASTA FILES FOR A GIVEN GENE:
with open('preprocessed.csv') as csvfile:
    gene_probe_sets = csv.reader(csvfile)
    next(gene_probe_sets)
    for data_set in gene_probe_sets:
        print(data_set)

        # Write 3' probe to three.seq
        with open("three.seq", "w") as seq_file_three:
            three_constant_type_A = 'GTAATTAACCCTCACTAAAGGG'
            three_constant_type_B = 'CTAATTAACCCTCACTAAAGGG'
            three_end = data_set[1].replace(three_constant_type_A, '').replace(three_constant_type_B, '')
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

        # For each transcript found for a given gene:
        for row in query.rows():
            print("")
            print("TRANSCRIPT ID: ", row["transcripts.primaryIdentifier"], "||", "TRANSCRIPT LENGTH: ",
                  row["transcripts.length"])
            transcript_primary_identifiers.append(row["transcripts.primaryIdentifier"])

        print("")
        print("")
        print("LIST OF TRANSCRIPT PRIMARY IDENTIFIERS: ", transcript_primary_identifiers)

        # Need to fetch the associated sequence for each of the gene transcript ID's that were returned from the above query:
        for transcript in transcript_primary_identifiers:

            print("")
            print("<----- NEW TRANSCRIPT DATA BELOW ----->", transcript)
            print("")

            service = Service("http://www.flymine.org/flymine/service")
            query = service.new_query("Transcript")
            query.add_view("Transcript.sequence.residues")
            query.add_constraint("Transcript", "LOOKUP", transcript, "D. melanogaster", code="A")

            # Write each of the returned transcript sequences to a new FASTA file:
            for row in query.rows():
                with open("transcription.FASTA", "w") as trans:
                    trans.write(">test")
                    trans.write("\n")
                    trans.write(row["Transcript.sequence.residues"])
                os.system("makeblastdb -in transcription.FASTA -dbtype nucl")

            bounds = []
            end = []
            res = []
            restwo = []

            # Three prime probe blast:
            three_query = "blastn -query three.seq -db transcription.FASTA -task blastn-short -strand plus -out three.txt -outfmt 10"
            os.system(three_query)

            print("<----- THREE PRIME BEST MATCH ----->")
            print("")

            # Open the result of the 3' Blast on the current transcript and check that their was indeed a result by checking if the file was empty.
            with open('three.txt') as blast_results_file:
                if not (os.stat("three.txt").st_size == 0):
                    results = blast_results_file.readline()
                    formatted_match = results.split(",")
                    print("<< RESULT >>")
                    print("Percentage Match: " + formatted_match[2])
                    print("Match Start: " + formatted_match[8])
                    print("Match End: " + formatted_match[9])
                    print("<< RESULT >>")
                    print("")

                    bounds.append(formatted_match[8])  # Log the base position of start of BLAST match.
                    end.append(formatted_match[9])  # Log the base position of the end of the BLAST match.

                    low = int(formatted_match[8])  # IS THIS NEEDED?
                    high = int(formatted_match[9])
                    res.append(low)
                    res.append(high)

            # Five prime probe blast:
            five_query = "blastn -query five.seq -db transcription.FASTA -task blastn-short -strand plus -out five.txt -outfmt 10"
            os.system(five_query)

            print("<----- FIVE PRIME BEST MATCH ----->")
            print("")

            # Open the result of the 5' Blast on the current transcript and check that their was indeed a result by checking if the file was empty.
            with open('five.txt') as blast_results_file:
                if not (os.stat("five.txt").st_size == 0):
                    results = blast_results_file.readline()
                    formatted_match = results.split(",")
                    print("<< RESULT >>")
                    print("Percentage Match: " + formatted_match[2])
                    print("Match Start: " + formatted_match[8])
                    print("Match End: " + formatted_match[9])
                    print("<< RESULT >>")
                    print("")

                    bounds.append(formatted_match[9])  # Log the end of the sequence Helen and the team studied.
                    end.append(formatted_match[8])

                    low = int(formatted_match[8])  # DO WE NEED THIS?
                    high = int(formatted_match[9])
                    restwo.append(low)
                    restwo.append(high)

            # Read original transcript sequence from FASTA file into string:
            gene_sequence = ""
            line1 = True

            with open('transcription.FASTA', 'r') as myfile:

                for line in myfile:
                    if not line1:
                        gene_sequence += line
                    if line1:
                        line1 = False
            print("")
            print("---GENE SEQUENCE---")
            print("")
            print(textwrap.fill(gene_sequence, 50))
            print("")
            print("LENGTH: ")
            print(len(gene_sequence))

            print("")
            print("------------------ FINAL DATA ---------------------")
            print("")

            if (len(res) > 1):
                print("---3'---")
                print(gene_sequence[res[0] - 1:res[1]].replace('\n', '').replace('\r', ''))
                print("")
            if (len(restwo) > 1):
                print("---5'---")
                print(gene_sequence[restwo[0] - 1:restwo[1]].replace('\n', '').replace('\r', ''))
                print("")

            if (len(restwo) > 1 and len(res) > 1 and (res[1] - restwo[0] + 1) > 0):
                # Calculate subsequence that was the focus sequence studied by Helen and the team:
                print("---TARGET SEQUENCE--")
                target_sequence = colored(textwrap.fill(gene_sequence[restwo[0] - 1:res[1]], 50), 'magenta', 'on_grey',
                                          attrs=['bold'])
                print(target_sequence)
                print("")

                print("---TARGET LOCATION ON GENE---")
                print("START: " + str(restwo[0] - 1))
                print("END: " + str(res[1]))
                print("")

                print("---TARGET LENGTH---")
                print("LENGTH: " + str((res[1] - restwo[0]) + 1))
                print("QUESTION: DO WE INLCUDE THE PROBE SEQUENCES IN THE TARGET SEQUENCE?")

                temp = [data_set[0], gene_sequence[res[0] - 1:res[1]].replace('\n', '').replace('\r', ''),
                        gene_sequence[restwo[0] - 1:restwo[1]].replace('\n', '').replace('\r', ''),
                        gene_sequence[restwo[0] - 1:res[1]], transcript]
                final.append(temp)

# WRITE RESULTS TO FINAL FILE:
with open('probes.csv', 'w') as csvfile:
    fieldnames = ['PROBE', 'TRANSCRIPT', '3-PRIME-SEQUENCE', '5-PRIME-SEQUENCE', 'TARGET-SEQUENCE']
    tester = csv.DictWriter(csvfile, fieldnames=fieldnames)
    for row in final:
        tester.writerow({'PROBE': row[0], 'TRANSCRIPT': row[4], '3-PRIME-SEQUENCE': row[1], '5-PRIME-SEQUENCE': row[2],
                         'TARGET-SEQUENCE': row[3]})

#
#
#
# import os
# #
# #    [1] TEST WITH DUMMY TRANSCRIPT SEQUENCE --> STATUS: PASSED
# # # [2] FETCH TRANSCRIPT FROM FLYMINE USING API --> STATUS: PASSED
# # # [3] PIPE DATA INTO Probe_BLAST_for_target_sequence.py --> STATUS: PASSED
# # # [4] CREATE Probe_BLAST_for_target_sequence.py FILE FOR TRANSCRIPT SEQUENCE --> STATUS:
# #
# #
# # # Script to find sub-sequence that was target of expression pattern recording.
# # # [1] Need to find the actual sub-sequence.
# # # [2] Need to find the length of the sub-sequence.
# #
# # # The constant sequence used in the 3' probe is:
# # #-----------------------------------------------
# # # [GTAATTAACCCTCACTAAAGG]
# # #-----------------------------------------------
# #
# #
# # # Used Bio complement method to find complement of the 3' before beforming blastn on it
# # #--------------------------------------------------------------------------------------
# # # []? from Bio.Seq import Seq
# # #
# # # seq = Seq("TCGGGCCC")
# # #
# # # print (seq.reverse_complement())
# #
#
# print("")
# print("---------------------------------------------------")
# print("")
# import textwrap
# bounds = []
# end = []
# res = []
# restwo = []
#
# # Three prime probe Probe_BLAST_for_target_sequence.py:
# three_query = "blastn -query three.seq -db transcription.FASTA -task blastn-short -out three.txt -outfmt 10"
# os.system(three_query)
#
# print("THREE PRIME BEST MATCH:")
# print("")
#
# with open('three.txt') as blast_results_file:
#     results = blast_results_file.readline()
#     formatted_match = results.split(",")
#     print("<< RESULT >>")
#     print("Percentage Match: " + formatted_match[2])
#     print("Match Start: " + formatted_match[8])
#     print("Match End: " + formatted_match[9])
#     print("<< RESULT >>")
#     print("")
#
#     bounds.append(formatted_match[8])   # Log the start of the sequence Helen and the team studied.
#     end.append(formatted_match[9])
#
#     low = int(formatted_match[8])
#
#     high = int(formatted_match[9])
#
#     res.append(low)
#     res.append(high)
#
#
# # # Five prime probe Probe_BLAST_for_target_sequence.py:
# five_query = "blastn -query five.seq -db transcription.FASTA -task blastn-short -out five.txt -outfmt 10"
# os.system(five_query)
#
# print("FIVE PRIME BEST MATCH:")
# print("")
#
# with open('five.txt') as blast_results_file:
#     results = blast_results_file.readline()
#     formatted_match = results.split(",")
#     print("<< RESULT >>")
#     print("Percentage Match: " + formatted_match[2])
#     print("Match Start: " + formatted_match[8])
#     print("Match End: " + formatted_match[9])
#     print("<< RESULT >>")
#     print("")
#
#     bounds.append(formatted_match[9])  # Log the end of the sequence Helen and the team studied.
#     end.append(formatted_match[8])
#
#
#     low = int(formatted_match[8])
#
#     high = int(formatted_match[9])
#
#     restwo.append(low)
#     restwo.append(high)
#
# # # Read in the original sequence:
# gene_sequence = ""
# line1 = True
#
# with open('transcription.FASTA', 'r') as myfile:
#
#     for line in myfile:
#         if not line1:
#             gene_sequence += line
#         if line1:
#             line1 = False
# print("")
# print("---GENE SEQUENCE---")
# print("")
# print(textwrap.fill(gene_sequence, 50))
# print("")
# print("LENGTH: ")
# print(len(gene_sequence))
# lower = int(bounds[0])
#
#
# # Calulate subsequence that was the focus sequence studied by Helen and the team:
# print("")
# print("------------------ FINAL DATA ---------------------")
# print("")
#
# print("---3'---")
# print(gene_sequence[res[0] - 1 :res[1]].replace('\n', '').replace('\r', ''))
# print("")
# print("---5'---")
# print(gene_sequence[restwo[0] - 1:restwo[1]].replace('\n', '').replace('\r', ''))
# print("")
# print("---TARGET SEQUENCE--")
# target_sequence = colored(textwrap.fill(gene_sequence[restwo[0] - 1:res[1]], 50), 'magenta', 'on_grey', attrs=['bold'])
# print(target_sequence)
# print("")
# print("---TARGET LOCATION ON GENE---")
# print("START: " + str(restwo[0] - 1))
# print("END: " + str(res[1]))
# print("")
# print("---TARGET LENGTH---")
# print("LENGTH: " + str((res[1]-restwo[0]) + 1))
#
# print("QUESTION: DO WE INLCUDE THE PROBE SEQUENCES IN THE TARGET SEQUENCE?")
