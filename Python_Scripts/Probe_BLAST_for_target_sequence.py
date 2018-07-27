import sys
from termcolor import colored, cprint

# import os
#
# # [1] TEST WITH DUMMY TRANSCRIPT SEQUENCE --> STATUS: PASSED
# # [2] FETCH TRANSCRIPT FROM FLYMINE USING API --> STATUS: PASSED
# # [3] PIPE DATA INTO Probe_BLAST_for_target_sequence.py --> STATUS: PASSED
# # [4] CREATE Probe_BLAST_for_target_sequence.py FILE FOR TRANSCRIPT SEQUENCE --> STATUS:
#
#
# # Script to find sub-sequence that was target of expression pattern recording.
# # [1] Need to find the actual sub-sequence.
# # [2] Need to find the length of the sub-sequence.
#
# # The constant sequence used in the 3' probe is:
# #-----------------------------------------------
# # [GTAATTAACCCTCACTAAAGG]
# #-----------------------------------------------
#
#
# # Used Bio complement method to find complement of the 3' before beforming blastn on it
# #--------------------------------------------------------------------------------------
# # []? from Bio.Seq import Seq
# #
# # seq = Seq("TCGGGCCC")
# #
# # print (seq.reverse_complement())
#
#
#
# bounds = []
# end = []
# res = []
# restwo = []
#
# # Three prime probe Probe_BLAST_for_target_sequence.py:
# three_query = "blastn -query three.seq -db CG5762.FASTA -task blastn-short -out three.txt -outfmt 10"
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
#     res.append(low -1 + (int(low/70)))
#     res.append(high + (int(high / 70)))
#
#
# # Five prime probe Probe_BLAST_for_target_sequence.py:
# five_query = "blastn -query five.seq -db CG5762.FASTA -task blastn-short -out five.txt -outfmt 10"
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
#     restwo.append(low -1 + (int(low/70)))
#     restwo.append(high + (int(high / 70)))
#
# # Read in the original sequence:
# gene_sequence = ""
# line1 = True
# with open('CG5762.FASTA', 'r') as myfile:
#
#     for line in myfile:
#         if not line1:
#             gene_sequence += line
#         if line1:
#             line1 = False
# print("")
# print("---GENE SEQUENCE---")
# print("")
# print(gene_sequence)
# print("LENGTH: ")
# print(len(gene_sequence))
# lower = int(bounds[0])
#
#
# # Calulate subsequence that was the focus sequence studied by Helen and the team:
#
# print("------------------ FINAL ---------------------")
# print("")
#
# upper = int(bounds[0])
# lower = int(bounds[1])
# lower_adjust = int(lower/70)
# upper_adjust = int(upper/70)
#
# print("---5'---")
# print(gene_sequence[restwo[0]:restwo[1]].replace('\n', '').replace('\r', ''))
# print("")
# print("---3'---")
# print(gene_sequence[res[0]:res[1]].replace('\n', '').replace('\r', ''))
# print("")
# print("---TARGET--")
# print(gene_sequence[restwo[0]:res[1]])  # Does it include the probes or not?
# print("")
# print("---TARGET LOCATION ON GENE---")
# print("START: " + str(restwo[1]))
# print("END: " + str(res[0]))
# print("")
# print("---TARGET LENGTH---")
# print("LENGTH: " + str(res[0]-restwo[1]))
#
# print("QUESTION: DO WE INLCUDE THE PROBE SEQUENCES IN THE TARGET SEQUENCE?")
#
#
# # TEST FLYMINE API:
# print("ATTEMPT TO RETIEVE TRANSCRIPT IDENTIFIERS FOR A GIVEN GENE: ")
# print("")
#
# from intermine.webservice import Service
# service = Service("http://www.flymine.org/flymine/service")
# query = service.new_query("Gene")
# query.add_view("transcripts.primaryIdentifier", "transcripts.length")
# query.add_sort_order("Gene.secondaryIdentifier", "ASC")
# query.add_constraint("Gene", "LOOKUP", "CG8564", code = "A")
#
# for row in query.rows():
#     print ("TRANSCRIPT ID: ", row["transcripts.primaryIdentifier"], "||", "TRANSCRIPT LENGTH: ", row["transcripts.length"])
#



# PIPE FETCHED TRANSCRIPT DATA AND 3', 5' PROBES INTO Probe_BLAST_for_target_sequence.py:

#NEED:
# [1] --> 5' PROBE FROM EXCEL PLACED IN .SEQ FILE
# [2] --> 3' PROBE FROM EXCEL PLACED IN .SEQ FILE
# [3] --> TRANSCRIPT(LONGEST) FROM FLYMINE FORMATTED TO Probe_BLAST_for_target_sequence.py DB

# [1] READ THE EXCEL DOCUMENT:
#
final = []
import csv
row_count = 0
with open('probes_disc_1.csv') as csvfile:
    rows = csv.reader(csvfile, delimiter=',', quotechar='|')
    row_count = sum(1 for row in rows)
    #print(row_count)

#[2] WRITE PROBE COUPLES TO TEXT FILE:
with open('preprocessed.csv', 'w') as csvfile:
    fieldnames = ['PROBE', '3-END', '5-END']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    with open('probes_disc_1.csv') as csvfile:
        rows = csv.reader(csvfile, delimiter=',', quotechar='|')
        counter = 0
        for row in rows:
                # PRINT THE PROBE:
                # print("PROBE: ", row[0].replace('-3T3','').replace('-5',''))
            probe_id = row[0].replace('-3T3', '').replace('-5', '')
                # PRINT THE 3':
                # print("3' SEQUENCE:",row[1])
            three_end = row[1]

                # PRINT THE 5'
            if (counter < (row_count / 2)):
                # print("5' SEQUENCE:", next(rows)[1])
                five_end = next(rows)[1]
                from Bio.Seq import Seq

            #print("")
            counter = counter + 1


            # REMOVE CONSTANT FROM THE 3' END:
            three_constant_type_A = 'GTAATTAACCCTCACTAAAGGG'
            three_constant_type_B = 'CTAATTAACCCTCACTAAAGGG'
            three_end = three_end.replace(three_constant_type_A, '').replace(three_constant_type_B, '')
            writer.writerow({'PROBE': probe_id, '3-END': three_end, '5-END': five_end})
        print("FINISHED THE WRITE.")

import os
from Bio.Seq import Seq
print("")
print("")
# READ PROBE SETS LINE BY LINE AND PERFORM Probe_BLAST_for_target_sequence.py ON THEM:
with open('preprocessed.csv') as csvfile:
    # Write answers to text file:


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
            #seq_file_five.write("TTCAGTTCGTTACACGCCGCCTCTTTGGAAGCAGGCTA")
            seq_file_five.write(data_set[2])

        # Need to fetch transcript sequences:
        from intermine.webservice import Service
        service = Service("http://www.flymine.org/flymine/service")
        query = service.new_query("Gene")
        query.add_view("transcripts.primaryIdentifier", "transcripts.length")
        query.add_sort_order("transcripts.length", "DESC")
        query.add_constraint("Gene", "LOOKUP", data_set[0], code = "A")
        transcript_primary_identifiers = []
        for row in query.rows():
            print("")
            print ("TRANSCRIPT ID: ", row["transcripts.primaryIdentifier"], "||", "TRANSCRIPT LENGTH: ", row["transcripts.length"])
            transcript_primary_identifiers.append(row["transcripts.primaryIdentifier"])
        print("")
        print("")
        print("LIST OF TRANSCRIPT PRIMARY IDENTIFIERS: ", transcript_primary_identifiers)
        #transcript_primary_identifiers = transcript_primary_identifiers[0]  # FETCH LONGEST
        # Need to fetch sequence for each of the transcripts:
        for transcript in transcript_primary_identifiers:
            print("")
            print("NEXT TRANSCRIPT DATA BELOW:", transcript)
            print("")
            from intermine.webservice import Service
            service = Service("http://www.flymine.org/flymine/service")
            query = service.new_query("Transcript")
            query.add_view("Transcript.sequence.residues")
            query.add_constraint("Transcript", "LOOKUP",transcript, "D. melanogaster", code="A")
            for row in query.rows():
                #print("****", row["Transcript.sequence.residues"])
                with open("transcription.FASTA", "w") as trans:
                    trans.write(">test")
                    trans.write("\n")
                    trans.write(row["Transcript.sequence.residues"])
                os.system("makeblastdb -in transcription.FASTA -dbtype nucl")

            import textwrap

            bounds = []
            end = []
            res = []
            restwo = []

            # Three prime probe Probe_BLAST_for_target_sequence.py:
            three_query = "blastn -query three.seq -db transcription.FASTA -task blastn-short -strand plus -out three.txt -outfmt 10"
            os.system(three_query)

            print("THREE PRIME BEST MATCH:")
            print("")

            with open('three.txt') as blast_results_file:
                if not(os.stat("three.txt").st_size == 0):
                    results = blast_results_file.readline()
                    formatted_match = results.split(",")
                    print("<< RESULT >>")
                    print("Percentage Match: " + formatted_match[2])
                    print("Match Start: " + formatted_match[8])
                    print("Match End: " + formatted_match[9])
                    print("<< RESULT >>")
                    print("")

                    bounds.append(formatted_match[8])  # Log the start of the sequence Helen and the team studied.
                    end.append(formatted_match[9])

                    low = int(formatted_match[8])

                    high = int(formatted_match[9])

                    res.append(low)
                    res.append(high)

            # Five prime probe Probe_BLAST_for_target_sequence.py:
            five_query = "blastn -query five.seq -db transcription.FASTA -task blastn-short -strand plus -out five.txt -outfmt 10"
            os.system(five_query)

            print("FIVE PRIME BEST MATCH:")
            print("")

            with open('five.txt') as blast_results_file:
                if not(os.stat("five.txt").st_size == 0):
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

                    low = int(formatted_match[8])

                    high = int(formatted_match[9])

                    restwo.append(low)
                    restwo.append(high)

            # Read in the original sequence:
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

            # Calulate subsequence that was the focus sequence studied by Helen and the team:
            print("")
            print("------------------ FINAL DATA ---------------------")
            print("")


            if(len(res) > 1):
                print("---3'---")
                print(gene_sequence[res[0] - 1:res[1]].replace('\n', '').replace('\r', ''))
                print("")
            if (len(restwo) > 1):
                print("---5'---")
                print(gene_sequence[restwo[0] - 1:restwo[1]].replace('\n', '').replace('\r', ''))
                print("")

            if(len(restwo) > 1 and len(res) > 1 and (res[1] - restwo[0] + 1) > 0):
                print("---TARGET SEQUENCE--")
                #print(textwrap.fill(gene_sequence[restwo[0] - 1:res[1]], 50))  # Does it include the probes or not?
                target_sequence = colored(textwrap.fill(gene_sequence[restwo[0] - 1:res[1]], 50), 'magenta', 'on_grey', attrs=['bold'])
                print(target_sequence)
                print("")
                print("---TARGET LOCATION ON GENE---")
                print("START: " + str(restwo[0] - 1))
                print("END: " + str(res[1]))
                print("")
                print("---TARGET LENGTH---")
                print("LENGTH: " + str((res[1] - restwo[0]) + 1))
                print("QUESTION: DO WE INLCUDE THE PROBE SEQUENCES IN THE TARGET SEQUENCE?")
                temp = [data_set[0], gene_sequence[res[0] - 1:res[1]].replace('\n', '').replace('\r', ''), gene_sequence[restwo[0] - 1:restwo[1]].replace('\n', '').replace('\r', ''),  gene_sequence[restwo[0] - 1:res[1]], transcript]
                final.append(temp)


# WRITE RESULTS TO FINAL FILE:
with open('probes.csv', 'w') as csvfile:
    fieldnames = ['PROBE', 'TRANSCRIPT', '3-PRIME-SEQUENCE', '5-PRIME-SEQUENCE', 'TARGET-SEQUENCE']
    tester = csv.DictWriter(csvfile, fieldnames=fieldnames)
    for row in final:
        tester.writerow({'PROBE': row[0], 'TRANSCRIPT': row[4], '3-PRIME-SEQUENCE':row[1], '5-PRIME-SEQUENCE':row[2], 'TARGET-SEQUENCE':row[3]})

#with open('news.FASTA', 'w') as file:
 #   file.write(">test")
  #  file.write("\n")
   # file.write("ACCGAGTTATGCATGGTTAAACCCG")

#os.system("makeblastdb -in news.FASTA -dbtype nucl")



















































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
