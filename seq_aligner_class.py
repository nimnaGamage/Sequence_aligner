'''
Author: Nimna Gamage
Date: 24/04/2023
Final_Software_Project
Project_No. 2: Sequence aligner
A class and it's methods to align biological sequences (DNA, RNA, and amino acid) using Biopython and other 3rd party packages.
Following methods take a FASTA file containing multiple DNA, RNA, or amino acid sequences.
This class consists with 8 methods ;
    1-To output the type of a given FASTA sequence
    2-To calculate the global pairwise similarity score between all sequence pairs in the input file and write the output in a text file
    3-To generate a dot plot(as a scatter plot) to represent pairwise similarity for two given sequences
    4-To run BLAST searches on several sequence inputs in a single FASTA file
      (run an appropriate BLAST algorithm type, depending on the type of the sequence)
    5-To find a sequence from GenBank when its accession number is given, and save it in FASTA format
    6-To calculate the melting temperature of a given nucleotide sequence
    7-To calculate the isoelectric point of a given amino acid sequence
    8-To calculate the aromaticity of a given amino acid sequence
'''

# Import python packages

import re
import pylab
from collections import OrderedDict
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import pairwise2, SeqIO, Entrez
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

'''Sequence_aligner Class'''

class Seq_aligner:

    # Constructor method with input parameters
    def __init__(self, fastaFile):
        # Instance variables
        self.seq_records = list(SeqIO.parse(fastaFile, "fasta"))


    # 1
    # method to output the type of a given FASTA sequence
    def get_Seq_Type(self, sequence):

        # input sequence is converted to a string
        seq = str(sequence)

        # re module search for DNA sequence
        if re.search(r"[^ACGT]", seq, re.IGNORECASE) is None:
            return "DNA Sequence"
        # re module search for RNA sequence
        if re.search(r"[^ACGU]", seq, re.IGNORECASE) is None:
            return "RNA Sequence"
        # re module search for Protein sequence
        elif re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", seq, re.IGNORECASE) is None:
            return "Amino Acid Sequence"


    # 2
    # method to calculate the global pairwise similarity score between all sequence pairs in the input file and write the output in a text file
    def cal_similarity_score(self):

        # Create a dictionary to store global pairwise similarity scores
        dictSimilarityScore = {}

        # Loop through all sequence pairs in the input fasta file
        for seq1 in range(len(self.seq_records)):
            for seq2 in range(len(self.seq_records)):
                seq_1 = self.seq_records[seq1].seq
                seq_2 = self.seq_records[seq2].seq
                if seq_1 != seq_2 and self.get_Seq_Type(seq_1) == self.get_Seq_Type(seq_2):
                    # perform sequence alignment and calculate the global pairwise similarity score
                    alignment_score = pairwise2.align.globalxx(seq_1, seq_2, score_only=True)
                    # Create the key of the dictionary
                    dictSimilarityScore_key = (self.seq_records[seq1].id, self.seq_records[seq2].id)
                    # 'dictSimilarityScore' dictionary
                    dictSimilarityScore[dictSimilarityScore_key] = alignment_score

        # get the scores in descending order
        sorted_Score_dict = OrderedDict(sorted(dictSimilarityScore.items(), key=lambda t: t[1], reverse=True))

        # write the results into a text file
        with open ("2_Similarity_Scores.txt", 'w') as file:
            # write the header
            file.write("Sequence_combination\t\tSmiliarity_Score\n")
            # writing the similarity score for all sequence pairs in the input file in descending order
            for key, value in sorted_Score_dict.items():
                file.write(key[0] + " and " + key[1] + "\t-->\t" + str(value) + "\n")

        # Return the text file
        return file


    # 3
    # method to generate a dot plot(as a scatter plot) to represent pairwise similarity for two given sequences
    def dot_Plot_Scatter(self, seq1, seq2):

        # Get the two sequences and transform both into an one format(upper())
        seq_1 = seq1.seq.upper()
        seq_2 = seq2.seq.upper()
        # Calculate pairwise similarity score only between same type sequences
        if self.get_Seq_Type(seq_1) == self.get_Seq_Type(seq_2):
            # define the size of the window
            window = 7
            # Create two empty dictionaries
            dict_one = {}
            dict_two = {}
            # Iterate over two sequences and divide them into overlapping windows of size '7'
            # Subsequences of length window are stored in the 'section_dict' dictionaries along with their starting positions in the original sequence
            for (seq, section_dict) in [
                (seq_1, dict_one),
                (seq_2, dict_two),
            ]:
                for i in range(len(seq) - window):
                    section = seq[i: i + window]
                    try:
                        section_dict[section].append(i)
                    except KeyError:
                        section_dict[section] = [i]
            # get the sub-sequences found in both sequences
            matches = set(dict_one).intersection(dict_two)

            # Create lists of x and y coordinates  for scatter plot
            x = []
            y = []
            # Iterate over matching sub-sequences and retrieve the starting positions in both sequences
            for section in matches:
                for i in dict_one[section]:
                    for j in dict_two[section]:
                        # store starting positions of sub-sequences in 'x' and 'y'
                        x.append(i)
                        y.append(j)

            # Create the scatter plot using 'pylab' module
            # Clear any prior graph
            pylab.cla()
            pylab.gray()
            pylab.scatter(x, y)
            pylab.xlim(0, len(seq1) - window)
            pylab.ylim(0, len(seq2) - window)
            # Label of x-axis
            pylab.xlabel("%s (length %i bp)" % (seq1.id, len(seq1)))
            # Label of y-axis
            pylab.ylabel("%s (length %i bp)" % (seq2.id, len(seq2)))
            # Title of the Scatter plot
            pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
            # Save the plot
            pylab.savefig("3_Dot_Plot_Scatter")
            pylab.show()
            print("Figure saved successfully!")
        # Exception handling, if the two sequences are not in the same type
        else:
            print("Entered two sequences are not compatible!")
        return


    # 4
    # method to run BLAST searches on several sequence inputs in a single FASTA file
    # (run an appropriate BLAST algorithm type, depending on the type of the sequence)
    def run_Blast(self):

        # Loop through each and every sequence record
        for seqq in range(len(self.seq_records)):
            # Obtain the sequence
            sequence = self.seq_records[seqq].seq
            # Perform blastn for DNA and RNA Sequences
            if self.get_Seq_Type(sequence) == "DNA Sequence" or self.get_Seq_Type(sequence) == "RNA Sequence":
                # Run blastn over the internet
                result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
                # Save the blast results to an XML file
                with open(f"4_Blast_Result_{self.seq_records[seqq].id}.xml", "w") as out_handle:
                    out_handle.write(result_handle.read())
                #print(f"Blast result for {self.seq_records[seqq].id} is saved successfully!")
            # Perform blastp for Amino Sequences
            elif self.get_Seq_Type(sequence) == "Amino Acid Sequence":
                # Run blastp over the internet
                result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
                # Save the blast results to an XML file
                with open(f"4_Blast_Result_{self.seq_records[seqq].id}.xml", "w") as out_handle:
                    out_handle.write(result_handle.read())
                #print(f"Blast result for {self.seq_records[seqq].id} is saved successfully!")
            # Exception handling, if the sequence type is not recognizable
            else:
                print("Sequence type of the input sequence is not recognized!")
                continue

            # # Save the blast results to an XML file
            # with open(f"4_Blast_Result_{self.seq_records[seqq].id}.xml", "w") as out_handle:
            #     out_handle.write(result_handle.read())
            # # Close the result handle
            # result_handle.close()
            # # Return the file
            # return out_handle
            print(f"Blast result for {self.seq_records[seqq].id} is saved successfully!")


    # 5
    # method to find a sequence from GenBank when its accession number is given, and save it in FASTA format
    @staticmethod
    def find_GenBank_Seq(accession):

        # provide E-mail
        Entrez.email = "nimnagamage65@gmail.com"
        # retrieve the fasta sequence from Entrez
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        # close the handle
        handle.close()
        # open the fasta file
        file = open("5_Gb_Seq_{}.fasta".format(accession), 'w')
        # write the header and the sequence of the fasta file
        file.write(">{}\n{}".format(record.description, record.seq))
        # close the file
        file.close()
        return


    # 6
    # method to calculate the melting temperature of a given nucleotide sequence
    def cal_melting_temp(self, seq):

        # Melting temperature is calculated only for nucleotide sequences
        if self.get_Seq_Type(seq) == "DNA Sequence" or self.get_Seq_Type(seq) == "RNA Sequence":
            # Calculate melting temperature using the nearest-neighbor method with default parameters
            meltingTemp = mt.Tm_NN(seq)
            print(f"The melting temperature of the input sequence is {meltingTemp:.3f}°C")
            # Return the melting temperature
            return meltingTemp
        # for amino acid or other types of sequences
        else:
            print("Can not define a melting temperature for this sequence!")
            return


    # 7
    # method to calculate the isoelectric point of a given amino acid sequence
    def cal_isoelectric_point(self, seq):

        # Check the type of the input sequence
        seqType = self.get_Seq_Type(seq)

        # select only the protein sequences
        if seqType == "Amino Acid Sequence":
            # Creating a ProteinAnalysis object from the amino acid sequence
            analyzed_seq = ProteinAnalysis(seq)

            # Calculate the isoelectric point using the 'isoelectric_point()' method
            iso_point = analyzed_seq.isoelectric_point()

            print(f"The isoelectric point of the input sequence is {iso_point:.3f}")
            # return the isoelectric point
            return iso_point
        # for non-amino_acid sequences
        else:
            print("Can not define isoelectric point for this sequence!")
            return


    # 8
    # method to calculate the aromaticity of a given amino acid sequence
    def cal_aromaticity(self, seq):

        # Check the type of the input sequence
        seqType = self.get_Seq_Type(seq)

        # select only the protein sequences
        if seqType == "Amino Acid Sequence":
            # Creating a ProteinAnalysis object from the amino acid sequence
            analyzed_seq = ProteinAnalysis(seq)

            # Calculate the aromaticity using the 'aromaticity()' method
            aromaticity_no = analyzed_seq.aromaticity()

            print(f"The aromaticity of the input sequence is {aromaticity_no:.3f}")
            # return the aromaticity
            return aromaticity_no
        else:
            print("Can not define aromaticity for this sequence!")
            return





