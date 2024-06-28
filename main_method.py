'''
Author: Nimna Gamage
Date: 24/04/2023
Final_Software_Project
Project_No. 2: Sequence aligner
Separate main program to demonstrate the use of the 'Seq_aligner' class, class variables and methods
Input : 'Seq_records.fasta' file
Output : calling all methods of the 'seq_aligner' class
'''

# Import the python file containing the 'seq_aligner' class
import seq_aligner_class

# main method
if __name__ == "__main__":

    # Get the input FASTA file containing multiple DNA, RNA, or amino acid sequences as a user input
    print("Seq_aligner Class!")
    fasta_File = input("Enter the Fasta File Name(e.g.:'Seq_records.fasta') : ")
    # Create sequence objects from the input fasta file
    records = seq_aligner_class.Seq_aligner(fasta_File)
    print("For further user inputs, input the relevant fasta sequence header without the '>' mark(e.g.:'seq1p')!")
    # print(records.seq_records)


    # Method 1 - "get_Seq_Type"
    print("\n#Method 1\nGet the type of the fasta sequence!")
    # Get the user input
    type_Seq_Id = input("Enter the Fasta Sequence Id[e.g.:'seq1p'] : ")
    # loop through each record
    for seq_rec in range(len(records.seq_records)):
        # select the sequence with the same id entered as a user input
        if records.seq_records[seq_rec].id == type_Seq_Id:
            print(f"The {records.seq_records[seq_rec].id} is a {records.get_Seq_Type(records.seq_records[seq_rec].seq)}.")


    # Method 5 - "find_GenBank_Seq" - Static Method
    print("\n#Method 5\nFind GenBank Sequence for the Accession Number!")
    # Get the input from the user
    accession_No = input("Enter the Accession Number (with version) [e.g.:'NM_000188.3']: ")
    seq_aligner_class.Seq_aligner.find_GenBank_Seq(accession_No)
    print(f"The Genbank sequence of the {accession_No} is saved successfully in fasta format!")


    # Method 6 - "cal_melting_temp"
    print("\n#Method 6\nCalculate the melting temperature of any nucleotide sequence!")
    # Get the user input
    mt_Seq_Id = input("Enter the Fasta Sequence Id[e.g.:'seq2cds'] : ")
    # loop through each record
    for seq_rec in range(len(records.seq_records)):
        # select the sequence with the same id entered as a user input
        if records.seq_records[seq_rec].id == mt_Seq_Id:
            records.cal_melting_temp(records.seq_records[seq_rec].seq)


    # Method 7 - "cal_isoelectric_point"
    print("\n#Method 7\nCalculate the isoelectric point of any amino acid sequence!")
    # Get the user input
    ip_Seq_Id = input("Enter the Fasta Sequence Id[e.g.:'seq2p'] : ")
    # loop through each record
    for seq_rec in range(len(records.seq_records)):
        # select the sequence with the same id entered as a user input
        if records.seq_records[seq_rec].id == ip_Seq_Id:
            records.cal_isoelectric_point(records.seq_records[seq_rec].seq)


    # Method 8 - "cal_aromaticity"
    print("\n#Method 8\nCalculate the aromaticity of any amino acid sequence!")
    # Get the user input
    aro_Seq_Id = input("Enter the Fasta Sequence Id[e.g.:'seq2p'] : ")
    # loop through each record
    for seq_rec in range(len(records.seq_records)):
        # select the sequence with the same id entered as a user input
        if records.seq_records[seq_rec].id == aro_Seq_Id:
            records.cal_aromaticity(records.seq_records[seq_rec].seq)


    # Method 2 - "cal_similarity_score"
    print("\n#Method 2\nGet the global pairwise similarity scores between all sequence pairs in the input file into a text file!")
    records.cal_similarity_score()
    print("File saved successfully!")


    # Method 3 - "dot_Plot_Scatter"
    print("\n#Method 3\nGenerate a dot plot(as a scatter plot) to represent pairwise similarity!")
    # Get the user input
    seq1_Id = input("Enter The Fasta Sequence 1 (x-axis) [e.g.:'seq1p'] : ")
    seq2_Id = input("Enter The Fasta Sequence 2 (y-axis) [e.g.:'seq2p'] : ")
    # loop through each record
    for seq_rec in range(len(records.seq_records)):
        # select the sequence with the same id entered as a user input
        if records.seq_records[seq_rec].id == seq1_Id:
            seq_one = records.seq_records[seq_rec]
        elif records.seq_records[seq_rec].id == seq2_Id:
            seq_two = records.seq_records[seq_rec]
    records.dot_Plot_Scatter(seq_one, seq_two)


    # Method 4 - "run_Blast"
    print("\n#Method 4\nRun an appropriate BLAST algorithm type, depending on the type of the sequence!")
    print("Running Blast Search...")
    records.run_Blast()






