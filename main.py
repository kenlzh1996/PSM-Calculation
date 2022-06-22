"""
Calculates the percentage of amino acid residues that are labeled relative to the number that should be labeled
"""
import pandas as pd

# functions
def Read_table(File):
    """Load the user input file path"""
    DF = pd.read_table(File)
    peptide_col_index = DF.columns.get_loc('Peptide sequence')
    sequence_list = DF.iloc[:,peptide_col_index].to_list()
    if is_by_fraction == "y":
        num_of_fraction = DF['Fraction'].unique()
    else:
        num_of_fraction = None
    num_of_sequence = len(sequence_list)

    return DF, peptide_col_index, sequence_list, num_of_sequence, num_of_fraction

def Check_Sequences(Sequence_List):
    for string in Sequence_List:
        if string.startswith('+304.207'):
            if string.find("K") != -1: # found K, but is the K tagged?
                if string.count("K+304.207") == string.count("K"): # found K with tag --> correct
                    correctly_labeled.append(string)
                else: # no tag --> in correct
                    incorrectly_labeled.append(string)
            else: # start with tag but no K --> correct
                correctly_labeled.append(string)
        else: # not start with tag --> incorrect
            incorrectly_labeled.append(string)

    return correctly_labeled

def At_Least_One_Label(Sequence_List):
    at_least_one_label = 0
    for string in Sequence_List:
        if string.find("+304.207") != -1:
            at_least_one_label += 1

    return at_least_one_label

def Amino_Residues_Label(Sequence_List, Correctly_Labeled_List):
    """
    total # of labelable amino acids residues =
    + # of N-terminal
    - # of sequence start with K
    - # of sequence start with +304.207K
    + # of K found in sequence
    """
    total_N_terminal = len(Sequence_List)
    total_seq_startwith_K = 0
    total_seq_start_with_tagged_K = 0
    total_K_residues = 0
    for string in Sequence_List:
        if string.startswith("K"):
            total_seq_startwith_K += 1
        if string.startswith("+304.207K"):
            total_seq_start_with_tagged_K += 1
        total_K_residues += string.count("K")

    total_possible_label_residues = total_N_terminal - total_seq_startwith_K - total_seq_start_with_tagged_K + total_K_residues

    correctly_labeled_N_terminal = len(Correctly_Labeled_List)
    correctly_labeled_K_residues = 0
    for string in Correctly_Labeled_List:
        correctly_labeled_K_residues += string.count("K")
        if string.startswith("+304.207K"):
            correctly_labeled_K_residues -= 1

    total_correctly_labeled_residues = correctly_labeled_N_terminal + correctly_labeled_K_residues

    return total_possible_label_residues, total_correctly_labeled_residues, correctly_labeled_N_terminal, correctly_labeled_K_residues


def Print_Report(Correct, Incorrect, At_Least_One, Total, Total_Possible_Residues, Total_Correct_Residues, Correct_N, Correct_Lysine):
    print("\n**************************** Report **************************")
    print(f"Total number of peptide sequence: {Total}")
    print(f"Correctly labeled peptide sequence: {len(Correct)}")
    print(f"Incorrectly labeled peptide sequence: {len(Incorrect)}")
    print(f"% of peptide sequence that has at least one TMT label: {round(At_Least_One / Total * 100, 2)}%")
    print(f"% of correctly labeled peptide sequence: {round(len(Correct) / Total * 100, 2)}%")
    print(f"% of incorrectly labeled peptide sequence: {round(100 - len(Correct) / Total * 100, 2)}%")
    print("---------------------------------------------------------------")
    print(f"Total number of possible labeled amino acid residues: {Total_Possible_Residues}")
    print(f"Total number of correctly labeled amino acid residues: {Total_Correct_Residues}")
    print(f"Total number of correctly labeled N-terminus: {Correct_N}")
    print(f"Total number of correctly labeled Lysine residues: {Correct_Lysine}")
    print(f"% of correctly labeled amino acid residues: {round(Total_Correct_Residues / Total_Possible_Residues *100, 2)}%")
    print("**************************** End *****************************\n")

def Print_Report_By_Fraction(Fraction_Num, Correct, Incorrect, At_Least_One, Total):
    print("\n**************************** Report **************************")
    print(f"Fraction Number: {Fraction_Num}")
    print(f"Total number of peptide sequence: {Total}")
    print(f"Correctly labeled peptide sequence: {len(Correct)}")
    print(f"Incorrectly labeled peptide sequence: {len(Incorrect)}")
    print(f"Percentage of peptide sequence that has at least one TMT label: {round(At_Least_One / Total * 100, 2)}%")
    print(f"Percentage of correctly labeled peptide sequence: {round(len(Correct) / Total * 100, 2)}%")
    print(f"Percentage of incorrectly labeled peptide sequence: {round(100 - len(Correct) / Total * 100, 2)}%")
    print("**************************** End *****************************\n")

def Save_Sequence_List(Correct_List, Incorrect_List):
    with open('Correctly_labeled_sequences.txt', 'w') as f:
        for item in Correct_List:
            f.write("%s\n" % item)
    with open('Incorrectly_labeled_sequences.txt', 'w') as f:
        for item in Incorrect_List:
            f.write("%s\n" % item)

# call function
# path
is_running = True
while is_running:
    # initiate empty lists
    correctly_labeled = []
    incorrectly_labeled = []

    # user input
    File_location = input("\nWelcome! Please drag your peptide table txt file to your terminal: ")

    # ask if you want to calculate by fraction
    is_by_fraction = input("Do you want to calculate the percentage based on fractions? Type Y/N: ").lower()

    # function calls
    [DF, peptide_col_index, sequence_list, total_num_of_sequence, number_of_fraction] = Read_table(File_location)


    if is_by_fraction == 'y':    # PSM is by fraction
        DF_report = pd.DataFrame({'Fraction': [],
                                  'Total number of peptides': [],
                                  "% of at least one label": [],
                                  "% of correct": [],
                                  "% of incorrect": []})
        for i in number_of_fraction:
            correctly_labeled = []
            incorrectly_labeled = []

            DF_by_fraction = DF.loc[DF['Fraction'] == i]
            sequence_list_by_fraction = DF_by_fraction.iloc[:,peptide_col_index].to_list()
            num_of_sequence_by_fraction = len(sequence_list_by_fraction)
            Check_Sequences(sequence_list_by_fraction)
            at_least_one_label = At_Least_One_Label(sequence_list_by_fraction)

            dict = {'Fraction': i,
                    'Total number of peptides': num_of_sequence_by_fraction,
                    "% of at least one label": round(at_least_one_label / num_of_sequence_by_fraction * 100, 2),
                    "% of correct": round(len(correctly_labeled) / num_of_sequence_by_fraction * 100, 2),
                    "% of incorrect": round(len(incorrectly_labeled) / num_of_sequence_by_fraction * 100, 2)}

            DF_report = DF_report.append(dict, ignore_index=True)
            Print_Report_By_Fraction(i, correctly_labeled, incorrectly_labeled, at_least_one_label, num_of_sequence_by_fraction)
        DF_report.to_csv("report_by_fraction.csv")
        print("Summary report is saved where your main.py is located.")

        correctly_labeled = []
        incorrectly_labeled = []
        Check_Sequences(sequence_list)

        is_saving = input(
            "Do you want to save your correctly and incorrectly labeled sequences into txt files? Type Y/N: ").lower()
        if is_saving == 'y':
            Save_Sequence_List(correctly_labeled, incorrectly_labeled)
            print("\nThe files are saved where your main.py is located.")
        else:
            pass

    else:           # PSM is NOT by fraction
        correctly_labeled = []
        incorrectly_labeled = []
        correctly_labeled = Check_Sequences(sequence_list)
        at_least_one_label = At_Least_One_Label(sequence_list)
        [total_possible_label_residues, total_correctly_labeled_residues, correctly_labeled_N_terminal, correctly_labeled_K_residues
            ] = Amino_Residues_Label(sequence_list, correctly_labeled)
        Print_Report(correctly_labeled, incorrectly_labeled, at_least_one_label, total_num_of_sequence, total_possible_label_residues, total_correctly_labeled_residues,
                     correctly_labeled_N_terminal, correctly_labeled_K_residues)


        is_saving = input("Do you want to save your correctly and incorrectly labeled sequences into txt files? Type Y/N: ").lower()
        if is_saving == 'y':
            Save_Sequence_List(correctly_labeled, incorrectly_labeled)
            print("\nThe files are saved where your main.py is located.")
        else:
            pass

    # ask if you want to continue checking
    is_continue = input("\nDo you want to check your next peptide table? Type Y/N: ").lower()
    if is_continue == 'y':
        is_running = True
    else:
        print("\nHave a nice day :)")
        break
