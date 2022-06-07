"""
Calculates the percentage of amino acid residues that are labeled relative to the number that should be labeled
"""
import pandas as pd


# functions
def Read_table(File):
    """Load the user input file path"""
    DF = pd.read_table(File)
    sequence_list = DF.iloc[:,0].to_list()
    num_of_sequence = len(sequence_list)

    return sequence_list, num_of_sequence

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
        else:
            incorrectly_labeled.append(string)

def Print_Report(Correct, Incorrect, Total):
    print("\n**************************** Report **************************")
    print(f"Total number of peptide sequences: {Total}")
    print(f"Correctly labeled peptide sequences: {len(Correct)}")
    print(f"Incorrectly labeled peptide sequences: {len(Incorrect)}")
    print(f"Percentage of correctly labeled peptide sequences: {round(len(Correct) / Total * 100, 2)}%")
    print(f"Percentage of incorrectly labeled peptide sequences: {round(100 - len(Correct) / Total * 100, 2)}%")
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

    # function calls
    [sequence_list, total_num_of_sequence] = Read_table(File_location)
    Check_Sequences(sequence_list)
    Print_Report(correctly_labeled, incorrectly_labeled, total_num_of_sequence)

    is_saving = input("Do you want to save your correctly and incorrectly labeled sequences into txt files? Type Y/N: ").lower()
    if is_saving == 'y':
        Save_Sequence_List(correctly_labeled, incorrectly_labeled)
        print("\nThe files are saved where your main.py is located.")
    else:
        pass

    is_continue = input("\nDo you want to check your next peptide table? Type Y/N: ").lower()
    if is_continue == 'y':
        is_running = True
    else:
        print("\nHave a nice day :)")
        break
