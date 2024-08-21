"""
The script provided takes in a user input of a protein sequence, minimum and maximum peptide lengths, and a list of HLA alleles from a csv file. It generates peptides of different lengths from the protein sequence, and then predicts the binding affinity of each peptide with each HLA allele. The predicted binding affinities are stored in a pandas DataFrame and filtered to remove peptides that are similar to each other based on their start positions. Lastly, the filtered DataFrame is converted to binary values based on a Kd_treshold of 500 nM and the relative HLA coverage is calculated.

To execute the script, you will need to install the mhcflurry package, which is used to predict binding affinities. You can install it by running pip install mhcflurry. You will also need to have pandas installed, which you can install by running pip install pandas.

Here is a step-by-step explanation of the code:

Import the necessary modules and packages, including pandas and mhcflurry.
Define a function to generate peptides of different lengths from a given protein sequence. The function takes in the protein sequence, minimum and maximum peptide lengths, and returns a list of peptides, their start positions, and their lengths.
Get user inputs for the protein sequence, minimum and maximum peptide lengths.
Check if the protein sequence contains only accepted amino acid characters and if the minimum and maximum length inputs are integers.
Convert the minimum and maximum length inputs to integers.
Generate peptides of different lengths from the protein sequence using the generate_peptides function.
Load the Class1AffinityPredictor from the mhcflurry package, which will be used to predict binding affinities.
Read in a csv file containing a list of HLA alleles, and store them in a list.
Loop through each HLA allele and each peptide, and predict the binding affinity of the peptide with the HLA allele using the predict method of the Class1AffinityPredictor. Store the predicted binding affinities in a list.
Create a pandas DataFrame from the predicted binding affinities, and pivot the DataFrame so that each row corresponds to a peptide and each column corresponds to an HLA allele.
Group the rows of the DataFrame by the start position of the peptide.
Loop through each group, and filter out peptides that are similar to each other based on their start positions. Keep only the peptide with the minimum binding affinity for each HLA allele.
Convert the binding affinities to binary values based on a Kd_treshold of 500 nM, and calculate the relative HLA coverage, which is the sum of the binary values in each column divided by the total number of peptides.
Select the peptide with maximum coverage score and eliminate the HLA alleles correspondent to those ones
loop though the rest of peptides to find next maximum coverage score untill the score is 100% coverage
Make a dataframe out of peptides and their correspondent coverge score
export the df into csv format
"""
#import packages
import pandas as pd
from mhcflurry import Class1AffinityPredictor
import numpy as np

##Making protein into peptide of different length
def generate_peptides(protein_sequence, min_length, max_length):
    peptides = []
    start_positions = []
    peptide_length = []
    for i in range(len(protein_sequence)):
        for j in range(i+min_length, min(i+max_length+1, len(protein_sequence)+1)):
            peptides.append(protein_sequence[i:j])
            peptide_length += [len(peptides)]
            start_positions.append(i)
    return peptides, start_positions, peptide_length

# Define the set of accepted amino acid characters
accepted_aas = set("ACDEFGHIKLMNPQRSTVWY")

# Get user inputs for protein sequence, minimum length, and maximum length
protein_sequence = input("Please enter your protein sequence here: ")
protein_sequence = protein_sequence.upper()
print(protein_sequence)
min_length = input("Please enter the minimum peptide length: ")
max_length = input("Please enter the maximum peptide length: ")

# Check if protein sequence contains only accepted amino acid characters
for i in protein_sequence:
    if set(protein_sequence.upper()) - accepted_aas:
        print("Warning: Your protein sequence contains invalid amino acid characters")
    break

# Check if minimum and maximum length inputs are integers
for i in min_length or max_length:
    if not (min_length.isdigit() and max_length.isdigit()):
     print("Warning: Minimum and maximum length inputs must be integers")
    break

# Convert minimum and maximum length inputs to integers
min_length = int(min_length)
max_length = int(max_length)

# Rest of the code
peptide, start , pe= generate_peptides(protein_sequence, min_length, max_length)

# print("peptides:", peptide)
# print("Start :", start)


# predict binding affinity throughout different HLA allele
predictor = Class1AffinityPredictor.load()

HLA_type_I = pd.read_csv("HLA type I.csv")
HLA_type_I_list = list(HLA_type_I.columns)
peptide_list = list(peptide)

kd_predictions = []
peptide_length = []
start_positions = []

for mhc_allele in HLA_type_I_list:
    for pep in peptide_list:
        peptide_length.append(len(pep))
        start_positions.append(start[peptide.index(pep)])
        try:
            predicted_affinity = predictor.predict([pep], [mhc_allele])
            kd_predictions += list(predicted_affinity)
        except ValueError:
            kd_predictions += ["NA"]
            continue

# Create DataFrame
df = pd.DataFrame(kd_predictions, index=pd.MultiIndex.from_product([peptide_list, HLA_type_I_list]),
                  columns=["Binding Affinity"])
df.index.names = ["Peptide", "MHC Allele"]
df = df.reset_index()
df = df.pivot(index="Peptide", columns="MHC Allele", values="Binding Affinity")
# df.insert(0, "id", range(1, len(df)+1))
# df.insert(1, "peptide_length", [len(p) for p in peptide_list])
# df.insert(2, "start", start_positions)
df["peptide"] = peptide
df["start"] = start
#df["peptide_enght"] =peptide_length
#df1["peptide_lenght"] = peptide_lenght
df1 = df.reindex(columns=["peptide", "start"] + list(df.columns[:-2]))
df1.to_csv("binding_affinities.csv", index=False)

df1 = pd.read_csv("binding_affinities.csv")
# Group rows by the same value in the 'start' column
grouped = df1.groupby(['start'])
print(grouped)
# Create an empty dataframe to store the filtered values
filtered_df = pd.DataFrame(columns=df1.columns)

# Loop through each group
for group_name, group_data in grouped:
    # Get the minimum value of each column (excluding peptide_length, start, and id)
    min_values = group_data.drop([ 'start',  'peptide'], axis=1).min()
    # Create a new dataframe with the same columns as group_data
    filtered_group = pd.DataFrame(columns=df1.columns)
    # Copy the peptide_length, start, and id columns to the filtered_group dataframe
    filtered_group['start'] = group_data['start']
    filtered_group['peptide'] = group_data['peptide']
    # Set values greater than the minimum to 1000
    for col in min_values.index:
        filtered_group[col] = group_data[col].apply(lambda x: x if x == min_values[col] else 1000)
    # Append the filtered_group dataframe to the filtered_df dataframe
    filtered_df = pd.concat([filtered_df, filtered_group], ignore_index=True)
df1 = filtered_df
print(df1)


Kd_treshold_nM = 500
exclude_cols = ["start", "peptide"]  # list of column names to exclude
numeric_cols = df1.columns.drop(exclude_cols)  # list of numeric column names

# convert numeric columns to numeric dtype
df1[numeric_cols] = df1[numeric_cols].apply(pd.to_numeric, errors="coerce")
# Convert Kd_treshold_nM to 1 and 0, excluding NaN values

# apply the condition to the binary numeric cells where it is  a number
#df1[numeric_cols] = np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0)
#df1[numeric_cols] = np.where(np.isan(df1[numeric_cols]) & np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0))
df1[numeric_cols] = np.where(~np.isnan(df1[numeric_cols]), np.where(df1[numeric_cols] <= Kd_treshold_nM, 1, 0), np.nan)
print(df1)  # print the updated dataframe

df1.to_csv("filterd for kd.csv", index=False)

# Initialize the total coverage to 0
relative_HLA_coverage = 0
# Loop through each column (excluding peptide_length, start, id, and peptide)
# Define the columns to exclude from the loop
exclude_cols = ["start", "peptides"]
# Initialize the relative_HLA_coverage
relative_HLA_coverage = 0

# Loop through each column (excluding peptides, start, id, and peptide_length)
for col in df1.columns[2:]:
    if col not in exclude_cols:
        # Sum the values of the column and set the sum value to 1 if the sum is greater than 0, and 0 otherwise
        sum_value = df1[col].sum()
        if sum_value > 0:
            sum_value = 1
        else:
            sum_value = 0

        # Add the sum_value to the total_coverage
        relative_HLA_coverage += sum_value

print("Relative HLA coverage:", relative_HLA_coverage)
df1.to_csv('filtered_e7_netmhcpan_2915_alleles_8_to_11_mers.csv', index=False)

#make the sume coverage value column for each peptide
df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
#print(df1['sum_weight'])
df1.to_csv('E7_clement_result_table.csv', index=False)
total1 = (df1.shape[1])-1
#print(total1)
df1 = df1.sort_values(by='sum_weight',axis=0, ascending =False, inplace=False, kind='quicksort', na_position='last', ignore_index=False, key=None)
maximum1 = max(df1['sum_weight'])
lead_peptide = (df1.iloc[0,0])
Lead_p = len(df1.iloc[0,0])

leadcoverage = (maximum1 / relative_HLA_coverage) * 100
row_number = df1[df1['sum_weight'] == maximum1].index
df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
# print ('totalscore',total1)
# print('lead pep score',df1.iloc[0]['sum_weight'])
# print('lead pep',lead_peptide)
# print('leadcoverage',leadcoverage)
# print('rowindex',row_number)
# print(df1)
totalcoverage = leadcoverage
total = total1
df1.to_csv('jojojo.csv', index=False)


####################bringing start to the second column
col_name = 'start' # Replace 'column_name' with the name of the column you want to move

# Remove the column from its current position and save it to a variable
col_to_move = df1.pop(col_name)

# Insert the column at the second position (index 1)
df1.insert(1, col_name, col_to_move)

print (df1)
max_num = maximum1
#print(maximum1)
#####going to second line
#####------------------------------------------------------------------------------------------------------------------
totalscore =[total1]
#lead_pep_score = [df1.iloc[0]['sum_weight']]
lead_pep_score = [max_num]
lead_pep = [lead_peptide]
Lead_peptide_lenght = [int(Lead_p)]
rowindex = [row_number]
accuulativcoverage = [leadcoverage]
start = [(df1.iloc[0, 1])]
start_data = df1['start']  # Extract the data from the start column
df1['start'] = df1['start'].astype(str)
while total > 0:
    df= df1.iloc[0]
    s = ((df.index))
    g = []
    for i in range (0, len(s)):
        h = df.iloc[i]
        if type(h) != str:
            if h == 0:
             g = [s[i]] + g
    ali = ['peptide', 'start'] + g
    df1 = df1[ali]
    df1.shape[0]

    df1 = df1[1:(df1.shape[0])]
    df1['sum_weight'] = df1.sum(axis=1, skipna=True, numeric_only=True)
    total = df1['sum_weight'].sum()
    #print(total)
    df1 = df1.sort_values(by='sum_weight',axis=0, ascending =False, inplace=False, kind='quicksort', na_position='last', ignore_index=False, key=None)
    maximum = max(df1['sum_weight'])
    max_num += maximum
    lead_peptide = (df1.iloc[0,0])

    Lead_peptide_lenght += [len(df1.iloc[0, 0])]
    leadcoverage = [(maximum / total1) * 100]
    #print (max_num)
    row_number = df1[df1['sum_weight'] == maximum].index
    #calculate accumulative coverage from all the HLA allotypes (including those that are not recognized by any peptide in peptide pool)
    absolut_acumulativelcoverage = (max_num/total1) * 100
    #calculate acumulative coverage from all the recognized hits of HLA allotypes (excluding those that are not recognized by any peptide in peptide pool)
    acumulativelcoverage = (max_num/relative_HLA_coverage) * 100
    #print(max_num)
    # print('rowindex', row_number[0])
    # print('lead pep score', df1.iloc[0]['sum_weight'])
    # print('lead pep', lead_peptide)
    # print('leadcoverage', leadcoverage)
    # print ('accuulativcoverage', acumulativelcoverage)
    # print(df1)
    accuulativcoverage += [acumulativelcoverage]
    lead_pep_score += [df1.iloc[0]['sum_weight']]

    #print(lead_pep_score)
    lead_pep += [lead_peptide]
    rowindex += [(row_number[0])+1]
    start += (df1.iloc[0, 1])

# Combine start_data and the results of the while loop into a single DataFrame
data = { "start position in protein" : start,"Peptide lenght": Lead_peptide_lenght, 'peptide sequence': lead_pep,'# HLA allel hits': lead_pep_score, 'accumulativcoverage(%)': accuulativcoverage }
df3 = pd.DataFrame(data)

print (df3)
df3.to_csv('MHCI-output-minimal-epitop-maxcoverage.csv', index=False)
