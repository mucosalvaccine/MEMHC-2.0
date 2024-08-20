## License

This software is licensed under the [CC BY-NC 4.0 License](LICENSE) for non-commercial use.
2019

For commercial use, please contact Mohammad Arabpour at mohammad.arabpour.sinior@gmail.com to obtain a commercial license.

"""
The script provided takes in a user input of a protein sequence, minimum and maximum peptide lengths, and a list of HLA alleles from a csv file. It generates peptides of different lengths from the protein sequence, and then predicts the binding affinity of each peptide with each HLA allele. The predicted binding affinities are stored in a pandas DataFrame and filtered to remove peptides that are similar to each other based on their start positions. Lastly, the filtered DataFrame is converted to binary values based on a Kd_treshold of 500 nM and the relative HLA coverage is calculated.

To execute the script, you will need to install the mhcflurry package, which is used to predict binding affinities. You can install it by running pip install mhcflurry. You will also need to have pandas installed, which you can install by running pip install pandas.

Here is what this code do:

Import the necessary modules and packages, including pandas and mhcflurry, NetMHC, etc,.
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

Work procedure ;

This script is multi potential and could be ran by diffrent MHC affinity predicting based on machin learning. Here we use mhcflurry as an example, but it could be used by any of the following : 

MHCflurry
NetMHC
NetMHCpan
NetMHCII
NetMHCIIpan
DeepHLApan
MixMHCpred
MARIA
MHC-Bench
HLA-CNN
Depen

Running Script:

For running this script you will need to install mhcflurry.
you would then need to input a protein sequence of intrest, a minimum and a maximum lenght for the peptide in comand-prompt . 
You will also need to give a HLA file including all the inquire HLA that are intended for prediction. An HLA file with exhusative list of HLA (MHCI) is enlisted in HLA type I.csv file as an example.
An out put would be a CSV file (MHCI-output-minimal-epitop-maxcoverage.csv)


Please, refer to Author (Mohammd Arabpour) for any further questions.











