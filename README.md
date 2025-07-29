# Introduction
This version of rintc is software designed to estimate energy parameters for modified base substructures. 
It performs parameter estimation by taking ground truth base pairing probabilities and the parameters you wish to estimate as input.

# Technologies Used
The source code is written in C++. It can be built using a Makefile, and the g++ compiler is utilized. 
The compiler version we used was GCC 11.4.0.

# build
We've provided a Makefile, so you can build the executable rintc file using the make command. 
We compiled it on Ubuntu 22.04 within WSL V.2.

# sample
A run.sh file is provided as an example.
When you run this script, "rintc" will start in the background to estimate the energy parameters of modified bases using gradient descent.


# Argument Descriptions

#arg1:"train" only
Currently, only "train" is accepted.

#arg2:nsample
This argument is currently inactive. You can skip it by providing any value.

#arg3:seq_size
This argument is currently inactive. You can skip it by providing any value.

#arg4:max_span
This argument is currently inactive. You can skip it by providing any value.

#arg5:learning_rate
Specifies the learning rate. Please adjust it based on the noise level and the amount of data.

#arg6:epoch num
Specifies the number of epochs for training.

#arg7:energy_param initial value
Specifies the initial value for the energy parameters you want to estimate. Provide the value in 100×kcal/mol. For example, -100 corresponds to -1.00 kcal/mol.

#arg8: weight_central
This argument is currently fixed at 1.0.

#arg9:param_file:correct energy_param
Specifies the file containing the energy parameters. For the parameters you want to estimate, please enter any arbitrary value. We've provided the following sample file for reference:
energy_param/modified_nuculeobases.m6a_20240718.par
It uses the same format as ViennaRNA's .par files. Parameters marked with a 1 in the update_mask.txt file will be subject to estimation.

#arg10:data_file: .bppm file (seq and correct bppm)
Specifies the file containing the Ground Truth base pairing probabilities. We've provided the following sample file for format reference:
data/m6A_random_135.bppm

# Update Mask Settings
You need an update_mask.txt file in the same directory as the rintc executable. This file specifies which parameters are subject to updates. 0 indicates the parameter is not an update target, while 1 indicates it is an update target.

# Interpreting the Results
Parameters during optimization are output to standard output. We recommend redirecting this output to a log file.
The parameter table is printed after each [logstack parameter] line, every time an update calculation is performed. This data is output with every calculation, even if no parameters are updated as a result of the calculation.
Therefore, you can track the update status by parsing and plotting this table with a log analysis script, along with the number of updates.

# How to Calculate Base Pair Probabilities
If you want to calculate the base pair probabilities, specify "calc" as the first argument to rintc. For the second argument, provide the FASTA file containing the sequences you want to calculate.
The calculation results will be saved in a new file where the .fa extension is replaced with .bppm.

Example: ./rintc calc random.fa

The resulting base pair probability file can then be directly used for rintc's training.