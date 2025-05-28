#rintc Version for estimating energy parameters of modified bases using gradient descent.
$ make
This operation compiles the necessary components to create the executable "rintc."
The compilation was successful on Ubuntu 22.04 within WSL V2.


./run.sh
When you run this script, "rintc" will start in the background to estimate the energy parameters of modified bases using gradient descent.

Here's a brief explanation of the "rintc" arguments:
#arg1:"train" only
#arg2:nsample
#arg3:seq_size
#arg4:max_span
#arg5:learning_rate
#arg6:epoch num
#arg7:energy_param initial value
#arg8: weight_central
#arg9:param_file:correct energy_param
#arg10:data_file: .bppm file (seq and correct bppm)
