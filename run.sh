#!/bin/bash -eux

mkdir logs -p

#arg1:"train" only
#arg2:nsample(no use)
#arg3:seq_size(no use)
#arg4:max_span(no use)
#arg5:learning_rate
#arg6:epoch num
#arg7:energy_param initial value
#arg8: weight_central(?)
#arg9:param_file:correct energy_param
#arg10:data_file: .bppm file (seq and correct bppm)



./rintc train 100 200 50 50.0 300 -100.0 1.0 modified_nuculeobases.m6a_20240718 data/m6A_random_135.bppm > data/m6A_random135_00_new_gyaku.log &

