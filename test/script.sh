#!/bin/bash

../bin/linux64/bayescenv sim_small.txt -env env.txt -out_pilot -all_trace -unif_pr 10 -pr_jump 0.1 -pr_pref 0.5 -o test -burn 3000 -thin 10 -nbp 10 -pilot 2000
