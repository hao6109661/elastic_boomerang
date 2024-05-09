#! /bin/bash


# Build the executable
make hao

if [ -e RESLT ]; then
    echo "RESLT directory already exists. Please remove it."
    exit
fi

    # Make the result directory
    mkdir RESLT

    q=0.39
    alpha_in_degrees=22.5

    # run the code
    # --alpha_in_degrees 45  > OUTPUT_alpha$alpha_in_degree
    ./hao --q $q --alpha_in_degrees $alpha_in_degrees --Initial_value_for_theta_eq 0.00 --I_max 0.03 --I_increment_default 1.0e-4

    ./hao --q $q --alpha_in_degrees $alpha_in_degrees --Initial_value_for_theta_eq -0.10 --I_max 0.03 --I_increment_default 1.0e-4

    #./hao --q $q --alpha_in_degrees $alpha_in_degrees --Initial_value_for_theta_eq 0.10 --I_max 2.0 --I_increment_default 1.0e-2

    #./hao --q $q --alpha_in_degrees $alpha_in_degrees --Initial_value_for_theta_eq -1.62 --I_max 2.0 --I_increment_default 1.0e-2

    #--I_max --I_increment_default

    # Move the results
    mv RESLT RESLT_q_${q}_alpha_${alpha_in_degrees}_four_branches
