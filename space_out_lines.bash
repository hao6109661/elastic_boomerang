#! /bin/bash

# Check if the number of arguments passed to the script is not equal to 1
if [ $# -ne 1 ]; then
    echo "Need filename as single argument" # If not, print an error message
    exit # Exit the script
fi

echo "Filename = " $1 # Print the filename provided as an argument

# AWK command to process the content of the file provided as an argument
awk 'BEGIN{tol=0.1; i=0}{I_current=$1; theta_current=$3; if (i>1){ds = sqrt((I_current-I_prev)*(I_current-I_prev)+(theta_current-theta_prev)*(theta_current-theta_prev)); if (ds > tol){print "\n"$0}else{print $0} }else{print $0};
theta_prev=theta_current; I_prev=I_current; i++}' $1
