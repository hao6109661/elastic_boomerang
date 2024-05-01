#! /bin/bash

if [ $# -ne 1 ]; then
    echo "Need filename as single argument"
    exit
fi

echo "Filename = " $1


awk 'BEGIN{tol=0.01; i=0}{I_current=$1; theta_current=$3; if (i>1){ds = sqrt((I_current-I_prev)*(I_current-I_prev)+(theta_current-theta_prev)*(theta_current-theta_prev)); if (ds > tol){print "\n"$0}else{print $0} }else{print $0};
theta_prev=theta_current; I_prev=I_current; i++}' $1
