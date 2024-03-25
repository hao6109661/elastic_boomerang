reset

# Set the format of the output file 
set term png

set grid back

# Normally, we use beam0.dat as the reference because there is no traction applied to the beam
file_ref=sprintf("beam0.dat")

# Loop over all files 
do for [i=1:1]{

# Non-dimensional coefficient (FSI)
# Note that this is not the real value of Q, instead the real one is "Q=i*10^{-7}"
Q=i

# Set the output file 
# Note that the output files will not pop up, instead they will be saved in the directory automatically 
set output sprintf("beam_with_two_arm%d.png",i)

# Set the plot range
    set yrange [-1:5]
    set xrange [-2:6]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Traction"

# Assign file to the file handle 
file = sprintf("beam%d.dat",i)
file2 = sprintf("beam_second_arm%d.dat",i)

# Write the additional information to the plot 
#set label 1 sprintf("{/symbol q} = π/6,V = 0.5, U_0 = 0.8, X_0 = 1.5, Y_0 = 2.5") at screen 0.98, 0.02 right front font ",10"
set label 2 sprintf("Q=%.1f*10^{-7}",Q) at screen 0.2, 0.02 right front font ",12"

# Plot the figures 
plot file using 1:2 with lines linewidth 2 lc rgb "red" title "Beam", \
     file using 3:4 with lines linewidth 2 lc rgb "blue" title "Beam move",\
     file2 using 1:2 with lines linewidth 2 lc rgb "pink" title "Beam 2",\
     file2 using 3:4 with lines linewidth 2 lc rgb "brown" title "Beam move 2"
}

