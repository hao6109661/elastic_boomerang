reset

# Set the format of the output file 
set term png

set size ratio -1

set grid back

# Normally, we use beam0.dat as the reference because there is no traction applied to the beam
file_ref1=sprintf("beam_first_arm0.dat")
file_ref2=sprintf("beam_second_arm0.dat")

# Loop over all files 
do for [i=0:10]{

# Non-dimensional coefficient (FSI)
# Note that this is not the real value of Q, instead the real one is "Q=i*10^{-6}"
Q=i

# Set the output file 
# Note that the output files will not pop up, instead they will be saved in the directory automatically 
set output sprintf("beam_with_two_arm%d.png",i)

# Set the plot range
set yrange [-0.5:2.5]
set xrange [-0.5:2.5]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Beam with two arms"

# Assign file to the file handle 
file = sprintf("beam_first_arm%d.dat",i)
file2 = sprintf("beam_second_arm%d.dat",i)

# Write the additional information to the plot 
#set label 1 sprintf("{/symbol q} = π/6,V = 0.5, U_0 = 0.8, X_0 = 1.5, Y_0 = 2.5") at screen 0.98, 0.02 right front font ",10"
#set label 2 sprintf("Q=%.1f*10^{-6}",Q) at screen 0.2, 0.02 right front font ",12"
set label 1 sprintf("Q = %.1f*10^{-4}, {/symbol a} = π/4, q = 0.4, X_0 = 1.0, Y_0 = 1.5", Q) at screen 0.95, 0.02 right front font ",12"

# Plot the figures 
plot file using 1:2 with lines linewidth 2 lc rgb "red" title "first arm", \
     file using 3:4 with lines linewidth 2 lc rgb "blue" title "first arm (move)",\
     file2 using 1:2 with lines linewidth 2 lc rgb "green" title "second arm",\
     file2 using 3:4 with lines linewidth 2 lc rgb "orange" title "second arm (move)"
     #file_ref1 using 1:2 with lines linewidth 1 lc rgb "black" title "reference",\
     #file_ref1 using 3:4 with lines linewidth 1 lc rgb "black" title "reference",\
     #file_ref2 using 1:2 with lines linewidth 1 lc rgb "black" title "reference",\
     #file_ref2 using 3:4 with lines linewidth 1 lc rgb "black" title "reference"
}

