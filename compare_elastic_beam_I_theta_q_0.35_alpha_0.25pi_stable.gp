reset
# Set the size of the plot
set terminal pngcairo size 1200,1200


# Give the value to the q, alpha and initial values
q=0.35
alpha=0.25*pi
initial1=3.14
initial2=-0.10
initial3=0.00
initial4=-3.24
    
    
# Loop over all the files 
do for [i=0:99] {

# Set the output file 
set output sprintf("compare_combine_elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_stable_%d.png", q, alpha/pi, i)

# Set up the multiplot layout
set multiplot layout 2,2


# Plot the first subplot
set xlabel 'I'
set ylabel "{/symbol q}_{eq}"
set xrange [0:0.001]
set yrange [-4:1]
set title "Fixed Points"

set grid back
    
# Set file handle
file1 = sprintf("elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f.dat", q, alpha/pi, initial1)
file2 = sprintf("elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f.dat", q, alpha/pi, initial2)
file3 = sprintf("elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f.dat", q, alpha/pi, initial3)
file4 = sprintf("elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f.dat", q, alpha/pi, initial4)
    
# Define label size and offset
label_size = 15
label_offset_x = 1
label_offset_y = 0.7

# Extract I value from the file
I_value = real(system(sprintf("awk 'NR==%d+1{print $1}' \"%s\"", i, file1)))

set label 1 sprintf("I = %.9f", I_value) at screen 0.97, 0.01 right front font ",12"

# Plot the picture
plot file1 using 1:(strcol(3) eq "#" ? 1/0 : $3) with points pointsize 0.5 pointtype 7 lc rgb "red" title sprintf('q=%.2f, {/symbol a}=%.3fπ', q, alpha/pi), \
     file2 using 1:(strcol(3) eq "#" ? 1/0 : $3) with points pointsize 0.5 pointtype 7 lc rgb "red" title sprintf('q=%.2f, {/symbol a}=%.3fπ', q, alpha/pi) , \
     file3 using 1:(strcol(3) eq "#" ? 1/0 : $3) with points pointsize 0.5 pointtype 7 lc rgb "orange" title "" , \
     file4 using 1:(strcol(3) eq "#" ? 1/0 : $3) with points pointsize 0.5 pointtype 7 lc rgb "orange" title "" , \
     file1 using 1:(strcol(3) eq "#" ? 1/0 : $3) every ::i::i with points pointsize 1.2 pointtype 7 lc rgb "blue" title "point with iteration (upper branch)", \
     file1 using 1:3:(stringcolumn(8)) every ::i::i with labels offset label_offset_x,label_offset_y font ",label_size" tc rgb "blue" title "" , \
     file3 using 1:(strcol(3) eq "#" ? 1/0 : $3) every ::i::i with points pointsize 1.2 pointtype 7 lc rgb "black" title "point with iteration (lower branch)", \
     file3 using 1:3:(stringcolumn(8)) every ::i::i with labels offset label_offset_x,label_offset_y font ",label_size" tc rgb "black" title ""


# Plot the second subplot

set grid back

# Set the plot range
set xrange [-1:0.3]
set yrange [-0.6:0.7]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Beam with two arms"

# Assign file to the file handle 
file5 = sprintf("beam_first_arm_initial_%.2f_%d.dat",initial1, i)
file6 = sprintf("beam_second_arm_initial_%.2f_%d.dat",initial1, i)
file7 = sprintf("beam_first_arm_initial_%.2f_%d.dat",initial3, i)
file8 = sprintf("beam_second_arm_initial_%.2f_%d.dat",initial3, i)


# Plot the figure 
plot file5 using (-column(3)):(-column(4)) with points pointsize 1.5 pointtype 8 lc rgb "blue" title "elastic boomerang (upper branch) rotates -π",\
     file6 using (-column(3)):(-column(4)) with points pointsize 1.5 pointtype 8 lc rgb "blue" title "" ,\
     file7 using 3:4 with lines linewidth 2 lc rgb "black" title "elastic boomerang (lower branch)",\
     file8 using 3:4 with lines linewidth 2 lc rgb "black" title ""
     
     
# Plot the third subplot     
     
# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set the plot range
set xrange [-0.3:1]
set yrange [-0.6:0.7]

# Set title for the plot
set title "Plot of Beam with two arms (upper branch)"
   
# Since the magnitude of the traction is too small, multiplying a scale so that it can be shown clearly
scale=10.0

# Plot the figure   
plot file5 using 3:4 with lines linewidth 2 lc rgb "blue" title "elastic boomerang (upper branch)",\
     file6 using 3:4 with lines linewidth 2 lc rgb "blue" title "" ,\
     file5 every 5::1 using 3:4:(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 1.2 lc "red" title "10*traction f" , \
     file6 every 15::1 using 3:4:(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 1.2 lc "red" title ""
     
     
# Plot the fourth subplot
# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set the plot range
set xrange [-1:0.3]
set yrange [-0.6:0.7]

# Set title for the plot
set title "Plot of Beam with two arms (lower branch)"

# Plot the figure 
plot file7 using 3:4 with lines linewidth 2 lc rgb "black" title "elastic boomerang (lower branch)",\
     file8 using 3:4 with lines linewidth 2 lc rgb "black" title "" ,\
     file7 every 5::1 using 3:4:(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 1.2 lc "orange" title "10*traction f" , \
     file8 every 15::1 using 3:4:(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 1.2 lc "orange" title ""


# Remove label for the next frame
unset label 1

unset multiplot
}

# Save the final multiplot layout to the output file
set output
