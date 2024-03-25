reset
set term png

set grid back

file_ref=sprintf("beam0.dat")

do for [i=1:1]{
Q=i

set output sprintf("beam_with_two_arm%d.png",i)
# Set the plot range
    set yrange [-1:5]
    set xrange [-2:6]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Traction"


#set label 1 sprintf("~{/symbol g}{1\\.} =%.2f, Scale=10^{-4}.", (1.0e-2)*i) at screen 0.97, 0.02 right front font ",15"

# Plot the first subplot
file = sprintf("beam%d.dat",i)
file2 = sprintf("beam_second_arm%d.dat",i)

#set label 1 sprintf("{/symbol q} = π/6,V = 0.5, U_0 = 0.8, X_0 = 1.5, Y_0 = 2.5") at screen 0.98, 0.02 right front font ",10"
set label 2 sprintf("Q=%.1f*10^{-7}",Q) at screen 0.2, 0.02 right front font ",12"

# Calculate the endpoint coordinates of the vectors
#set datafile separator ","  # Assuming data file is comma-separated
plot file using 1:2 with lines linewidth 2 lc rgb "red" title "Beam", \
     file using 3:4 with lines linewidth 2 lc rgb "blue" title "Beam move",\
     file2 using 1:2 with lines linewidth 2 lc rgb "pink" title "Beam 2",\
     file2 using 3:4 with lines linewidth 2 lc rgb "brown" title "Beam move 2"
}

