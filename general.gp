reset
# Set the size of the plot
set terminal pngcairo size 1200,1200

set xlabel font 'Arial,21'
set ylabel font 'Arial,17'
set title font 'Arial,17'

number_q=2


initial1=-4.80

array step[2]
step[1]=36
step[2]=39

# Set the output file 
set output sprintf("RESLT_gernal_plot/I_plot.png")

# Set up the multiplot layout
set multiplot

# Loop over different I value within the interval [0, I_max], I_max is "max_value" in the previous text
do for [k=0:1] {

I=k*0.1

# Set the output file 
#set output sprintf("RESLT_gernal_plot/I_%.2f.png", I)

filelist = ""
 do for [z=1:number_q] {
     # Create file names using suffix and add them to filelist
      handle = sprintf("data%d.dat", z)
      filelist = filelist . "\"" . handle . "\" "}
      
# print filelist

set size 0.35,0.35
set origin 0+k*0.3,0


do for [i=1:number_q] {

q=0.38+(i-1)*0.02

filename = sprintf("data%d.dat", i)
set print filename

do for [j=1:step[i]] {

# Give the value to the q,alpha and initial value
alpha=((1.8*j)/180.0)*pi

# Set file handle for the general information which contains the I, theta_eq and step number 
file1 = sprintf("../elastic_boomerang_animation/RESLT_q_%.3f_alpha_%.3f_element_20/elastic_beam_I_theta_q_%.3f_alpha_%.3fpi_initial_%.2f.dat", q, 1.8*j, q, alpha/pi, initial1)

# Got the max value of I from the general file 
max_value = system(sprintf("awk '{if(NR==1 || $1 > max) max=$1} END {print max}' \"%s\"", file1))
# print "Max value in the first column:", max_value


if(max_value>=I){
  print q, alpha
}

}
unset print
}

# Set the plot range
set xrange [0:0.5]
set yrange [0:pi]
set grid back

set ytics( "0" 0 , "0.125π" 0.125*pi, "0.25π" 0.25*pi , "0.375π" 0.375*pi, "0.5π" 0.5*pi , "0.625π" 0.625*pi, "0.75π" 0.75*pi , "0.875π" 0.875*pi, "π" pi)

# Set labels for the axes
set xlabel "q"
set ylabel "{/symbol a}"

# Set title for the plot
set title "Plot of Fixed Points"

plot for [file in filelist] \
             file using 1:2 with points title ""
     

             

# Delte the temporary files
#do for [i=1:number_q] {
    #system(sprintf("rm data%d.dat", i))
#}


}









##############################################################################################################################################################

##############################################################################################################################################################

# Exit the mode of multiplot
unset multiplot

# Save the final multiplot layout to the output file
set output
