reset
# Set the size of the plot
set terminal pngcairo size 1200,1200

set xlabel font 'Arial,14'
set ylabel font 'Arial,14'
set title font 'Arial,14'

number_q=11


initial1=-4.80

array step[number_q]
step[1]=5
step[2]=9
step[3]=18
step[4]=24
step[5]=29
step[6]=33
step[7]=36
step[8]=39
step[9]=41
step[10]=44
step[11]=46

# Set the output file 
set output sprintf("RESLT_gernal_plot/I_plot_1.png")

# Set up the multiplot layout
set multiplot


# Loop over different I value within the interval [0, I_max], I_max is "max_value" in the previous text
do for [k=0:9] {

I=k*0.0001

# Set the output file 
#set output sprintf("RESLT_gernal_plot/I_%.2f.png", I)

filelist = ""
 do for [z=1:number_q] {
     # Create file names using suffix and add them to filelist
      handle = sprintf("data%d.dat", z)
      filelist = filelist . "\"" . handle . "\" "}
      
# print filelist


do for [i=1:number_q] {

if(i==1){
q=0.275
}
else{
q=0.28+(i-2)*0.02
}

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

set ytics( "0" 0 , "" 0.125*pi, "0.25π" 0.25*pi , "" 0.375*pi, "0.5π" 0.5*pi , "" 0.625*pi, "0.75π" 0.75*pi , "" 0.875*pi, "π" pi)

# Set labels for the axes
set xlabel "q"
set ylabel "{/symbol a}"
set title sprintf("𝓘=%.4f",I)

set size 0.237,0.20
set origin 0,0.58

if(k>0){
unset tics
unset xlabel
unset ylabel
set size 0.175,0.17
set origin k*0.086+0.061,k*0.024+0.61
}

plot for [file in filelist] \
             file using 1:2 with points pointtype 7 pointsize 1 lc "red" title ""
     
    
# Delte the temporary files
do for [i=1:number_q] {
    system(sprintf("rm data%d.dat", i))
}


}




#######################################################################################################################
# Loop over different I value within the interval [0, I_max], I_max is "max_value" in the previous text
do for [k=10:19] {

I=k*0.0001

# Set the output file 
#set output sprintf("RESLT_gernal_plot/I_%.2f.png", I)

filelist = ""
 do for [z=1:number_q] {
     # Create file names using suffix and add them to filelist
      handle = sprintf("data%d.dat", z)
      filelist = filelist . "\"" . handle . "\" "}
      
# print filelist


do for [i=1:number_q] {

if(i==1){
q=0.275
}
else{
q=0.28+(i-2)*0.02
}

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

set title sprintf("𝓘=%.4f",I)

minus=0.2

set origin (k-10)*0.086+0.061,(k-10)*0.024+0.61-minus


plot for [file in filelist] \
             file using 1:2 with points pointtype 7 pointsize 1 lc "red" title ""
     
    
# Delte the temporary files
do for [i=1:number_q] {
    system(sprintf("rm data%d.dat", i))
}


}




#######################################################################################################################
# Loop over different I value within the interval [0, I_max], I_max is "max_value" in the previous text
do for [k=20:29] {

I=k*0.0001

# Set the output file 
#set output sprintf("RESLT_gernal_plot/I_%.2f.png", I)

filelist = ""
 do for [z=1:number_q] {
     # Create file names using suffix and add them to filelist
      handle = sprintf("data%d.dat", z)
      filelist = filelist . "\"" . handle . "\" "}
      
# print filelist


do for [i=1:number_q] {

if(i==1){
q=0.275
}
else{
q=0.28+(i-2)*0.02
}

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

set title sprintf("𝓘=%.4f",I)

minus=0.2*2

set origin (k-20)*0.086+0.061,(k-20)*0.024+0.61-minus


plot for [file in filelist] \
             file using 1:2 with points pointtype 7 pointsize 1 lc "red" title ""
     
    
# Delte the temporary files
do for [i=1:number_q] {
    system(sprintf("rm data%d.dat", i))
}


}





#######################################################################################################################
# Loop over different I value within the interval [0, I_max], I_max is "max_value" in the previous text
do for [k=30:39] {

I=k*0.0001

# Set the output file 
#set output sprintf("RESLT_gernal_plot/I_%.2f.png", I)

filelist = ""
 do for [z=1:number_q] {
     # Create file names using suffix and add them to filelist
      handle = sprintf("data%d.dat", z)
      filelist = filelist . "\"" . handle . "\" "}
      
# print filelist


do for [i=1:number_q] {

if(i==1){
q=0.275
}
else{
q=0.28+(i-2)*0.02
}

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

set title sprintf("𝓘=%.4f",I)

minus=0.2*3

set origin (k-30)*0.086+0.061,(k-30)*0.024+0.61-minus


plot for [file in filelist] \
             file using 1:2 with points pointtype 7 pointsize 1 lc "red" title ""
     
    
# Delte the temporary files
do for [i=1:number_q] {
    system(sprintf("rm data%d.dat", i))
}


}


##############################################################################################################################################################

##############################################################################################################################################################

# Exit the mode of multiplot
unset multiplot

# Save the final multiplot layout to the output file
set output
