#!/bin/bash

q=0.35
alpha_without_pi=0.250


# Set frame rate and input image file name pattern
framerate=9

# Combine the first set of images into the first MP4 animation
ffmpeg -framerate $framerate -i "compare_combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_stable_%d.png" -c:v libx264 -r $framerate -pix_fmt yuv420p first.mp4
#ffmpeg -framerate $framerate -i "combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_initial_${initial1}_I_max_0.1_%d.png" -c:v libx264 -r $framerate -pix_fmt yuv420p first.mp4
#ffmpeg -framerate $framerate -i "compare_combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_%d.png" -c:v libx264 -r $framerate -pix_fmt yuv420p first.mp4

# Combine the second set of images into the second MP4 animation
ffmpeg -framerate $framerate -i "compare_combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_unstable_%d.png" -c:v libx264 -r $framerate -pix_fmt yuv420p second.mp4
#ffmpeg -framerate $framerate -i "combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_initial_${initial2}_I_max_0.1_%d.png" -c:v libx264 -r $framerate -pix_fmt yuv420p second.mp4
#ffmpeg -framerate $framerate -i "sec_combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_%d.png" -c:v libx264 -r $framerate -pix_fmt yuv420p second.mp4

# Create a file list containing the names of the two MP4 files
echo "file 'first.mp4'" > filelist.txt
echo "file 'second.mp4'" >> filelist.txt

# Merge both MP4 animations into one MP4 file
#ffmpeg -f concat -safe 0 -i filelist.txt -c copy combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_four_branches.mp4
#ffmpeg -f concat -safe 0 -i filelist.txt -c copy combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_I_max_0.1.mp4
ffmpeg -f concat -safe 0 -i filelist.txt -c copy compare_combine_elastic_beam_I_theta_q_${q}_alpha_${alpha_without_pi}pi_four_branches.mp4

# Clean up intermediate files
rm first.mp4 second.mp4 filelist.txt
