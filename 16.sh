## REARRANGE PLOTS IN 16 PANELS: "convert +/-append" ##
## AND ANIMATE: mencoder ##

# CREDIT: user fmw42 on
# http://www.imagemagick.org/discourse-server/viewtopic.php?t=14743
# convert \( image1 image2 -append \) \
# \( image3 image4 -append \) +append result

spin=93.75 # 92
th=170 # 157 # 196 # 137
f=230 # 102 # 230
N=151 # 199
iter_start=2140 # 2200
iter_end=$iter_start
# case_start=4235484 # 418
case_start=4442140 # 45502140 
case_end=$case_start

FILENAME_BASE=shotimag$spin\th$th\f$f\fn

for snapshot in `seq $iter_start 1 $iter_end`; do
  for case in `seq $case_start 1 $case_end`; do

    upper_left_panel=$FILENAME_BASE$snapshot\case$case\_$N\_IQUV_xy.png
    echo $upper_left_panel
    upper_right_panel=$FILENAME_BASE$snapshot\case$case\_$N\_I-LP-EVPA-CP_xy.png
    lower_left_panel=$FILENAME_BASE$snapshot\case$case\_$N\_IQUV_uv.png
    lower_right_panel=$FILENAME_BASE$snapshot\case$case\_$N\_miniversion_4panel_uv.png

    output_filename=16panel_$snapshot\_$case\.png

    convert \( $upper_left_panel $upper_right_panel -append \) \( $lower_left_panel $lower_right_panel -append \) +append $output_filename

######
done # 
done #
######

#############
## ANIMATE ##
#############

# credit: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
if hash mencoder 2>/dev/null; 
  then mencoder mf://16panel_????_???.png -mf fps=10:type=png -o movie.mp4 -ovc x264 -of lavf -lavfopts format=mp4
fi
#############
