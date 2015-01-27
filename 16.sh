## REARRANGE PLOTS IN 16 PANELS: "convert +/-append" ##
## AND ANIMATE: mencoder ##

# CREDIT: user fmw42 on
# http://www.imagemagick.org/discourse-server/viewtopic.php?t=14743
# convert \( image1 image2 -append \) \
# \( image3 image4 -append \) +append result

for snapshot in `seq 6100 1 6100`; do
  for case in `seq 100 1 199`; do

    upper_left_panel=shotimag93.75th140f230fn$snapshot\_299_$case\_IQUV_xy.png
    upper_right_panel=shotimag93.75th140f230fn$snapshot\_299_$case\_I-LP-EVPA-CP_xy.png
    lower_left_panel=shotimag93.75th140f230fn$snapshot\_299_$case\_IQUV_uv.png
    lower_right_panel=shotimag93.75th140f230fn$snapshot\_299_$case\_4panel_uv.png
    output_filename=16panel_$snapshot\_$case\.png

    convert \( $upper_left_panel $upper_right_panel -append \) \( $lower_left_panel $lower_right_panel -append \) +append $output_filename

    #convert \( shotimag93.75th140f230fn$snapshot\_299_IQUV_xy.png shotimag93.75th140f230fn$snapshot\_299_I-LP-EVPA-CP_xy.png -append \) \( shotimag93.75th140f230fn$snapshot\_299_IQUV_uv.png shotimag93.75th140f230fn$snapshot\_299_4panel_uv.png -append \) +append 16panel_$snapshot\.png
######
done # 
done #
######

#############
## ANIMATE ##
mencoder mf://16panel_????_???.png -mf fps=10:type=png -o movie.mp4 -ovc x264 -of lavf -lavfopts format=mp4
#############
