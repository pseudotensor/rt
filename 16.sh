## REARRANGE PLOTS IN 16 PANELS: "convert +/-append" ##
## AND ANIMATE: mencoder ##

# CREDIT: user fmw42 on
# http://www.imagemagick.org/discourse-server/viewtopic.php?t=14743
# convert \( image1 image2 -append \) \
# \( image3 image4 -append \) +append result

for snapshot in `seq 100 1 199`; do

    convert \( shotimag93.75th140f230fn6100_299_$snapshot\_IQUV_xy.png shotimag93.75th140f230fn6100_299_$snapshot\_I-LP-EVPA-CP_xy.png -append \) \( shotimag93.75th140f230fn6100_299_$snapshot\_IQUV_uv.png shotimag93.75th140f230fn6100_299_$snapshot\_4panel_uv.png -append \) +append 16panel_$snapshot\.png
######
done # 
######

#############
## ANIMATE ##
mencoder mf://16panel_???.png -mf fps=10:type=png -o thickdisk7-r-slices-fn6100.mp4 -ovc x264 -of lavf -lavfopts format=mp4
#############
