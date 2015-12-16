import vlbi_imaging_utils_v3 as vb

# Load the image and array
im = vb.load_fits('roman.fits') #for a fits image like the one attached
#im = vb.load_textim('roman.txt') #for a text file like the one attached

# DEFAULT
# eht = vb.load_array('EHT2017.txt') #see the attached array text file
eht = vb.load_array('EHT-SMT.txt') #see the attached array text file

# Look at the image
im.display()

# Observe the image
# tint_sec is the integration time in seconds, and tadv_sec is the advance time between scans
# tstart_hr is the GMST time of the start of the observation and tstop_hr is the GMST time of the end
# bw_hz is the  bandwidth in Hz
# sgrscat=True blurs the visibilities with the hardcoded Sgr A* scattering kernel for the appropriate image frequency
tint_sec = 60
tadv_sec = 600
tstart_hr = 0
tstop_hr = 24
bw_hz = 4e9
obs = im.observe(eht, tint_sec, tadv_sec, tstart_hr, tstop_hr, bw_hz, sgrscat=True)

# There are some simple plots you can check
obs.plotall('u','v') # uv coverage
obs.plotall('uvdist','amp') # amplitude with baseline distance'
# obs.plot_bl('SMA','ALMA','phase') # visibility phase on a baseline over time
# obs.plot_cphase('SMA', 'SMT', 'ALMA') # closure phase on a triangle over time

tight_layout() # ROMAN: AVOID OVERLAPPING LABELS, ETC

# You can deblur the visibilities by dividing by the scattering kernel
obs = vb.deblur(obs)

# Export the visibility data to uvfits/text
obs.export_uvfits('obs.uvfits') # this requires the template.UVP file to be in the directory you're working in

#DEFAULT
obs.export_txt('obs.txt') # exports a text file with the visibilities
# obs.export_txt('obs-SMT-SMA.txt') # exports a text file with the visibilities
