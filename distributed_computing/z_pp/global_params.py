#Commented parameters are set in naturerc.

#### For Figs. 1-3.

#dpi = 150 # 300 for final
#format = '.jpg' # TIFF CMYK for final

#width = 89. / 25.4 #3.503937
#hieght = 247. * 4./5. / 25.4 # 7.779527559 # In inches.  Leaving 1/5 of the printable page for captions.

file_name_minus_num = 'Noble_fig'

label_fontsize = 8
max_fontsize = 7
min_fontsize = 5

leg_prop_size = 8
capsize = 2
elinewidth = 2.0

# Flat design colors #1.
# https://kuler.adobe.com/Flat-design-colors-%231-color-theme-3044245/edit/?copy=true&base=1&rule=Custom&selected=1&name=Copy%20of%20Flat%20design%20colors%20%231&mode=hsv&rgbvalues=0.2,0.30196078431372547,0.3607843137254902,0.27058823529411763,0.6980392156862745,0.615686274509804,0.9372549019607843,0.788235294117647,0.2980392156862745,0.8862745098039215,0.47843137254901963,0.24705882352941178,0.8745098039215686,0.35294117647058826,0.28627450980392155&swatchOrder=0,1,2,3,4
colors = \
    {128: '#334D5C', 64: '#45B29D', 32: '#EFC94C', 16: '#E27A3F', 8: '#AB5350'}
# An not-quite-black.
kp = '#262626'

x_min = 0.50
x_max = 1.50
Phi_min = -0.55
Phi_max = 0.55

#stagger = 0.0025
#stagger_2 = 0.001


#### For EDF.  The maximum page dimension is 183mm by 240mm.

#dpi_edf = 150 # 300 for final with total file size < 10MB
#format_edf = '.jpg' # JPG RGB for final

max_width_edf = 183 / 25.4 # In inches.
max_hight_edf = 247 / 25.4 # In inches.

file_name_minus_num_edf = 'Noble_EDfig'

# For animations.
'''
Note:  If imshow plots a 2D array called X,
then X[0,0] is the value of the top-left
square; X[-1,-1] the values of the bottom-right.
(See http://dept.astro.lsa.umich.edu/~msshin/science/code
/matplotlib_cm/ to select a color map.)
'''
interval_in_ms = 1  # Changing this parameter does nothing in my experience.
blit_q = True       # True only updates parts altered parts of image.
fps = 2             # Frames per second.
num_frames = 20
num_frames_2 = 50
num_frames_3 = 20
anim_vmin = x_min
anim_vmax = x_max


# For reading and writing files.
pickle_protocol = 2
