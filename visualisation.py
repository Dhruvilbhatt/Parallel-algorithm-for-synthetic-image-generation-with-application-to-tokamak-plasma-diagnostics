import numpy as np
import matplotlib.pyplot as plt

intensity_data = [] # Array to store the intensity data of each pixel

dimnension = 0 # Dimension of the image

with open('pixels_intensity_txt.txt') as f: 
    for line in f:
        intensity_data.append(float(line))

    intensity_data = np.array(intensity_data)
    dimnension = int(np.sqrt(len(intensity_data)))
    intensity_data = intensity_data.reshape(dimnension, dimnension)

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 9})
plt.imshow(intensity_data, cmap = 'viridis', extent = [0,dimnension,0,dimnension]) # Visualising the plasma image
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(13)
    label.set_fontweight("bold") 
plt.xlabel('\n\n(a)', fontsize = 16, fontweight="bold")

plt.show()
