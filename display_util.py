#!Python
# ===================================================================
#ScriptName     Image Display utilities
# Purpose:      Use matplotlib for diagonosis 
# Author:       Anchi Cheng
# ===================================================================
import numpy as np

image_buffer = []

def addImage(arr, peak=None):
    """
    Add image to image_buffer for display
    # peak is correlation shift with 0,0 unshifted.
    # when peak is specified. arr is wrapped correlation image with unshifted
    # at the center of the correlation image
    """
    global image_buffer
    print(arr.shape)
    image_buffer.append(arr.copy())

    print(len(image_buffer))
    print(peak)
    print(peak is not None)
    if peak is not None:
        arr_min_shape = min(arr.shape)
        a = 0.01 * arr_min_shape
        b = 0.04 * arr_min_shape
        c = np.array(arr.shape)//2
        fill = (arr.max()-arr.min())*2 + arr.max()
        image_buffer[-1][int(c[0]-a):int(c[0]+a),int(c[1]-a):int(c[1]+a)] = fill 
        image_buffer[-1][int(-peak[0]+c[0]-b):int(-peak[0]+c[0]+b),int(-peak[1]+c[1]-b):int(-peak[1]+c[1]+b)] = fill

def showImages():
    import matplotlib.pyplot as plt
    number_of_buffer_images = len(image_buffer)
    print(number_of_buffer_images)
    width = number_of_buffer_images * 5
    fig, ax = plt.subplots(1,number_of_buffer_images, figsize=(width,4))
    for i in range(number_of_buffer_images):
        ax[i].imshow(image_buffer[i])
    plt.tight_layout()
    plt.show()
