from PIL import Image
import os
import numpy as np
import cv2




data = np.zeros((400,640,4))
count = 0
for filename in os.listdir("./results"):
    img = Image.open("./results/"+filename)
    data += np.array(img)
    count+=1

data/=count

output = Image.fromarray(data.astype(np.uint8))
output.save("output.ppm")


