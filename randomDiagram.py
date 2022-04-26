from codecs import xmlcharrefreplace_errors
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(8000)/1000 - 4

plt.figure(figsize=(9, 5))


values = np.sin(x)*1
values = values - np.floor(values)
plt.subplot(221)
plt.scatter(x, values,s = 0.01)
plt.xlabel('y = fract(sin(x)*1.0);')


values = np.sin(x)*2
values = values - np.floor(values)
plt.subplot(222)
plt.scatter(x, values,s = 0.01)
plt.xlabel('y = fract(sin(x)*2.0);')

values = np.sin(x)*10
values = values - np.floor(values)
plt.subplot(223)
plt.scatter(x, values,s = 0.01)
plt.xlabel('y = fract(sin(x)*10.0);')

values = np.sin(x)*100
values = values - np.floor(values)
plt.subplot(224)
plt.scatter(x, values,s = 0.01)
plt.xlabel('y = fract(sin(x)*100.0);')

plt.tight_layout()
plt.show()