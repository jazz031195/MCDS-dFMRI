import numpy as np
import matplotlib.pyplot as plt
dist = np.random.gamma(3.5, 0.1, 25000)
plt.hist(dist)
plt.show()
for i in range(50):
    file = open("gamma"+str(i)+".txt","w")
    for d in dist[i:i+500]:
        if (d < 0.05) :
            d = 0.05
        file.write(str(d) + "\n")

    file.close()
