import numpy as np
import matplotlib.pyplot as plt

for i in range(0,2503):
    # Step 1: Read data from file
    first = 'spinImage'+str(i)

    txt = '.txt'

    data = np.loadtxt(first+txt)

    # Step 2: Prepare data for plotting

    # Step 3: Plotting
    plt.figure(figsize=(6, 4))  # Adjust the figure size as needed
    plt.imshow(data, cmap='hot', interpolation='nearest')
    plt.colorbar()  # Add color bar indicating the scale
    plt.title('Heatmap of Data')
    plt.xlabel('X-axis label')
    plt.ylabel('Y-axis label')
    plt.savefig('SpinImages/'+first+'.png')
    # plt.show()
    plt.close()
