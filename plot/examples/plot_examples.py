import matplotlib.pyplot as plt
import numpy as np

# Load your custom style (replace with your actual style file path if needed)
plt.style.use('../mystyle.mplstyle')

# Create a 2x2 subplot
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle('Matplotlib Style Showcase')

# Line plot
x = np.linspace(0, 10, 100)
y = np.sin(x)
axs[0, 0].plot(x, y, label='sin(x)')
axs[0, 0].set_title('Line Plot')
axs[0, 0].legend()

# Scatter plot
x = np.random.rand(50)
y = np.random.rand(50)
colors = np.random.rand(50)
sizes = 100 * np.random.rand(50)
axs[0, 1].scatter(x, y, c=colors, s=sizes, alpha=0.6, edgecolors='w')
axs[0, 1].set_title('Scatter Plot')

# Bar plot
categories = ['A', 'B', 'C', 'D']
values = [3, 7, 2, 5]
axs[1, 0].bar(categories, values, color='skyblue')
axs[1, 0].set_title('Bar Plot')

# Error bar plot
x = np.arange(5)
y = np.exp(-x/3)
error = 0.1 + 0.1 * x
axs[1, 1].errorbar(x, y, yerr=error, fmt='o-', capsize=5)
axs[1, 1].set_title('Error Bar Plot')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('test.png', dpi=300)
plt.show()

