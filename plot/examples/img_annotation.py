import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import numpy as np
import os

if os.path.exists("mystyle.mplstyle"):
    plt.style.use("mystyle.mplstyle")

# Create figure and gridspec
fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.5, 1], height_ratios=[1, 1])

# Generate sample data for scatter plot
np.random.seed(42)
x = np.random.randn(100)
y = np.random.randn(100)
colors = np.random.rand(100)

# Right block: Scatter plot (spans both rows of the right column)
ax_scatter = fig.add_subplot(gs[:, 0])
scatter = ax_scatter.scatter(x, y, c=colors, alpha=0.7, s=50)
ax_scatter.set_xlabel('X values')
ax_scatter.set_ylabel('Y values')
ax_scatter.set_title('Main Scatter Plot')
ax_scatter.grid(True, alpha=0.3)

# Add colorbar
# plt.colorbar(scatter, ax=ax_scatter, shrink=0.8)

# Highlight two specific points for annotation
highlight_indices = [25, 75]  # indices of points to highlight
for i, idx in enumerate(highlight_indices):
    ax_scatter.scatter(x[idx], y[idx], c='red', s=200, marker='o',
                      facecolors='none', edgecolors='red', linewidth=3)
    ax_scatter.annotate(f'Point {i+1}',
                       xy=(x[idx], y[idx]),
                       xytext=(20, 20),
                       textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', color='black', connectionstyle='arc3,rad=0'))

# Left top block: First PNG image
ax_img1 = fig.add_subplot(gs[0, 1])
try:
    # Replace 'image1.png' with your actual image path
    img1 = mpimg.imread('image1.png')
    ax_img1.imshow(img1)
    ax_img1.set_title('Annotation for Point 1')
except FileNotFoundError:
    # Fallback: create a placeholder
    ax_img1.text(0.5, 0.5, 'PNG Image 1\n(Point 1 Detail)',
                ha='center', va='center', transform=ax_img1.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    ax_img1.set_title('Annotation for Point 1')
ax_img1.set_xticks([])
ax_img1.set_yticks([])

# Left bottom block: Second PNG image
ax_img2 = fig.add_subplot(gs[1, 1])
try:
    # Replace 'image2.png' with your actual image path
    img2 = mpimg.imread('image2.png')
    ax_img2.imshow(img2)
    ax_img2.set_title('Annotation for Point 2')
except FileNotFoundError:
    # Fallback: create a placeholder
    ax_img2.text(0.5, 0.5, 'PNG Image 2\n(Point 2 Detail)',
                ha='center', va='center', transform=ax_img2.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    ax_img2.set_title('Annotation for Point 2')
ax_img2.set_xticks([])
ax_img2.set_yticks([])

for ax in [ax_img1, ax_img2]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

# Adjust layout to prevent overlapping
plt.tight_layout()
plt.show()

# Optional: Save the figure
# plt.savefig('gridspec_figure.png', dpi=300, bbox_inches='tight')
