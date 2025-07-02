import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.transforms import ScaledTranslation
import numpy as np

def add_subplot_label(fig, ax, label, offset_x=40, offset_y=30, fontsize=14, fontweight='bold', color='black'):
    """
    Add an alphabetical label to the upper left corner outside of a subplot.
    Position is independent of figure size and aspect ratio.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to label
    label : str
        The label text (e.g., 'A', 'B', 'C', 'D')
    offset_x : float, default=10 pixels
        Horizontal offset as pixels from the left spine
        (positive = left of spine, negative = right of spine)
    offset_y : float, default=10 pixels
        Vertical offset as pixels above the top spine
        (positive = above spine, negative = below spine)
    fontsize : int, default=14
        Font size for the label
    fontweight : str, default='bold'
        Font weight for the label
    color : str, default='black'
        Color of the label text
    """
    # Use axes coordinates for size-independent positioning
    # Position relative to axes: (0,0) = bottom-left, (1,1) = top-right
    # Negative values place labels outside the axes boundaries


    offset_x_inches, offset_y_inches = - offset_x / fig.dpi, offset_y / fig.dpi
    offset = ScaledTranslation(offset_x_inches, offset_y_inches, fig.dpi_scale_trans)

    # Add text using axes coordinates (transform=ax.transAxes)
    ax.text(0, 1, label,
            transform=ax.transAxes + offset,  # Use axes coordinates
            fontsize=fontsize, fontweight=fontweight,
            color=color, ha='center', va='center',
            bbox=dict(boxstyle='square,pad=0.3', facecolor='white',
                     edgecolor='none', linewidth=1),
            clip_on=False)  # Allow text to appear outside axes boundaries

# Example usage with a 2x2 grid
def create_labeled_gridspec_example():
    """Create an example with 2x2 gridspec layout with alphabetical labels"""

    # Create figure and gridspec
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.3)

    # Labels in order: left to right, top to bottom
    labels = ['A', 'B', 'C', 'D']

    # Create subplots and add labels
    axes = []

    # Top row
    ax1 = fig.add_subplot(gs[0, 0])  # Top-left -> A
    ax2 = fig.add_subplot(gs[0, 1])  # Top-right -> B

    # Bottom row
    ax3 = fig.add_subplot(gs[1, 0])  # Bottom-left -> C
    ax4 = fig.add_subplot(gs[1, 1])  # Bottom-right -> D

    axes = [ax1, ax2, ax3, ax4]

    # Add some sample content to each subplot
    np.random.seed(42)

    # Subplot A: Line plot
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    ax1.plot(x, y, 'b-', linewidth=2)
    ax1.set_title('Sine Wave')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.grid(True, alpha=0.3)

    # Subplot B: Scatter plot
    x_scatter = np.random.randn(50)
    y_scatter = np.random.randn(50)
    ax2.scatter(x_scatter, y_scatter, alpha=0.6, c='red')
    ax2.set_title('Random Scatter')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.grid(True, alpha=0.3)

    # Subplot C: Bar chart
    categories = ['Cat 1', 'Cat 2', 'Cat 3', 'Cat 4']
    values = [23, 45, 56, 78]
    ax3.bar(categories, values, color='green', alpha=0.7)
    ax3.set_title('Bar Chart')
    ax3.set_ylabel('Values')

    # Subplot D: Histogram
    data = np.random.normal(0, 1, 1000)
    ax4.hist(data, bins=30, color='purple', alpha=0.7)
    ax4.set_title('Histogram')
    ax4.set_xlabel('Value')
    ax4.set_ylabel('Frequency')

    # Add alphabetical labels to each subplot
    for ax, label in zip(axes, labels):
        add_subplot_label(fig, ax, label)

    plt.show()
    return fig, axes

# Example with custom gridspec layout (like your original scenario)
def create_custom_labeled_layout():
    """Create the custom layout from your first question with alphabetical labels"""

    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[1, 1])

    # Left top block -> A
    ax_A = fig.add_subplot(gs[0, 0])
    ax_A.text(0.5, 0.5, 'Image 1\nPlaceholder', ha='center', va='center',
              transform=ax_A.transAxes, bbox=dict(boxstyle='round', facecolor='lightblue'))
    ax_A.set_title('Top Left Block')
    ax_A.set_xticks([])
    ax_A.set_yticks([])

    # Right block (spans both rows) -> B
    ax_B = fig.add_subplot(gs[:, 1])
    np.random.seed(42)
    x = np.random.randn(100)
    y = np.random.randn(100)
    ax_B.scatter(x, y, alpha=0.6)
    ax_B.set_title('Main Scatter Plot')
    ax_B.set_xlabel('X values')
    ax_B.set_ylabel('Y values')
    ax_B.grid(True, alpha=0.3)

    # Left bottom block -> C
    ax_C = fig.add_subplot(gs[1, 0])
    ax_C.text(0.5, 0.5, 'Image 2\nPlaceholder', ha='center', va='center',
              transform=ax_C.transAxes, bbox=dict(boxstyle='round', facecolor='lightgreen'))
    ax_C.set_title('Bottom Left Block')
    ax_C.set_xticks([])
    ax_C.set_yticks([])

    # Add labels (left to right, top to bottom order)
    add_subplot_label(fig, ax_A, 'A')  # Top-left
    add_subplot_label(fig, ax_B, 'B')  # Right (spans both rows, so gets B)
    add_subplot_label(fig, ax_C, 'C')  # Bottom-left

    plt.tight_layout()
    plt.show()
    return fig, [ax_A, ax_B, ax_C]

# Example with different figure sizes to demonstrate size-independence
def test_size_independence():
    """Test the labeling function with different figure sizes"""

    sizes = [(8, 6), (12, 8), (16, 10), (6, 12)]  # Different aspect ratios and sizes

    for i, (width, height) in enumerate(sizes):
        fig = plt.figure(figsize=(width, height))
        gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.3)

        # Create 2x2 grid
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        ax4 = fig.add_subplot(gs[1, 1])

        axes = [ax1, ax2, ax3, ax4]
        labels = ['A', 'B', 'C', 'D']

        # Add simple content
        for j, ax in enumerate(axes):
            ax.plot([0, 1], [0, 1], 'o-')
            ax.set_title(f'Subplot {labels[j]}')

        # Add labels - position should be consistent across all figure sizes
        for ax, label in zip(axes, labels):
            add_subplot_label(fig, ax, label)

        fig.suptitle(f'Figure Size: {width}"×{height}"')
        plt.show()

# Run examples
if __name__ == "__main__":
    print("Creating 2x2 grid example...")
    fig1, axes1 = create_labeled_gridspec_example()

    print("Creating custom layout example...")
    fig2, axes2 = create_custom_labeled_layout()

    print("Testing size independence...")
    test_size_independence()
