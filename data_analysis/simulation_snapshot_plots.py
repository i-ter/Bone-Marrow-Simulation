import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def get_step_data(df: pd.DataFrame, step: int) -> pd.DataFrame:
    d = df[df.step == step].copy()
    return d



def plot_cells(step_data: pd.DataFrame):
    colors = ['orange', 'red', 'green', 'blue', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

    plt.figure(figsize=(5, 5), dpi=500, facecolor='black')
    ax = plt.gca()
    ax.set_facecolor('black')
    
    for _, row in step_data.iterrows():
        # color = colors[row['clone_id'].astype(int)] if row['clone_id'] >= 0 else 'gray'
        color = 'blue'
        size = 2
        alpha = 0.6

        if row['cell_type'] <= 9:
            color = 'red'
            alpha = 0.9
            size = 3
        elif row['cell_type'] == 40:
            size = 6
            color = 'gray'
        elif row['cell_type'] == 27:
            size = 8
            color = 'gray'
        circle = plt.Circle((row['x'], row['y']), size, color=color, alpha=alpha)
        ax.add_patch(circle)
    # Load and plot blood vessels first (so they appear behind cells)
    vessels_df = pd.read_csv(f'./simulation/data/results_pres/mot_seed_{run_id}_vessels.csv')
    
    for _, vessel in vessels_df.iterrows():
        # Calculate the vessel as a rectangle
        start_x, start_y = vessel['start_x'], vessel['start_y']
        end_x, end_y = vessel['end_x'], vessel['end_y']
        radius = vessel['radius']
        
        # Calculate vessel length and angle
        length = ((end_x - start_x)**2 + (end_y - start_y)**2)**0.5
        
        rect = Rectangle(
            (min(start_x, end_x)-70, start_y-radius),
            abs(end_x - start_x)+140,
            2 * (radius),
            # facecolor='green',
            # alpha=0.5,  # Translucent
            edgecolor='green',
            linewidth=3,
            linestyle='--',
            fill=False,
            hatch='/////',
        )
        ax.add_patch(rect)

    ax.set_xlim(0, 500)
    ax.set_ylim(0, 500)
    ax.set_aspect('equal')
    ax.set_xticks([])  # Remove x-axis ticks
    ax.set_yticks([])  # Remove y-axis ticks
    # plt.title(f'Cells at step={step_data.step.values[0]}')
    plt.savefig(f'./figures/mot_seed_{run_id}_frame_{step_data.step.values[0]}.png',
                 facecolor='black', bbox_inches='tight', pad_inches=0)
    # plt.show()



def plot_cells_clonal(step_data: pd.DataFrame, smc: bool):

    plt.figure(figsize=(5, 5), dpi=500, facecolor='black')
    ax = plt.gca()
    ax.set_facecolor('black')

    clone_colors = {-1: 'gray', 0: 'yellow', 1: 'blue', 2: 'green'}
    
    for _, row in step_data.iterrows():
        # color = colors[row['clone_id'].astype(int)] if row['clone_id'] >= 0 else 'gray'
        color = clone_colors[row['clone_id']]
        size = 2
        alpha = 0.6

        if row['cell_type'] <= 9:
            color = 'red'
            alpha = 0.9
            size = 3
        elif row['cell_type'] == 40:
            size = 6
            color = 'gray'
        elif row['cell_type'] == 27:
            size = 8
            color = 'gray'
        circle = plt.Circle((row['x'], row['y']), size, color=color, alpha=alpha)
        ax.add_patch(circle)
    # Load and plot blood vessels first (so they appear behind cells)
    vessels_df = pd.read_csv(f'./simulation/data/results_pres/{'smc' if smc else 'mot'}_seed_{run_id}_vessels.csv')
    
    for _, vessel in vessels_df.iterrows():
        # Calculate the vessel as a rectangle
        start_x, start_y = vessel['start_x'], vessel['start_y']
        end_x, end_y = vessel['end_x'], vessel['end_y']
        radius = vessel['radius']
        
        # Calculate vessel length and angle
        length = ((end_x - start_x)**2 + (end_y - start_y)**2)**0.5
        
        rect = Rectangle(
            (min(start_x, end_x)-70, start_y-radius),
            abs(end_x - start_x)+140,
            2 * (radius),
            # facecolor='green',
            # alpha=0.5,  # Translucent
            edgecolor='green',
            linewidth=3,
            linestyle='--',
            fill=False,
            hatch='/////',
        )
        ax.add_patch(rect)

    ax.set_xlim(0, 500)
    ax.set_ylim(0, 500)
    ax.set_aspect('equal')
    ax.set_xticks([])  # Remove x-axis ticks
    ax.set_yticks([])  # Remove y-axis ticks
    # plt.title(f'Cells at step={step_data.step.values[0]}')
    plt.savefig(f'./simulation/figures/{'smc' if smc else 'mot'}_seed_{run_id}_frame_{step_data.step.values[0]}_clonal.png',
                 facecolor='black', bbox_inches='tight', pad_inches=0)
    # plt.show()

def plot_clonal_frequency(df: pd.DataFrame):
    clones_count = df.groupby(['step', 'clone_id']).size().reset_index(name='count')
    clones_count = clones_count[clones_count['clone_id'] >= 0]
    clones_count['percentage'] = clones_count.groupby('step')['count'].transform(lambda x: x / x.sum())

    plt.figure(figsize=(8, 5), facecolor='white')
    
    # Get the Set1 color palette
    set1_colors = sns.color_palette('Set1', n_colors=len(clones_count['clone_id'].unique()))
    clone_color_map = dict(zip(sorted(clones_count['clone_id'].unique()), set1_colors))
    
    # Plot the lines
    sns.lineplot(
        data=clones_count,
        x='step',
        y='percentage',
        hue='clone_id',
        palette='Set1',
        legend='full',
        linewidth=2,
        alpha=0.8
    )
    
    # Add colored bands above threshold for each clone
    threshold = 2/3
    for clone_id in clones_count['clone_id'].unique():
        clone_data = clones_count[clones_count['clone_id'] == clone_id].sort_values('step')
        above_threshold = (clone_data['percentage'] > threshold) & (clone_data['step'] >= 500)
        
        if above_threshold.any():
            plt.fill_between(
                clone_data['step'], 
                0, 
                1,
                where=above_threshold,
                color=clone_color_map[clone_id],
                alpha=0.3,
                interpolate=True,
                edgecolor=clone_color_map[clone_id],
                linewidth=2
            )

    plt.axhline(y=threshold, color='black', linestyle='--', linewidth=1.5, alpha=0.8)
    
    # Add text label on the threshold line
    # plt.text(140000, threshold, 'Dominance Threshold', fontsize=10, 
    #          ha='center', va='bottom', color='black', fontweight='bold')
    
    # plt.title('Clonal Dynamics in a Simulation', fontweight='bold', fontsize=15)
    plt.xlabel('Time (simulation steps)', fontsize=17)
    plt.ylabel('Clone size (%)', fontsize=17)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=13)
    
    # Update legend title and position it at the bottom
    legend = plt.legend(title='Clone ID', fontsize=15, loc='upper right')
    legend.get_title().set_fontweight('bold')
    plt.tight_layout()
    plt.savefig(f'./figures/mot_seed_{run_id}_clonal_frequency.png', pad_inches=0.3, dpi=200)

    # plt.show()

run_id = 11
smc = True  

df = pd.read_csv(f'./simulation/data/results_pres/{'smc' if smc else 'mot'}_seed_{run_id}_all_steps.csv')
df = df[df.step == 200000]
# plot_cells(get_step_data(50000))
plot_cells_clonal(df, smc)
# plot_clonal_frequency(df)