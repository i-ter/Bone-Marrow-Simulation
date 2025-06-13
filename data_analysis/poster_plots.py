import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


run_id = 26
df = pd.read_csv(f'./data/results/smc_seed_{run_id}_all_steps.csv')


def get_step_data(step: int) -> pd.DataFrame:
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
        size = 4
        alpha = 0.3

        if row['cell_type'] <= 9:
            color = 'red'
            alpha = 0.7
        elif row['cell_type'] == 40:
            size = 6
            color = 'gray'
        elif row['cell_type'] == 27:
            size = 8
            color = 'gray'
        circle = plt.Circle((row['x'], row['y']), size, color=color, alpha=alpha)
        ax.add_patch(circle)

    ax.set_xlim(0, 500)
    ax.set_ylim(0, 500)
    ax.set_aspect('equal')
    ax.set_xticks([])  # Remove x-axis ticks
    ax.set_yticks([])  # Remove y-axis ticks
    # plt.title(f'Cells at step={step_data.step.values[0]}')
    plt.savefig(f'./figures/mot_seed_{run_id}_frame_{step_data.step.values[0]}.png',
                 facecolor='black', bbox_inches='tight', pad_inches=0)
    # plt.show()


def plot_clonal_frequency(df: pd.DataFrame):
    clones_count = df.groupby(['step', 'clone_id']).size().reset_index(name='count')
    clones_count = clones_count[clones_count['clone_id'] >= 0]
    clones_count['percentage'] = clones_count.groupby('step')['count'].transform(lambda x: x / x.sum())

    plt.figure(figsize=(8, 5), dpi=200, facecolor='white')
    
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
    plt.text(140000, threshold, 'Dominance Threshold', fontsize=10, 
             ha='center', va='bottom', color='black', fontweight='bold')
    
    # plt.title('Clonal Dynamics in a Simulation', fontweight='bold', fontsize=15)
    plt.xlabel('Time (simulation steps)', fontsize=17)
    plt.ylabel('Clone size (%)', fontsize=17)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=13)
    
    # Update legend title and position it at the bottom
    legend = plt.legend(title='Clone ID', fontsize=15, loc='upper right')
    legend.get_title().set_fontweight('bold')
    plt.tight_layout()
    plt.savefig(f'./figures/smc_seed_{run_id}_clonal_frequency.png', facecolor='white', pad_inches=0.3)

    # plt.show()

# plot_cells(get_step_data(50000))
plot_clonal_frequency(df)