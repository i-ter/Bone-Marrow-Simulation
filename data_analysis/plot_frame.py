import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import enum
import numpy as np
import matplotlib as mpl
from itertools import cycle


run_id = 24
df = pd.read_csv(f'./data/results/smc_seed_{run_id}_all_steps.csv')


def get_step_data(step: int):
    d = df[df.step == step].copy()
    return d

colors = ['orange', 'red', 'green', 'blue', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'] # ten colors


def plot_cells(step_data: pd.DataFrame):
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


# plot_cells(get_step_data(50000))

def plot_clonal_frequency(df):
    clones_count =df.groupby(['step', 'clone_id']).size().reset_index(name='count')
    clones_count = clones_count[clones_count['clone_id'] >= 0]
    clones_count['percentage'] = clones_count.groupby('step')['count'].transform(lambda x: x / x.sum())
    d=clones_count.pivot(index='step', columns='clone_id', values='percentage')

    # Calculate total number of cells per step
    total_cells = df.groupby('step').size()

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100, facecolor='white')
    d.plot(cmap='Set1', ax=ax)

    # # Add secondary y-axis for total cell count
    # ax2 = ax.twinx()
    # total_cells.plot(ax=ax2, color='black', linewidth=2, alpha=0.7, label='Total Cells')
    # ax2.set_ylabel('Total Number of Cells', fontsize=14, color='black')
    # ax2.tick_params(axis='y', labelsize=12, colors='black')
    # ax2.grid(False)

    # Add title and axis labels for poster printing
    ax.set_title('Fraction of Cells Belonging to Each Clone Over Time', fontsize=16, fontweight='bold')
    ax.set_xlabel('Simulation Step', fontsize=14)
    ax.set_ylabel('Fraction of Cells per Clone', fontsize=14)

    # Improve grid and tick parameters for clarity
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # # Optionally add a legend for the secondary axis
    # lines, labels = ax.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=12)

    plt.tight_layout()
    plt.show()

def plot_clonal_frequency_sns(df):
    clones_count = df.groupby(['step', 'clone_id']).size().reset_index(name='count')
    clones_count = clones_count[clones_count['clone_id'] >= 0]
    clones_count['percentage'] = clones_count.groupby('step')['count'].transform(lambda x: x / x.sum())

    plt.figure(figsize=(7, 5), dpi=200, facecolor='white')
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
    plt.title('Fraction of Cells Belonging to Each Clone Over Time', fontweight='bold', fontsize=15)
    plt.xlabel('Simulation Step', fontsize=12)
    plt.ylabel('% of cell population per clone', fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(f'./figures/smc_seed_{run_id}_clonal_frequency.png', facecolor='white', pad_inches=0.3)

    # plt.show()

# plot_clonal_frequency(df)
plot_clonal_frequency_sns(df)