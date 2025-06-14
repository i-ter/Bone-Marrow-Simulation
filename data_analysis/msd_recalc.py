import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

# Load the data
mpp_data = pd.read_excel('./data_analysis/MPP_CellTracks_RawValues.xlsx')
mpp_data['time'] = mpp_data['Time increment (3-min interval)']*3

print("Data shape:", mpp_data.shape)
print("\nFirst 10 rows:")
print(mpp_data.head(10))
print("\nColumn names:")
print(mpp_data.columns.tolist())
print("\nData info:")
print(mpp_data.info())

# Check if there are any identifying columns for different trajectories
print("\nLast 10 rows:")
print(mpp_data.tail(10))

# Look for patterns in time increments to identify trajectory boundaries
print("\nTime increment patterns:")
time_resets = mpp_data[mpp_data['Time increment (3-min interval)'] == 0].index.tolist()
print(f"Number of trajectory starts (time=0): {len(time_resets)}")
print(f"Trajectory start indices: {time_resets[:10]}")  # Show first 10

# Identify trajectory boundaries (where time increment resets to 0)
time_resets = mpp_data[mpp_data['Time increment (3-min interval)'] == 0].index.tolist()
time_resets.append(len(mpp_data))  # Add end index
print(f"Number of trajectories: {len(time_resets)-1}")

# Extract individual trajectories
trajectories = []
for i in range(len(time_resets)-1):
    start_idx = time_resets[i]
    end_idx = time_resets[i+1]
    traj = mpp_data.iloc[start_idx:end_idx].copy()
    trajectories.append(traj)

print(f"Trajectory lengths: min={min(len(t) for t in trajectories)}, max={max(len(t) for t in trajectories)}, median={np.median([len(t) for t in trajectories])}")

def calculate_msd(trajectory):
    """Calculate MSD for a single trajectory"""
    x = trajectory['x-coordinate (um)'].values
    y = trajectory['y-coordinate (um)'].values  
    z = trajectory['z-coordinate (um)'].values
    t = trajectory['time'].values
    
    msd_values = []
    lag_times = []
    
    # Calculate MSD for different lag times
    max_lag = len(trajectory) # Use up to half the trajectory length
    
    for lag in range(1, max_lag):
        squared_displacements = []
        
        for i in range(len(trajectory) - lag):
            dx = x[i + lag] - x[i]
            dy = y[i + lag] - y[i]
            dz = z[i + lag] - z[i]
            
            squared_displacement = dx**2 + dy**2 + dz**2
            squared_displacements.append(squared_displacement)
        
        if squared_displacements:
            msd = np.mean(squared_displacements)
            msd_values.append(msd)
            lag_times.append(t[lag] - t[0])
    
    return np.array(lag_times), np.array(msd_values)

# Calculate MSD for each trajectory
all_lag_times = []
all_msd_values = []

print("Calculating MSD for each trajectory...")
for i, traj in enumerate(trajectories):
    if len(traj) > 10:  # Only use trajectories with sufficient length
        lag_times, msd_values = calculate_msd(traj)
        if len(lag_times) > 0:
            all_lag_times.append(lag_times)
            all_msd_values.append(msd_values)

print(f"Successfully calculated MSD for {len(all_msd_values)} trajectories")

# Ensemble average MSD
# Find common time points and calculate ensemble average
max_common_time = min(max(lag_times) for lag_times in all_lag_times)
common_times = np.arange(3, int(max_common_time) + 1, 3)  # 3-minute intervals

ensemble_msd = []
msd_std = []

for t in common_times:
    msd_at_t = []
    for lag_times, msd_values in zip(all_lag_times, all_msd_values):
        # Find closest time point
        idx = np.argmin(np.abs(lag_times - t))
        if abs(lag_times[idx] - t) < 1.5:  # Within 1.5 minutes
            msd_at_t.append(msd_values[idx])
    
    if len(msd_at_t) > 0:
        ensemble_msd.append(np.mean(msd_at_t))
        msd_std.append(np.std(msd_at_t))

ensemble_msd = np.array(ensemble_msd)
msd_std = np.array(msd_std)

# Fit power law: MSD = D * t^alpha
def power_law(t, D, alpha):
    return D * np.power(t, alpha)

# Fit the ensemble MSD
valid_points = (ensemble_msd > 0) & (common_times > 0)
popt, pcov = curve_fit(power_law, common_times[valid_points], ensemble_msd[valid_points])
D_fit, alpha_fit = popt

print(f"\nMSD Scaling Analysis:")
print(f"Fitted power law: MSD = {D_fit:.4f} * t^{alpha_fit:.4f}")
print(f"Scaling exponent α = {alpha_fit:.4f} ± {np.sqrt(pcov[1,1]):.4f}")

# Interpret the scaling exponent
if abs(alpha_fit - 1.0) < 0.1:
    motion_type = "Normal diffusion"
elif alpha_fit < 1.0:
    motion_type = "Subdiffusive motion"
else:
    motion_type = "Superdiffusive motion"

print(f"Motion type: {motion_type}")

# Plotting
plt.figure(figsize=(12, 8))

# Plot 1: Individual trajectories (subset)
plt.subplot(2, 2, 1)
for i in range(min(10, len(all_lag_times))):
    plt.loglog(all_lag_times[i], all_msd_values[i], alpha=0.3, color='gray', linewidth=0.5)
plt.xlabel('Lag time (min)')
plt.ylabel('MSD (μm²)')
plt.title('Individual Trajectory MSDs')
plt.grid(True, alpha=0.3)

# Plot 2: Ensemble MSD with fit
plt.subplot(2, 2, 2)
plt.errorbar(common_times, ensemble_msd, yerr=msd_std, fmt='bo', alpha=0.7, label='Ensemble MSD')
plt.loglog(common_times, power_law(common_times, D_fit, alpha_fit), 'r-', 
           label=f'Fit: MSD ∝ t^{alpha_fit:.2f}', linewidth=2)
plt.xlabel('Lag time (min)')
plt.ylabel('MSD (μm²)')
plt.title('Ensemble MSD and Power Law Fit')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 3: Linear scale
plt.subplot(2, 2, 3)
plt.errorbar(common_times, ensemble_msd, yerr=msd_std, fmt='bo', alpha=0.7)
plt.plot(common_times, power_law(common_times, D_fit, alpha_fit), 'r-', linewidth=2)
plt.xlabel('Lag time (min)')
plt.ylabel('MSD (μm²)')
plt.title('Ensemble MSD (Linear Scale)')
plt.grid(True, alpha=0.3)

# Plot 4: Histogram of trajectory lengths
plt.subplot(2, 2, 4)
traj_lengths = [len(t) for t in trajectories]
plt.hist(traj_lengths, bins=30, alpha=0.7, edgecolor='black')
plt.xlabel('Trajectory Length (time points)')
plt.ylabel('Frequency')
plt.title('Distribution of Trajectory Lengths')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Save results
results_df = pd.DataFrame({
    'lag_time_min': common_times,
    'ensemble_msd_um2': ensemble_msd,
    'msd_std_um2': msd_std
})

results_df.to_csv('./data_analysis/mpp_msd_results.csv', index=False)

print(f"\nResults saved to:")
print(f"- Plot: './data_analysis/msd_analysis.png'")
print(f"- Data: './data_analysis/msd_results.csv'")
print(f"\nSummary statistics:")
print(f"- Number of trajectories analyzed: {len(all_msd_values)}")
print(f"- Maximum lag time: {max(common_times)} min")
print(f"- MSD at 30 min: {ensemble_msd[np.argmin(np.abs(common_times - 30))]:.2f} μm²")
print(f"- MSD at 60 min: {ensemble_msd[np.argmin(np.abs(common_times - 60))]:.2f} μm²" if max(common_times) >= 60 else "- MSD at 60 min: Not available")







