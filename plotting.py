import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.patches as patches

def plot_flow_field(csv_path='logs/results.csv'):

    df = pd.read_csv(csv_path)

    x = df['x'].values
    y = df['y'].values
    p = df['p'].values
    u = df['u'].values
    v = df['v'].values

    grid_x, grid_y = np.mgrid[min(x):max(x):400j, min(y):max(y):100j]
    grid_p = griddata((x, y), p, (grid_x, grid_y), method='linear')
    grid_u = griddata((x, y), u, (grid_x, grid_y), method='linear')
    grid_v = griddata((x, y), v, (grid_x, grid_y), method='linear')

    velocity_magnitude = np.sqrt(grid_u**2 + grid_v**2)
    obstacle_mask = np.isclose(velocity_magnitude, 0.0, atol=1e-5)

    grid_p_masked = np.ma.masked_where(obstacle_mask, grid_p)
    grid_u_masked = np.ma.masked_where(obstacle_mask, grid_u)
    grid_v_masked = np.ma.masked_where(obstacle_mask, grid_v)

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(14, 16), constrained_layout=True)
    fig.suptitle('Flow Field Visualization (Re: 100, U_inlet: 1 m/s)', fontsize=20, weight='bold')


    def plot_field(ax, data, title, label, cmap):
        cf = ax.contourf(grid_x, grid_y, data, levels=100, cmap=cmap, extend='both')
        
        ax.contour(grid_x, grid_y, data, levels=15, colors='k', linewidths=0.5, alpha=0.5)
        ax.contourf(grid_x, grid_y, obstacle_mask, levels=[0.5, 1.5], colors='black')

        ax.set_title(title, fontsize=14, weight='bold')
        ax.set_ylabel('Y (m)', fontsize=12)
        ax.set_aspect('equal', adjustable='box')
        
        cbar = fig.colorbar(cf, ax=ax, aspect=30, pad=0.02)
        cbar.set_label(label, fontsize=12)
        return cf


    # Pressure Field
    plot_field(axes[0], grid_p_masked, 'Pressure Field Distribution', 'Pressure (Pa)', 'RdBu_r')
    # U-Velocity 
    plot_field(axes[1], grid_u_masked, 'Horizontal Velocity (U)', 'Velocity U (m/s)', 'viridis')
    # V-Velocity
    plot_field(axes[2], grid_v_masked, 'Vertical Velocity (V)', 'Velocity V (m/s)', 'viridis')

    axes[2].set_xlabel('X (m)', fontsize=12)

    output_filename = csv_path.replace('.csv', '_flow_fields.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to: {output_filename}")
    plt.close()

if __name__ == '__main__':
    print("Enter the path of simulation output CSV file:")
    csv_input = input().strip()
    plot_flow_field(csv_input)