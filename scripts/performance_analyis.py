import matplotlib.pyplot as plt
import os

# Data for weak scaling
weak_domain_size = [40, 80, 160, 240]
weak_processors = [1, 2, 4, 8]
weak_runtime = [2, 31, 380, 1590]

# Data for strong scaling
strong_processors = [1, 2, 4, 8, 12]
strong_runtime = [114, 75, 51, 50, 34]

# Calculate speedup and efficiency for weak scaling
weak_speedup = [weak_runtime[0] / rt for rt in weak_runtime]
weak_efficiency = [sp / proc for sp, proc in zip(weak_speedup, weak_processors)]

# Calculate speedup and efficiency for strong scaling
strong_speedup = [strong_runtime[0] / rt for rt in strong_runtime]
strong_efficiency = [sp / proc for sp, proc in zip(strong_speedup, strong_processors)]

img_dir = 'imgs'
# Plotting weak scaling results
# Plotting strong scaling results
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Plot strong scaling runtime
axs[0].plot(strong_processors, strong_runtime, marker='o', color='blue')
axs[0].set_title('Strong Scaling: Runtime vs. Processors')
axs[0].set_xlabel('Number of Processors')
axs[0].set_ylabel('Runtime (s)')
axs[0].grid(True)

# Plot strong scaling speedup
axs[1].plot(strong_processors, strong_speedup, marker='o', color='red')
axs[1].set_title('Strong Scaling: Speedup vs. Processors')
axs[1].set_xlabel('Number of Processors')
axs[1].set_ylabel('Speedup')
axs[1].grid(True)

# Plot strong scaling efficiency
axs[2].plot(strong_processors, strong_efficiency, marker='o', color='darkgreen')
axs[2].set_title('Strong Scaling: Efficiency vs. Processors')
axs[2].set_xlabel('Number of Processors')
axs[2].set_ylabel('Efficiency')
axs[2].grid(True)

# Adjust layout and save the plot for strong scaling
plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'strong_scaling.png'))
plt.show()


fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Plot weak scaling runtime
axs[0].plot(weak_processors, weak_runtime, marker='o', color='blue')
axs[0].set_title('Weak Scaling: Runtime vs. Processors')
axs[0].set_xlabel('Number of Processors')
axs[0].set_ylabel('Runtime (s)')
axs[0].grid(True)

# Plot weak scaling speedup
axs[1].plot(weak_processors, weak_speedup, marker='o', color='red')
axs[1].set_title('Weak Scaling: Speedup vs. Processors')
axs[1].set_xlabel('Number of Processors')
axs[1].set_ylabel('Speedup')
axs[1].grid(True)

# Plot weak scaling efficiency
axs[2].plot(weak_processors, weak_efficiency, marker='o', color='darkgreen')
axs[2].set_title('Weak Scaling: Efficiency vs. Processors')
axs[2].set_xlabel('Number of Processors')
axs[2].set_ylabel('Efficiency')
axs[2].grid(True)


# Adjust layout and save the plot for weak scaling
plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'weak_scaling.png'))
plt.show()

