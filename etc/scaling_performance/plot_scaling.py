import re
import numpy as np
import matplotlib.pyplot as plt

# My log list
log_files = ['log_1.lammps', 'log_2.lammps', 'log_3.lammps', 'log_4.lammps']

def parse_log_file(file_content):
    # Extract data from log files
    performance_pattern = r'Performance:\s+(\d+\.\d+)\s+ns/day,'
    cpu_use_pattern = r'(\d+\.\d+)%\s+CPU use'
    
    performances = re.findall(performance_pattern, file_content)
    cpu_uses = re.findall(cpu_use_pattern, file_content)
    
    performances = [float(p) for p in performances]
    cpu_uses = [float(c) for c in cpu_uses]
    
    return performances, cpu_uses

def analyze_logs(log_files):
    sizes = []
    performance_data = []
    cpu_use_data = []
    for log_file in log_files:
        with open(log_file, 'r') as file:
            content = file.read()
            performances, cpu_uses = parse_log_file(content)
            if performances and cpu_uses:
                size = int(re.search(r'\d+', log_file).group(0))
                sizes.append(size)
                performance_data.append(performances)
                cpu_use_data.append(cpu_uses)
    return sizes, performance_data, cpu_use_data

def plot_performance(sizes, performance_data, cpu_use_data):
    performance_means = [np.mean(data) for data in performance_data]
    performance_stds = [np.std(data) for data in performance_data]
    cpu_use_means = [np.mean(data) for data in cpu_use_data]
    cpu_use_stds = [np.std(data) for data in cpu_use_data]
    
    plt.figure(figsize=(12, 6))
    sizes = np.array(sizes)
    performance_means = np.array(performance_means)
    
    plt.subplot(1,2,1)
    plt.plot(sizes**3, 1/performance_means, c='C0', marker='o', ls='-')
    plt.xlabel('System Size')
    plt.ylabel('Computing time (day/ns)')
    plt.title('Time vs System Size')
    
    plt.subplot(1,2,2)
    plt.errorbar(sizes, performance_means, yerr=performance_stds, c='C0', fmt='o', capsize=5)
    plt.ylabel('Performance (day/ns)')
    plt.title('Performance vs System Size')

    plt.tight_layout()
    plt.show()

# Read log file
sizes, performance_data, cpu_use_data = analyze_logs(log_files)

# Plot
plot_performance(sizes, performance_data, cpu_use_data)
