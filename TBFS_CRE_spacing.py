import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

# Load file
spacing_file_path = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/D1-D4_Rep1+2_TFBS/D1+D4_Rep1_2_TFBS_CRM_spacing_input.bed"
df = pd.read_csv(spacing_file_path, sep='\t', header=None)

# Name columns
df.columns = ['TFBS_chr', 'TFBS_start', 'TFBS_end', 'CRE_chr', 'CRE_start', 'CRE_end', 'CRE']

# Sort DataFrame by CRE and TFBS_start
df.sort_values(by=['CRE', 'TFBS_start'], inplace=True)


# Calculate spacing within each CRE
def calculate_spacing(group):
    group = group.sort_values('TFBS_start').reset_index(drop=True)
    spacing_list = [None]  # First TFBS has no preceding TFBS

    for i in range(1, len(group)):
        prev_end = group.loc[i - 1, 'TFBS_end']
        curr_start = group.loc[i, 'TFBS_start']

        # Calculate spacing
        if curr_start <= prev_end:
            spacing = -1  # Overlapping TFBS
        else:
            spacing = curr_start - prev_end  # Non-overlapping spacing

        spacing_list.append(spacing)

    group['TFBS_spacing'] = spacing_list
    return group


# Apply the function to each CRE group
df_with_spacing = df.groupby('CRE', group_keys=False).apply(calculate_spacing)

# Filter valid spacing values (ignore None)
valid_spacing = df_with_spacing['TFBS_spacing'].dropna()

# Convert to float for calculation
only_spacing = np.array([x for x in valid_spacing if isinstance(x, (int, float))])

# Calculate and print mean
spacing_mean = np.mean(only_spacing)
spacing_median = np.median(only_spacing)
spacing_range = max(only_spacing)
print(f"Mean TFBS Spacing: {spacing_mean}")
print(f"Median TFBS Spacing: {spacing_median}")
print(f"Max TFBS Spacing: {spacing_range}")
# Plot violin plot


only_spacing_modded: list = []
modded_value: int = 0
for i in only_spacing:
    modded_value = math.log(2+i) # plus 2 because we have -1 values and we can't do log on negative values
    only_spacing_modded.append(modded_value)
    modded_value = 0

plt.boxplot(only_spacing_modded)
plt.xlabel('TFBS in known CREs')
plt.ylabel('Log Spacing (bp)')
plt.show()

'''plt.hist(only_spacing,bins = 10,color = '#cc0000', edgecolor = 'white', log = True)
plt.title('TFBS Spacing')
plt.xlabel('TFBS spacing (bp)')
plt.axvline(np.median(only_spacing), color='k', linestyle='dashed', linewidth=1)
plt.ylabel('Frequency')
plt.show()'''
