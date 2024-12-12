# Import necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set random seed for reproducibility
np.random.seed(725)

# make two fair dice (red and blue) ----
red_dice = np.random.choice(range(1, 7), size=50000, replace=True)
blue_dice = np.random.choice(range(1, 7), size=50000, replace=True)

# Combine into a DataFrame
dice_data = pd.DataFrame({'red': red_dice, 'blue': blue_dice})

# plot probability density function (uniform) ----
plt.hist(red_dice, bins=np.arange(1, 8) - 0.5, color='red', alpha=0.5, density=True)
plt.title('Probability Density of Red Dice')
plt.show()

plt.hist(blue_dice, bins=np.arange(1, 8) - 0.5, color='blue', alpha=0.5, density=True)
plt.title('Probability Density of Blue Dice')
plt.show()

# plot heatmap of outcomes ----
# Group by outcomes of red and blue dice and calculate probabilities
dice_grouped = dice_data.groupby(['red', 'blue']).size().reset_index(name='n')
dice_grouped['prob'] = 100 * dice_grouped['n'] / dice_grouped['n'].sum()

# Plot heatmap of outcomes
heatmap_data = dice_grouped.pivot(index='red', columns='blue', values='prob')


sns.heatmap(heatmap_data, annot=True, fmt=".1f", cmap="Reds", cbar_kws={'label': 'Probability (%)'})
plt.title('Heatmap of Red and Blue Dice Outcomes')
plt.show()

# make mirror dice ----
dice_data['mirror'] = dice_data['blue']

# Group by outcomes of mirror and blue dice and calculate probabilities
mirror_grouped = dice_data.groupby(['mirror', 'blue']).size().reset_index(name='n')
mirror_grouped['prob'] = 100 * mirror_grouped['n'] / mirror_grouped['n'].sum()

# Complete missing combinations of dice outcomes
mirror_complete = pd.merge(pd.MultiIndex.from_product([range(1, 7), range(1, 7)], names=['mirror', 'blue']).to_frame(index=False),
                           mirror_grouped, on=['mirror', 'blue'], how='left').fillna(0)

# Plot heatmap of mirror dice outcomes
mirror_heatmap_data = mirror_complete.pivot(index = 'mirror', columns = 'blue', values = 'prob')

sns.heatmap(mirror_heatmap_data, annot=True, fmt=".1f", cmap="Blues", cbar_kws={'label': 'Probability (%)'})
plt.title('Heatmap of Mirror and Blue Dice Outcomes')
plt.show()
