import matplotlib.pyplot as plt
import numpy as np

def get_color_dict():
    cmap = plt.get_cmap('Paired')
    
    colors = [cmap(index) for index in np.linspace(0, 1, num=12)]
    
    return {
        'LGWA': colors[1],
        'ET-HF': colors[2],
        'ET-LF': colors[3],
        'LISA': colors[7],
        'all': 'black'
    }