from gw_landscape.waterfall import plot_all
from settings import FIG_PATH

if __name__ == '__main__':
    plot_all(FIG_PATH / 'single_detector_waterfall.pdf', log=True)