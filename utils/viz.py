import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")


def set_display_params(config):
    print(config.context)
    if config.context == 'darktalk':
        sns.set(style=config.context_darktalk.style, context=config.context_darktalk.context,
                rc={"lines.linewidth": config.context_darktalk.llw, "legend.fontsize": config.context_darktalk.llw})
        plt.style.use(config.context_darktalk.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_darktalk.glw, "grid.alpha": config.context_darktalk.galpha,
                             'font.size': config.context_darktalk.fs})
        sns.set_palette(config.context_darktalk.palette)
    elif config.context == 'whitetalk':
        sns.set(style=config.context_whitetalk.style, context=config.context_whitetalk.context,
                rc={"lines.linewidth": config.context_whitetalk.llw, "legend.fontsize": config.context_whitetalk.llw})
        plt.style.use(config.context_whitetalk.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_whitetalk.glw, "grid.alpha": config.context_whitetalk.galpha,
                             'font.size': config.context_whitetalk.fs})
        sns.set_palette(config.context_paper.palette)
    elif config.context == 'paper':
        sns.set(style=config.context_paper.style, context=config.context_paper.context,
                rc={"lines.linewidth": config.context_paper.llw, "legend.fontsize": config.context_paper.llw})
        plt.style.use(config.context_paper.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_paper.glw, "grid.alpha": config.context_paper.galpha,
                             'font.size': config.context_paper.fs})
        sns.set_palette(config.context_paper.palette)
    else:
        raise ValueError(
            "unknown context {}. Should be either 'whitetalk', 'darktalk' or 'paper'".format(config.context))


if __name__ == "__main__":
    print('TODO')
