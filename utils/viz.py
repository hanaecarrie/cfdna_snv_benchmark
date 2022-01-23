import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")


def set_display_params(config):
    if config.context == 'darktalk':
        sns.set(style=config.context_talk.style, context=config.context_talk.context,
                rc={"lines.linewidth": config.context_talk.llw, "legend.fontsize": config.context_talk.llw})
        plt.style.use(config.context_talk.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_talk.glw, "grid.alpha": config.context_talk.galpha,
                             'font.size': config.context_talk.fs})
        sns.set_palette(config.context_talk.palette)
    elif config.context == 'whitetalk':
        print(config.context_paper)
        sns.set(style=config.context_paper.style, context=config.context_paper.context,
                rc={"lines.linewidth": config.context_paper.llw, "legend.fontsize": config.context_paper.llw})
        plt.style.use(config.context_paper.styleuse)
        plt.rcParams.update({"grid.linewidth": config.context_paper.glw, "grid.alpha": config.context_paper.galpha,
                             'font.size': config.context_paper.fs})
        sns.set_palette(config.context_paper.palette)
    elif config.context == 'paper':
        print(config.context_paper)
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
