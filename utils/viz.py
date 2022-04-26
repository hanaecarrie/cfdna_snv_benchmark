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


def function_to_split(hand, labl, dividor):
    Hand_L = []
    Hand_M = []
    Labl_L = []
    Labl_M = []
    for h, l in zip(hand, labl):
        co = h.get_color()
        ls = h.get_linestyle()
        lw = h.get_linewidth()
        mk = h.get_marker()
        mew = h.get_markeredgewidth()
        ms = h.get_markersize()
        mkf = h.get_fillstyle()
        LABS = l.split(dividor)
        if len(LABS) != 2:
            print('Split Legends Error: Only exactly 1 Dividor is accepted.')
            print('                     Currently ' + str(len(LABS)-1) + ' dividors were given')
            return hand, labl
        # Line and Color
        LICO = plt.Line2D((0, 1), (0, 0), color=co, marker='', linestyle=ls,linewidth=lw)
        # Marker
        MARK = plt.Line2D((0, 1), (0, 0), color='k', marker=mk, markeredgewidth=mew, markersize=ms, linestyle='', fillstyle=mkf)
        #if mkf == 'full':
        #    MARK = plt.Line2D((0, 1), (0, 0), color='k', marker=mk, markeredgewidth=mew, markersize=ms, linestyle='')
        #else:
        #    MARK = plt.Line2D((0, 1), (0, 0), markeredgecolor='k', markerfacecolor='white', marker=mk, markeredgewidth=mew, markersize=ms, linestyle='', fillstyle=mkf)
        if LABS[0] not in Labl_L:
            Hand_L.append(LICO)
            Labl_L.append(LABS[0])
        if LABS[1] not in Labl_M:
            Hand_M.append(MARK)
            Labl_M.append(LABS[1])
    return Hand_L+Hand_M, Labl_L+Labl_M

