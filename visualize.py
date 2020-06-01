import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def get_size(ch):
    sizes = {
        'chr1': 248956422,
        'chr2': 242193529,
        'chr3': 198295559,
        'chr4': 190214555,
        'chr5': 181538259,
        'chr6': 170805979,
        'chr7': 159345973,
        'chr8': 145138636,
        'chr9': 138394717,
        'chr10': 133797422,
        'chr11': 135086622,
        'chr12': 133275309,
        'chr13': 114364328,
        'chr14': 107043718,
        'chr15': 101991189,
        'chr16': 90338345,
        'chr17': 83257441,
        'chr18': 80373285,
        'chr19': 58617616,
        'chr20': 64444167,
        'chr21': 46709983,
        'chr22': 50818468,
        'chrX': 156040895,
        'chrY': 57227415
    }
    if ch == 'all':
        return sum(sizes.values())
    return sizes[ch]


def has_intersection(domain1, domain2):
    return domain1[0] <= domain2[1] and domain1[1] >= domain2[0]


def get_intersection(domain1, domain2):
    if not has_intersection(domain1, domain2):
        return 0
    lower = max(domain1[0], domain2[0])
    upper = min(domain1[1], domain2[1])
    return upper-lower


def get_fold_enrichment(df, tp="domain", num_div=5, resol=50):
    chromes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
               'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
               'chr21', 'chr22', 'chrX', 'chrY']
    df_tad = pd.read_csv("data/tad/hspc_{}kb.csv".format(resol))
    df_tad = df_tad.loc[df_tad['type'] == tp]

    # Calculate expected value
    length = sum(df['end'].values - df['start'].values)
    exp = length / get_size('all')

    # Calculate observed value
    obs = np.array(np.zeros(num_div), dtype='f')
    num_chrs = 0

    for ch in chromes:
        df_ch = df.loc[df['chrom'] == ch]
        df_tad_ch = df_tad.loc[df_tad['chrom'] == ch]

        if len(df_ch) > 1:
            list_signal = np.asarray([df_ch['start'], df_ch['end']]).T
        else:
            list_signal = np.asarray([df_ch['start'].values, df_ch['end'].values]).T

        if len(list_signal) == 0:
            continue
        num_chrs += 1

        list_tad = np.asarray([df_tad_ch['start'], df_tad_ch['end']]).T
        obs_ch = np.array(np.zeros(num_div), dtype='f')

        for i, domain in enumerate(list_tad):

            region = (domain[1] - domain[0]) / num_div
            obs_reg = np.array(np.zeros(num_div), dtype='f')
            j = -1

            for idx, signal in enumerate(list_signal):
                if has_intersection(signal, domain):
                    j = idx
                    break

            while has_intersection(list_signal[j], domain):
                for k in range(num_div):
                    start = domain[0] + k * region
                    obs_reg[k] += get_intersection(list_signal[j], [start, start+region])
                j += 1
                if j == len(list_signal):
                    j = -1
                    break

            obs_reg = obs_reg / region
            obs_ch += obs_reg
            if j < 0:
                continue

        obs += obs_ch / len(list_tad)
    obs /= num_chrs
    fold_enrichment = obs / exp

    return fold_enrichment


def plot_enrichment(div=5, resol=25):
    # Read the data-frame
    df_great_canyon = pd.read_excel("data/canyons/hcd34_canyons_sorted.xls", sheet_name="grant_canyons")
    df_short_canyon = pd.read_excel("data/canyons/hcd34_canyons_sorted.xls", sheet_name="short_canyons")
    df_ctcf = pd.read_csv("data/ctcf/ctcf.csv")

    # Obtain the fold enrichment
    # ['domain', 'boundary', 'gap']
    for tp in ['domain']:
        gc_fold_enrichment = get_fold_enrichment(df_great_canyon, tp, div, resol)
        sc_fold_enrichment = get_fold_enrichment(df_short_canyon, tp, div, resol)
        ctcf_fold_enrichment = get_fold_enrichment(df_ctcf, tp, div, resol)

        # Make the bar chart
        ind = np.arange(div)
        width = 0.25
        # plt.bar(ind, ctcf_fold_enrichment, width, label='CTCF')
        # plt.bar(ind + 2 * width, sc_fold_enrichment, width, label='Short Canyon')
        plt.bar(ind + width, gc_fold_enrichment, width, label='Grant Canyon')

        plt.ylabel('Fold Enrichment')
        plt.title('Fold Enrichment by Divided Normalized TAD {} Region'.format(tp.capitalize()))

        plt.legend(loc='upper left')
        plt.savefig("visualization_{}_{}.png".format(resol, div))
        plt.show()


def plot_length(resol=25, bins=50):
    # Read the data-frame
    df_tad = pd.read_csv("data/tad/hspc_{}kb.csv".format(resol))
    df_tad = df_tad.loc[df_tad['type'] == "domain"]
    length_tad = df_tad['end'].values - df_tad['start'].values

    # Obtain the histgram
    plt.hist(length_tad, bins=bins)
    plt.axvline(length_tad.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(length_tad.mean() * 1.1, max_ylim * 0.9, 'Mean: {:.2f}'.format(length_tad.mean()))

    plt.ylabel('Count')
    plt.title('Count of lengths of TAD')
    plt.savefig("length_distribution_{}kb.png".format(resol))
    plt.show()


if __name__ == '__main__':
    plot_enrichment(100, 25)
