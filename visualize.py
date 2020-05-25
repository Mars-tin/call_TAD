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


def get_fold_enrichment(df, tp="domain", num_div=5):
    chromes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5',
               'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
               'chr21', 'chr22', 'chrX', 'chrY']
    df_tad = pd.read_csv("data/tad/hspc_10kb.csv")
    df_tad = df_tad.loc[df_tad['type'] == tp]
    length_tad = sum(df_tad['end'].values - df_tad['start'].values)

    # Calculate expected value
    length = sum(df['end'].values - df['start'].values)
    exp = length / get_size('all')
    # exp = length / length_tad

    # Calculate observed value
    fold_enrichment = None
    if num_div == 5:
        obs = np.array([0, 0, 0, 0, 0], dtype='f')
        num_chrs = 0
        for ch in chromes:
            region = int(get_size(ch)/5)
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
            j = 0
            div = 0
            obs_ch = np.array([0, 0, 0, 0, 0], dtype='f')
            for i, boundary in enumerate(list_tad):
                while list_signal[j][1] < boundary[0]:
                    j += 1
                    if j == len(list_signal):
                        j = -1
                        break
                if j < 0:
                    break
                if list_signal[j][0] > boundary[1]:
                    continue
                lower = max(list_signal[j][0], boundary[0])
                upper = min(list_signal[j][1], boundary[1])
                if region * (div + 1) > upper:
                    obs_ch[div] += upper - lower
                elif region * (div + 1) < lower:
                    div += 1
                    obs_ch[div] += upper - lower
                else:
                    obs_ch[div] += region * (div + 1) - lower
                    div += 1
                    obs_ch[div] += upper - region * (div + 1)

            obs_ch = obs_ch/region
            obs += obs_ch
        obs /= num_chrs
        fold_enrichment = obs / exp

    elif num_div == 1:
        obs = 0
        num_chrs = 0
        for ch in chromes:
            region = int(get_size(ch))
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
            j = 0
            obs_ch = 0
            for i, boundary in enumerate(list_tad):
                while list_signal[j][1] < boundary[0]:
                    j += 1
                    if j == len(list_signal):
                        j = -1
                        break
                if j < 0:
                    break
                if list_signal[j][0] > boundary[1]:
                    continue
                lower = max(list_signal[j][0], boundary[0])
                upper = min(list_signal[j][1], boundary[1])
                obs_ch += upper - lower

            obs_ch = obs_ch/region
            obs += obs_ch
        obs /= num_chrs
        fold_enrichment = obs / exp

    return fold_enrichment


if __name__ == '__main__':
    # Read the data-frame
    df_great_canyon = pd.read_excel("data/canyons/hcd34_canyons_sorted.xls", sheet_name="grant_canyons")
    df_short_canyon = pd.read_excel("data/canyons/hcd34_canyons_sorted.xls", sheet_name="short_canyons")
    df_ctcf = pd.read_csv("data/ctcf/ctcf.csv")

    # Obtain the fold enrichment
    div = 1
    for tp in ['domain', 'boundary', 'gap']:
        gc_fold_enrichment = get_fold_enrichment(df_great_canyon, tp, div)
        sc_fold_enrichment = get_fold_enrichment(df_short_canyon, tp, div)
        ctcf_fold_enrichment = get_fold_enrichment(df_ctcf, tp, div)

        # Make the bar chart
        ind = np.arange(div)
        width = 0.25
        plt.bar(ind, ctcf_fold_enrichment, width, label='CTCF')
        plt.bar(ind + width, gc_fold_enrichment, width, label='Grant Canyon')
        plt.bar(ind + 2*width, sc_fold_enrichment, width, label='Short Canyon')

        plt.ylabel('Fold Enrichment')
        plt.title('Fold Enrichment by Divided Normalized TAD {} Region'.format(tp.capitalize()))
        if div == 5:
            plt.xticks(ind + width / 2, ('Border 20%', 'Mid 20%', 'Central 20%', 'Mid 20%', 'Border 20%'))
        elif div == 1:
            plt.xticks(ind + width / 2, ('Total'))

        plt.legend(loc='upper left')
        plt.savefig("visualization_{}_{}.png".format(tp, div))
        plt.show()
