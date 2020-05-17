import os
import numpy as np


def hic2txt(hic, txt, ch, resolution=10000):
    """
    Convert .hic file into txt (contact list)

    :param hic: str, path of input .hic file
    :param txt: str, path of output .txt file
    :param ch: str, chromosome
    :param resolution: int, default: 25000
    """
    juicer = 'juicer_tools_1.11.04_jcuda.0.8.jar'
    cmd = f'java -jar {juicer} dump observed NONE {hic} {ch} {ch} BP {resolution} {txt}'
    os.system(cmd)


def txt2matrix(txt, matrix, ch, length, resolution=10000):
    """
    Convert contact list into .matrix, size: N * (N + 3)

    :param txt: str, path of input .txt file
    :param matrix: str, path of output .matrix file
    :param ch: str, chromosome
    :param length: chromosome length
    :param resolution: int, default: 25000
    """
    f = open(txt)
    n_bins = length // resolution
    mat = np.zeros((n_bins, n_bins))

    for line in f:
        p1, p2, v = line.strip().split()
        p1, p2, v = int(p1), int(p2), float(v)

        if max(p1, p2) >= n_bins * resolution:
            continue

        mat[p1 // resolution, p2 // resolution] += v
        if p1 // resolution != p2 // resolution:
            mat[p2 // resolution, p1 // resolution] += v

    # np.savetxt('cd34_chr10_25kb_v0.matrix', mat, delimiter='\t', fmt='%i')

    f2 = open(matrix, 'w')
    for i in range(n_bins):
        f2.write(f'{ch}\t{i * resolution}\t{(i+1) * resolution}\t')
        lst = [str(int(elm)) for elm in mat[i]]
        f2.write('\t'.join(lst))
        f2.write('\n')
    f2.close()


sizes = {
    'chr1':	248956422,
    'chr2':	242193529,
    'chr3':	198295559,
    'chr4':	190214555,
    'chr5':	181538259,
    'chr6':	170805979,
    'chr7':	159345973,
    'chr8':	145138636,
    'chr9':	138394717,
    'chr10':	133797422,
    'chr11':	135086622,
    'chr12':	133275309,
    'chr13':	114364328,
    'chr14':	107043718,
    'chr15':	101991189,
    'chr16':	90338345,
    'chr17':	83257441,
    'chr18':	80373285,
    'chr19':	58617616,
    'chr20':	64444167,
    'chr21':	46709983,
    'chr22':	50818468,
    'chrX':	156040895,
    'chrY':	57227415
}

for ch, sz in sizes.items():
    hic2txt('../data/HSPC_downsample.hic', "../data/txt/hspc_{}_10kb.txt".format(ch), ch)
    txt2matrix('../data/txt/hspc_{}_10kb.txt'.format(ch), '../data/mat/hspc_{}_10kb.matrix'.format(ch), ch, sz, 10000)
