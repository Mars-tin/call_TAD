# README

1. Download HiC Data from: https://www.dropbox.com/s/3oooszaorhomo1z/HSPC_downsample.hic?dl=0, and put the .bed file in `data/`.
2. Run `python3 hic2matrix.py`. The output files are in `data/mat/`.
3. Start `R`. Select `TopDom.R` for this project.
> source("utils/TopDom.R")
> TopDom(matrix.file="data/mat/<filename>.matrix", window.size=10, outFile=<filename>)
4. There will be 3 output file. Usually the .bed file would be enough for our project.
