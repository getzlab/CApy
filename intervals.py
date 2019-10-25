import pandas as pd
import numpy as np

# Returns dataframe with chromosome lengths
def get_chrlen(assembly='hg19',primary=True):
    if (assembly=='hg19')|(assembly=="GRCh37"):
        ref_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"
    elif (assembly=='hg38')|(assembly=="GRCh38"):
        ref_url = "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
    else:
        # Also can provide a url directly
        ref_url = assembly

    try:
        C = pd.read_csv(ref_url,sep="\t",header=None)
        C.columns=["chr","len"]
    except:
       print("Error: reference url invalid")

    # Only return primary nuclear contigs
    if primary:
        C = C[C['chr'].str.match('chr(\d+|X|Y)$')]
        C['num'] = C['chr'].str.replace('chr','').str.replace('X','23').str.replace('Y','24').astype(int)
        C = C.sort_values('num').reset_index(drop=True)

    return(C)


# Creates a dataframe of intervals binning genome
def make_intervals(window,assembly='hg19'):

    C = get_chrlen()

    Ilist = list()
    for ind,row in C.iterrows():
        st = np.arange(1,row['len']-window,window).astype(int)
        en = np.arange(window,row['len'],window).astype(int)
        Ilist.append(pd.DataFrame({"chr" : row['chr'],"start" : st,"end" : en}))
    I = pd.concat(Ilist)

    return(I)


