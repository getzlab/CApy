import liftover
import tqdm
import pandas as pd
import numpy as np

def liftover_intervals(chrm, start, end, chainfile, from_name, to_name, funk_thresh = 60000):
    converter = liftover.ChainFile(chainfile, from_name, to_name)

    L = pd.DataFrame({ "chr" : chrm, "start" : start, "end" : end })
    L["chr_start_lift"] = 0
    L["chr_end_lift"] = 0
    L["start_lift"] = 0
    L["end_lift"] = 0
    L["start_strand_lift"] = 'x'
    L["end_strand_lift"] = 'x'
    for x in tqdm.tqdm(L.itertuples(), total = len(L)):
        conv = converter[x.chr][x.start]
        L.loc[x.Index, ["chr_start_lift", "start_lift", "start_strand_lift"]] = (-1, -1, '?') if len(conv) == 0 else conv[0]
        conv = converter[x.chr][x.end]
        L.loc[x.Index, ["chr_end_lift", "end_lift", "end_strand_lift"]] = (-1, -1, '?') if len(conv) == 0 else conv[0]

    # flip reverse strands
    idx = ((L["start_strand_lift"] == "-") & (L["end_strand_lift"] == "-")).values
    st_col = L.columns.get_loc("start_lift")
    en_col = L.columns.get_loc("end_lift")
    L.iloc[idx, [st_col, en_col]] = L.iloc[idx, [en_col, st_col]]
    L.loc[idx, ["start_strand_lift", "end_strand_lift"]] = "+"

    # remove funky intervals
    good_idx = (L["chr_start_lift"] == L["chr_end_lift"]) & (L["start_lift"] < L["end_lift"]) & (L["end_lift"] - L["start_lift"] < funk_thresh)

    return L.loc[good_idx]
