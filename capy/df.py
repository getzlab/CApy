import pandas as _pd

def multimap(df1, df2):
    """
    For each row in df1, finds the corresponding row in df2 (if it exists) and
    returns df2's index for that row.
    """
    assert "_multimap_" not in df2.columns
    orig_cols = list(df2.columns)
    return df1.merge(
      df2.rename_axis("_multimap_").reset_index(),
      left_on = list(df1.columns),
      right_on = orig_cols,
      how = "left"
    )["_multimap_"].astype(_pd.Int64Dtype())
