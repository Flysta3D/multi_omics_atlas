import pandas as pd
import numpy as np

def bin1tobinx(bin1_data, binx=None, save=None) -> pd.DataFrame:
    """
    Converts bin size from 1 to other bin size.
    Args:
        bin1_data: bin1 lasso data.
        binx: The size of bin after conversion.
        save: If save is not None, save is the path to save the lasso data.
    Returns:
        Binx lasso data.
    """

    bin1_data["x"] = (bin1_data["x"] / binx).astype(int) * binx
    bin1_data["y"] = (bin1_data["y"] / binx).astype(int) * binx
    bin1_data = bin1_data.astype(str)
    bin1_data["MIDCounts"] = bin1_data["MIDCounts"].astype(int)
    bin1_data["groupby"] = bin1_data["geneID"] + "%" + bin1_data["x"] + "%" + bin1_data["y"]

    groupby_data = bin1_data.groupby(by="groupby")["MIDCounts"].sum().to_frame("MIDCounts").reset_index()
    binx_data = pd.DataFrame(list(groupby_data["groupby"].str.split("%")), columns=["geneID", "x", "y"])
    binx_data["geneID"] = binx_data["geneID"].astype("category")
    binx_data["x"] = binx_data["x"].astype(np.uint32)
    binx_data["y"] = binx_data["y"].astype(np.uint32)
    binx_data["MIDCounts"] = groupby_data["MIDCounts"].astype(np.uint16)

    if save is not None:
        binx_data.to_csv(save, index=False, sep="\t")
    return binx_data


def binxtobiny(binx_data, bin1_data, binx=None, biny=1, save=None) -> pd.DataFrame:
    """
    Converts bin size from x to 1.
    Args:
        binx_data: binx lasso data.
        bin1_data：bin1 lasso data.
        binx: The size of bin before conversion.
        biny: The size of bin after conversion.
        save: If save is not None, save is the path to save the lasso data.
    Returns:
        Biny lasso data.
    """

    binx_coords = binx_data.loc[:, ["x", "y"]].drop_duplicates()
    binx_coords.index, binx_coords.columns = range(len(binx_coords.index)), [
        "binx_x",
        "binx_y",
    ]

    bin1_data["binx_x"] = (bin1_data["x"] / binx).astype(int) * binx
    bin1_data["binx_y"] = (bin1_data["y"] / binx).astype(int) * binx
    bin1_need = pd.merge(bin1_data, binx_coords, on=["binx_x", "binx_y"], how="inner")
    del bin1_need["binx_x"], bin1_need["binx_y"]

    if biny != 1:
        biny_data = bin1tobinx(bin1_need, binx=biny, save=save)
        return biny_data
    else:
        if save is not None:
            bin1_need.to_csv(save, index=False, sep="\t")
        return bin1_need
