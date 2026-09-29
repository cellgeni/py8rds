import pandas as pd
import numpy as np
from numpy.dtypes import StringDType

from .robj import INT_NA
from .parser import parse_rds


def factor2array(robj):
    """
    Converts an Robj with factor to list

    Parameters
    ----------
    robj : Robj

    Returns
    -------
    list
    """
    val = None
    cl = robj.get("class")
    if (cl is not None) and "factor" in cl.value:
        val = np.empty(len(robj.value), dtype=StringDType(na_object=np.nan))
        levels = robj.get("levels").value
        for i in range(len(robj.value)):
            if robj.value[i] != INT_NA:
                val[i] = levels[robj.value[i] - 1]
            else:
                val[i] = np.nan
    return val


def as_data_frame(robj):
    """
    Converts an Robj to a pandas DataFrame.

    Parameters
    ----------
    robj : Robj | str
        Robj or path to an RDS file.

    Returns
    -------
    DataFrame

    Examples
    --------
    as_data_frame(robj)
    as_data_frame("data.rds")
    """
    if isinstance(robj, str):
        robj = parse_rds(robj)
    cols = {}
    names = robj.get("names").value
    for i in range(len(names)):
        val = robj.value[i]
        cl = val.get("class")
        if (cl is not None) and "factor" in cl.value:
            val = factor2array(val)
        else:
            val = val.value
        cols[names[i]] = val

    r = pd.DataFrame(cols)
    row_names = robj.get("row.names")
    if (row_names is not None) and (row_names.value is not None):
        if (len(row_names.value) != 2) or (
            row_names.value[0] != INT_NA
        ):  # to treat NULL rownames
            r.index = row_names.value

    # set int NAs, cannot do in before as numpy int32 doesn't support NAs
    int32_cols = r.select_dtypes(include=["int32"]).columns

    r[int32_cols] = r[int32_cols].astype("Int32").replace(INT_NA, pd.NA)

    return r


def as_numpy(robj):
    """
    Converts an Robj with an array/matrix into a dense or sparse NumPy array.

    Parameters
    ----------
    robj : Robj | str
        Robj or path to an RDS file.

    Returns
    -------
    np.array

    Examples
    --------
    as_numpy(robj)
    as_numpy("data.rds")
    """
    if isinstance(robj, str):
        robj = parse_rds(robj)
    cl = robj.getClass()
    if cl == "dgCMatrix":
        return _dgCMatrix2numpy(robj)
    if cl == "REALSXP":
        return _array2numpy(robj)


def as_anndata(robj):
    """
    Converts an Robj with an array/matrix into dense/sparse AnnData, keeping dimnames if any.

    Parameters
    ----------
    robj : Robj | str
        Robj or path to an RDS file.

    Returns
    -------
    ad.AnnData

    Examples
    --------
    as_anndata(robj)
    as_anndata("data.rds")
    """
    import anndata as ad  # imported lazily as it is slow to import

    if isinstance(robj, str):
        robj = parse_rds(robj)
    X = as_numpy(robj).T
    adata = ad.AnnData(X=X)
    dimnames = robj.get("dimnames")
    if dimnames is None:
        dimnames = robj.get("Dimnames")  # for sparse Matrices

    if dimnames is not None:
        dimnames = dimnames.value
        if dimnames[0].value is not None:
            adata.var_names = dimnames[0].value
        if dimnames[1].value is not None:
            adata.obs_names = dimnames[1].value
    return adata


def as_dict(robj):
    """
    Converts a top-level Robj to a dict.
    Expects that robj has a `names` slot of the same length as the number of values.

    Doesn't pay any attention to lower-level objects (so they can still be Robj).
    """
    names = robj.get("names").value
    r = {names[i]: robj.get(i).value for i in range(len(names))}
    return r


def seurat2adata(robj, assay=0, layer="counts"):
    """
    Converts a Seurat Robj into AnnData.
    It loads:
    1. specified assay/layer as data
    2. cell metadata
    3. var metadata if any
    4. all reduced dimensions assotiated with given assay

    Parameters
    ----------
    robj : Robj or str (path to an RDS file)
    assay : int or string - assay index or name
    layer : str - name of the layer to use

    Returns
    -------
    AnnData

    Examples
    --------
    seurat2adata(robj)
    """
    if isinstance(robj, str):
        robj = parse_rds(robj)

    assay_names = robj.get(["assays", "names"]).value.tolist()

    if isinstance(assay, str):
        assay = assay_names.index(assay)

    cnts = robj.get(["assays", assay, layer])
    obs = as_data_frame(robj.get("meta.data"))

    if cnts is not None:
        var = as_data_frame(robj.get(["assays", assay, "meta.features"]))
    # try Assay5
    else:
        names = robj.get(["assays", assay, "layers", "names"]).value
        layer_idx = np.where(names == layer)[0]
        if layer_idx.size == 0:
            raise ValueError(
                f"Layer '{layer}' not found in assay {assay}.\n"
                f"Following layers are available: {names}.\n"
                f"You can try something like py8rds.seurat2adata(srds,layer='{names[0]}')"
            )
        layer_idx = int(layer_idx[0])
        cnts = robj.get(["assays", assay, "layers", layer_idx])
        var = pd.DataFrame(
            index=robj.get(["assays", 0, "features", "dimnames", 0]).value
        )

    adata = as_anndata(cnts)

    # try to get obs names from assay if they absent in layer
    if is_default_index(adata.obs):
        obs_names = robj.get(["assays", assay, "cells", "dimnames", 0])
        if (obs_names is not None) and (len(obs_names.value) == adata.shape[0]):
            adata.obs_names = obs_names.value

    # try to keep dimnames if they are missed in obs/var
    if is_default_index(obs) and (obs.shape[0] == adata.shape[0]):
        obs.index = adata.obs_names

    # one more place to find obs_names
    if is_default_index(obs):
        obs_names = robj.get(["active.ident", "names"]).value
        if obs.shape[0] == len(obs_names):
            obs.index = obs_names

    if is_default_index(var):
        var.index = adata.var_names

    if not is_default_index(obs) and not is_default_index(adata.obs):
        obs = obs.loc[adata.obs_names, :]

    adata.obs = obs
    adata.var = var

    # load reduced dims
    rdims = robj.get(["reductions"])
    rdims_names = rdims.get("names")
    if rdims_names is not None:
        rdims_names = rdims_names.value
        for i in range(len(rdims_names)):
            assay_used = rdims.get([i, "assay.used", 0])
            if (assay_used is None) or assay_used == assay_names[assay]:
                adata.obsm[rdims_names[i]] = _array2numpy(
                    rdims.get([i, "cell.embeddings"])
                )

    return adata


def seurat2adata_spatial(robj, assay=0, layer="counts"):
    """
    Converts a Visium Seurat Robj into spatial AnnData.
    It loads:
    1. specified assay/layer as data
    2. cell metadata
    3. var metadata if any
    4. spatial coordinates into adata.obsm['spatial']
    5. spatial metadata into adata.uns['spatial']

    Parameters
    ----------
    robj : Robj or str (path to an RDS file)
    assay : int - assay index
    layer : str - name of the layer to use

    Returns
    -------
    AnnData
    """
    if isinstance(robj, str):
        robj = parse_rds(robj)
    adata = seurat2adata(robj, assay=assay, layer=layer)

    images = robj.get("images")
    if images is None:  # not spatial
        return adata

    img_names = images.get("names").value.tolist()

    spatial_uns = {}
    spatial_coords = None
    for i, name in enumerate(img_names):
        img_obj = images.get(i)
        img = as_numpy(img_obj.get("image"))
        coords = as_data_frame(img_obj.get("coordinates"))
        spatial_coords = pd.concat([spatial_coords, coords], axis=0)

        scalefactors = as_dict(img_obj.get("scale.factors"))
        scalefactors = {
            "fiducial_diameter_fullres": scalefactors["fiducial"][0],
            "spot_diameter_fullres": scalefactors["fiducial"][0]
            * 0.55,  # Seurat messed up spot sizes, but this looks about right
            "tissue_hires_scalef": scalefactors["hires"][0],
            "tissue_lowres_scalef": scalefactors["lowres"][0],
        }
        spatial_uns[name] = {
            "images": {"lowres": img},  # Seurat uses lowres by default
            "scalefactors": scalefactors,
        }

    if spatial_uns:
        adata.uns["spatial"] = spatial_uns
    adata.obsm["spatial"] = spatial_coords.loc[
        adata.obs_names, ["imagecol", "imagerow"]
    ].to_numpy()
    return adata


def _array2numpy(robj):
    """
    Converts an Robj with an array into a dense NumPy array.
    Expects an R array (dense).

    Parameters
    ----------
    robj : Robj

    Returns
    -------
    np.array

    Examples
    --------
    _array2numpy(robj)
    """
    dim = robj.get("dim")
    if dim is None:
        dim = len(robj.value)
    else:
        dim = dim.value
    X = np.array(robj.value).reshape(dim, order="F")
    return X


def _dgCMatrix2numpy(robj):
    from scipy import sparse  # imported lazily as it is slow to import

    i = robj.get("i").value
    p = robj.get("p").value
    dim = robj.get("Dim").value
    x = robj.get("x").value
    dimnames = robj.get("Dimnames").value
    X = sparse.csc_matrix((x, i, p), dim)
    return X


def is_default_index(df):
    idx = df.index

    if isinstance(idx, pd.RangeIndex):
        return idx.start == 0 and idx.step == 1 and idx.stop == len(df)

    try:
        return pd.Index(idx.astype(int)).equals(pd.RangeIndex(len(df)))
    except (ValueError, TypeError):
        return False
