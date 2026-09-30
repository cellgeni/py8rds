import logging
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
    4. all reduced dimensions associated with given assay

    Parameters
    ----------
    robj : Robj or str (path to an RDS file)
    assay : int or string - assay index or name
    layer : str - name of the layer to use. For Assay5, if there is no such layer, split layers ({layer}.1, {layer}.2, ...)
        are concatenated, as Seurat::JoinLayers does

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
        adata = as_anndata(cnts)
        var = as_data_frame(robj.get(["assays", assay, "meta.features"]))
        if is_default_index(var):
            var.index = adata.var_names
    # try Assay5
    else:
        adata = _assay5_layer2adata(robj.get(["assays", assay]), layer)

        var = as_data_frame(robj.get(["assays", assay, "meta.data"]))
        var.index = robj.get(["assays", assay, "features", "dimnames", 0]).value
        if var.shape[0] != adata.shape[1]:
            logging.warning(
                f"Number of features in assay {assay} ({var.shape[0]}) does not match number of features in layer '{layer}' ({adata.shape[1]}). Subsetting features to match layer."
            )
        var = var.loc[adata.var_names, :]

    #try to keep dimnames if they are missed in obs. There are some weird cases...
    if is_default_index(obs) and (obs.shape[0] == adata.shape[0]):
        obs.index = adata.obs_names

    # one more place to find obs_names
    if is_default_index(obs):
        obs_names = robj.get(["active.ident", "names"]).value
        if obs.shape[0] == len(obs_names):
            obs.index = obs_names
    
    # make sure obs are the same
    if adata.shape[0] != obs.shape[0]:
        if is_default_index(obs) or is_default_index(adata.obs):
            raise ValueError("Incompatible shapes between AnnData and meta.data, and at least one of them has default index. Cannot subset.")
        logging.warning(
            f"Number of cells in layer '{layer}' ({adata.shape[0]}) does not match number of cells in meta.data ({obs.shape[0]}). Subsetting meta.data to match layer."
    )
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
                embedding = _array2numpy(rdims.get([i, "cell.embeddings"]))
                cell_names = rdims.get([i, "cell.embeddings", "dimnames", 0])
                if cell_names is not None:
                    embedding = (
                        pd.DataFrame(embedding, index=cell_names.value)
                        .loc[adata.obs_names]
                        .to_numpy()
                    )
                if embedding.shape[0] != adata.shape[0]:
                    raise ValueError(
                        f"Number of cells in reduced dimension '{rdims_names[i]}' ({embedding.shape[0]}) does not match number of cells in AnnData ({adata.shape[0]})."
                    )
                adata.obsm[rdims_names[i]] = embedding

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


def _logmap_names(logmap, name):
    """Returns row names of a SeuratObject LogMap that are TRUE in column `name`"""
    rownames, colnames = (d.value for d in logmap.get("dimnames").value)
    mask = np.asarray(logmap.value).reshape(logmap.get("dim").value, order="F")
    return rownames[mask[:, colnames.tolist().index(name)]]


def _assay5_layer2adata(assay, layer):
    """
    Converts a layer of a Seurat Assay5 into AnnData, taking cell and feature names from the assay.
    If there is no such layer but there are split layers ({layer}.1, {layer}.2, ...),
    they are concatenated (missing features are filled by zeros).
    """
    import anndata as ad  # imported lazily as it is slow to import

    names = assay.get(["layers", "names"]).value.tolist()
    if layer in names:
        layers = [layer]
    else:
        layers = [n for n in names if n.startswith(layer + ".")]
    if not layers:
        raise ValueError(
            f"Layer '{layer}' not found in assay.\n"
            f"Following layers are available: {names}.\n"
            f"You can try something like py8rds.seurat2adata(srds,layer='{names[0]}')"
        )

    layer_cells = [_logmap_names(assay.get("cells"), name) for name in layers]
    # split layers are made from one matrix, so each cell belongs to exactly one of them
    # calculate total number of cells in all layers and compare to number of unique cells
    n_cells = sum(len(c) for c in layer_cells)
    if len(layers) > 1 and len(np.unique(np.concatenate(layer_cells))) < n_cells:
        raise ValueError(
            f"Layer '{layer}' not found in assay. Supposedly split layers {layers} cannot be concatenated as they share cells.\n"
            f"Please specify one of them, for example py8rds.seurat2adata(srds,layer='{layers[0]}') or select different assay/layer"
        )

    adatas = []
    for name, cells in zip(layers, layer_cells):
        adata = as_anndata(assay.get(["layers", names.index(name)]))
        adata.obs_names = cells
        adata.var_names = _logmap_names(assay.get("features"), name)
        adatas.append(adata)
    if len(adatas) == 1:
        return adatas[0]

    logging.info(f"Concatenating split layers: {layers}")
    adata = ad.concat(adatas, join="outer", fill_value=0)
    return adata


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
