"""Native Python vertical extrapolation for ISMIP remapped fields."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import xarray as xr

__all__ = [
    'extrapolate_vertical',
    'vertical_extrapolation',
]


def extrapolate_vertical(
    *,
    in_path: str,
    out_path: str,
    variable: str,
    z_name: str = 'z_extrap',
    logger: logging.Logger | None = None,
) -> None:
    """Vertically extrapolate a remapped field and write a NetCDF output."""
    from i7aof.io import read_dataset, write_netcdf

    log = logger or logging.getLogger(__name__)
    log.info(f'Python vertical extrapolation: {in_path} -> {out_path}')

    with read_dataset(in_path) as ds_in:
        if variable not in ds_in:
            raise KeyError(
                f"Variable '{variable}' not found in {in_path}. "
                f'Available: {", ".join(sorted(ds_in.data_vars))}'
            )
        if z_name not in ds_in:
            raise KeyError(
                f"Vertical coordinate '{z_name}' not found in {in_path}."
            )

        ds_out = ds_in.load()
        da = ds_out[variable]
        miss = _get_missing_value(da)

        transposed = da.transpose(z_name, 'y', 'x', 'time')
        data = np.asarray(transposed.values, dtype=np.float32)
        data = np.where(np.isnan(data), np.float32(miss), data)

        result = vertical_extrapolation(data=data, miss=np.float32(miss))
        ds_out[variable] = _restore_data_array(
            template=da,
            values=result,
            dims=(z_name, 'y', 'x', 'time'),
        )

    write_netcdf(
        ds_out,
        out_path,
        has_fill_values=[variable],
        format='NETCDF3_64BIT',
        progress_bar=False,
    )


def vertical_extrapolation(
    *, data: np.ndarray, miss: np.float32
) -> np.ndarray:
    """Apply the legacy vertical fill in native Python.

    Parameters
    ----------
    data : numpy.ndarray
        Array with shape ``(nz, ny, nx, nt)``.
    miss : numpy.float32
        Missing value sentinel.
    """
    result = np.array(data, copy=True, dtype=np.float32)
    if result.shape[0] < 2:
        return result

    surface_mask = (result[0] == miss) & (result[1] != miss)
    result[0][surface_mask] = result[1][surface_mask]

    for z_index in range(1, result.shape[0]):
        fill_mask = (result[z_index] == miss) & (result[z_index - 1] != miss)
        result[z_index][fill_mask] = result[z_index - 1][fill_mask]

    return result


def _get_missing_value(da: xr.DataArray) -> np.float32:
    for mapping in (da.encoding, da.attrs):
        for name in ('_FillValue', 'missing_value'):
            if name in mapping:
                value = np.asarray(mapping[name]).reshape(-1)[0]
                return np.float32(value)
    return np.float32(-1.0e6)


def _restore_data_array(
    *, template: xr.DataArray, values: np.ndarray, dims: tuple[str, ...]
) -> xr.DataArray:
    import xarray as xr

    da = xr.DataArray(
        values,
        dims=dims,
        coords={dim: template.coords[dim] for dim in dims},
    ).transpose(*template.dims)
    da.attrs = dict(template.attrs)
    da.encoding = dict(template.encoding)
    return da.astype(template.dtype, copy=False)
