"""Native Python horizontal extrapolation for ISMIP remapped fields."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import xarray as xr

__all__ = [
    'extrapolate_horizontal',
    'horizontal_extrapolation',
]

_MIN_WEIGHT = np.float32(0.1)
_W1 = np.float32(0.946)
_W2 = np.float32(0.801)
_W11 = np.float32(0.895)
_W21 = np.float32(0.757)
_W22 = np.float32(0.641)
_W3 = np.float32(0.606)


def extrapolate_horizontal(
    *,
    in_path: str,
    out_path: str,
    basin_path: str,
    topo_path: str,
    variable: str,
    z_name: str = 'z_extrap',
    logger: logging.Logger | None = None,
    step3_max_iterations: int = 5,
) -> None:
    """Horizontally extrapolate a remapped field and write a NetCDF output."""
    import xarray as xr

    from i7aof.io import read_dataset, write_netcdf

    log = logger or logging.getLogger(__name__)
    log.info(f'Python horizontal extrapolation: {in_path} -> {out_path}')

    with (
        read_dataset(in_path) as ds_in,
        xr.open_dataset(basin_path) as ds_basin,
        xr.open_dataset(topo_path) as ds_topo,
    ):
        if variable not in ds_in:
            raise KeyError(
                f"Variable '{variable}' not found in {in_path}. "
                f'Available: {", ".join(sorted(ds_in.data_vars))}'
            )
        if z_name not in ds_in:
            raise KeyError(
                f"Vertical coordinate '{z_name}' not found in {in_path}."
            )
        if 'basinNumber' not in ds_basin:
            raise KeyError(
                f"'basinNumber' not found in basin file {basin_path}."
            )
        if 'bed' not in ds_topo:
            raise KeyError(f"'bed' not found in topography file {topo_path}.")
        ocean_name = _get_ocean_frac_name(ds_topo)

        ds_out = ds_in.load()
        da = ds_out[variable]
        miss = _get_missing_value(da)

        transposed = da.transpose(z_name, 'y', 'x', 'time')
        data = np.asarray(transposed.values, dtype=np.float32)
        data = np.where(np.isnan(data), np.float32(miss), data)

        z = np.asarray(ds_out[z_name].values, dtype=np.float32)
        basin_number = np.asarray(ds_basin['basinNumber'].values)
        bed = np.asarray(ds_topo['bed'].values, dtype=np.float32)
        ocean_frac = np.asarray(ds_topo[ocean_name].values, dtype=np.float32)

        result = horizontal_extrapolation(
            data=data,
            z=z,
            basin_number=basin_number,
            bed=bed,
            ocean_frac=ocean_frac,
            miss=np.float32(miss),
            logger=log,
            step3_max_iterations=step3_max_iterations,
        )

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


def horizontal_extrapolation(
    *,
    data: np.ndarray,
    z: np.ndarray,
    basin_number: np.ndarray,
    bed: np.ndarray,
    ocean_frac: np.ndarray,
    miss: np.float32,
    logger: logging.Logger | None = None,
    step3_max_iterations: int = 5,
) -> np.ndarray:
    """Apply the legacy 3-stage horizontal extrapolation in native Python.

    Parameters
    ----------
    data : numpy.ndarray
        Array with shape ``(nz, ny, nx, nt)``.
    z : numpy.ndarray
        Vertical coordinate values with shape ``(nz,)``.
    basin_number : numpy.ndarray
        IMBIE basin numbers with shape ``(ny, nx)``.
    bed : numpy.ndarray
        Bed topography with shape ``(ny, nx)``.
    ocean_frac : numpy.ndarray
        Ocean fraction with shape ``(ny, nx)``.
    miss : numpy.float32
        Missing value sentinel.
    logger : logging.Logger, optional
        Logger used for progress messages.
    step3_max_iterations : int, optional
        Maximum number of global fill passes away from the second level.
    """
    log = logger or logging.getLogger(__name__)
    result = np.array(data, copy=True, dtype=np.float32)
    basin_ids = range(int(np.max(basin_number)) + 1)
    stage3_limit = max(int(step3_max_iterations), 0)

    for z_index, z_level in enumerate(z):
        plane = result[z_index]
        below_level = bed < z_level

        for basin_id in basin_ids:
            basin_mask = basin_number == basin_id
            source_allowed = below_level & basin_mask
            target_allowed = source_allowed & (ocean_frac > np.float32(0.01))
            plane, iterations, filled = _iterate_stage(
                field=plane,
                source_allowed=source_allowed,
                target_allowed=target_allowed,
                miss=miss,
            )
            log.info(
                'Horizontal step1 z=%s basin=%s iterations=%s filled=%s',
                z_index,
                basin_id,
                iterations,
                filled,
            )

        for basin_id in basin_ids:
            basin_mask = basin_number == basin_id
            source_allowed = below_level & basin_mask
            plane, iterations, filled = _iterate_stage(
                field=plane,
                source_allowed=source_allowed,
                target_allowed=source_allowed,
                miss=miss,
            )
            log.info(
                'Horizontal step2 z=%s basin=%s iterations=%s filled=%s',
                z_index,
                basin_id,
                iterations,
                filled,
            )

        max_iterations = None if z_index == 1 else stage3_limit
        plane, iterations, filled = _iterate_stage(
            field=plane,
            source_allowed=np.ones(plane.shape[:2], dtype=bool),
            target_allowed=np.ones(plane.shape[:2], dtype=bool),
            miss=miss,
            max_iterations=max_iterations,
        )
        log.info(
            'Horizontal step3 z=%s iterations=%s filled=%s',
            z_index,
            iterations,
            filled,
        )
        result[z_index] = plane

    return result


def _iterate_stage(
    *,
    field: np.ndarray,
    source_allowed: np.ndarray,
    target_allowed: np.ndarray,
    miss: np.float32,
    max_iterations: int | None = None,
) -> tuple[np.ndarray, int, int]:
    iterations = 0
    total_filled = 0
    updated = np.array(field, copy=True, dtype=np.float32)

    while True:
        valid = updated[..., 0] != miss
        source_mask = source_allowed & valid
        target_mask = target_allowed & ~valid
        if not np.any(target_mask):
            break

        next_field, filled = _fill_from_neighbors(
            field=updated,
            source_mask=source_mask,
            target_mask=target_mask,
        )
        if filled == 0:
            break

        updated = next_field
        iterations += 1
        total_filled += filled
        if max_iterations is not None and iterations >= max_iterations:
            break

    return updated, iterations, total_filled


def _fill_from_neighbors(
    *,
    field: np.ndarray,
    source_mask: np.ndarray,
    target_mask: np.ndarray,
) -> tuple[np.ndarray, int]:
    supports = _build_supports(source_mask)
    denom = np.zeros(source_mask.shape, dtype=np.float32)
    active_supports = []
    for support, weight, dx, dy in supports:
        if not np.any(support):
            continue
        active_supports.append((support, weight, dx, dy))
        denom += support.astype(np.float32) * weight

    fillable = target_mask & (denom > _MIN_WEIGHT)
    filled = int(np.count_nonzero(fillable))
    if filled == 0:
        return field, 0

    numer = np.zeros_like(field, dtype=np.float32)
    for support, weight, dx, dy in active_supports:
        shifted = _shift_spatial(
            field,
            dx=dx,
            dy=dy,
            fill_value=np.float32(0.0),
        )
        numer += weight * support[..., np.newaxis].astype(np.float32) * shifted

    updated = np.array(field, copy=True, dtype=np.float32)
    updated[fillable] = numer[fillable] / denom[fillable, np.newaxis]
    return updated, filled


def _build_supports(
    mask: np.ndarray,
) -> list[tuple[np.ndarray, np.float32, int, int]]:
    im1_0 = _shift_spatial(mask, dx=-1, dy=0, fill_value=False)
    ip1_0 = _shift_spatial(mask, dx=1, dy=0, fill_value=False)
    i_0m1 = _shift_spatial(mask, dx=0, dy=-1, fill_value=False)
    i_0p1 = _shift_spatial(mask, dx=0, dy=1, fill_value=False)

    im2_0 = im1_0 & _shift_spatial(mask, dx=-2, dy=0, fill_value=False)
    ip2_0 = im1_0 & _shift_spatial(mask, dx=2, dy=0, fill_value=False)
    i_0m2 = i_0m1 & _shift_spatial(mask, dx=0, dy=-2, fill_value=False)
    i_0p2 = i_0m1 & _shift_spatial(mask, dx=0, dy=2, fill_value=False)

    im1_m1 = (im1_0 | i_0m1) & _shift_spatial(
        mask, dx=-1, dy=-1, fill_value=False
    )
    im1_p1 = (im1_0 | i_0p1) & _shift_spatial(
        mask, dx=-1, dy=1, fill_value=False
    )
    ip1_m1 = (ip1_0 | i_0m1) & _shift_spatial(
        mask, dx=1, dy=-1, fill_value=False
    )
    ip1_p1 = (ip1_0 | i_0p1) & _shift_spatial(
        mask, dx=1, dy=1, fill_value=False
    )

    im2_p1 = (im2_0 | im1_p1) & _shift_spatial(
        mask, dx=-2, dy=1, fill_value=False
    )
    im2_m1 = (im2_0 | im1_m1) & _shift_spatial(
        mask, dx=-2, dy=-1, fill_value=False
    )
    ip2_p1 = (ip2_0 | ip1_p1) & _shift_spatial(
        mask, dx=2, dy=1, fill_value=False
    )
    ip2_m1 = (ip2_0 | ip1_m1) & _shift_spatial(
        mask, dx=2, dy=-1, fill_value=False
    )
    im1_p2 = (im1_p1 | i_0p2) & _shift_spatial(
        mask, dx=-1, dy=2, fill_value=False
    )
    ip1_p2 = (ip1_p1 | i_0p2) & _shift_spatial(
        mask, dx=1, dy=2, fill_value=False
    )
    im1_m2 = (im1_m1 | i_0m2) & _shift_spatial(
        mask, dx=-1, dy=-2, fill_value=False
    )
    ip1_m2 = (ip1_m1 | i_0m2) & _shift_spatial(
        mask, dx=1, dy=-2, fill_value=False
    )

    im2_m2 = (im2_m1 | im1_m2) & _shift_spatial(
        mask, dx=-2, dy=-2, fill_value=False
    )
    im2_p2 = (im2_p1 | im1_p2) & _shift_spatial(
        mask, dx=-2, dy=2, fill_value=False
    )
    ip2_m2 = (ip2_m1 | ip1_m2) & _shift_spatial(
        mask, dx=2, dy=-2, fill_value=False
    )
    ip2_p2 = (ip2_p1 | ip1_p2) & _shift_spatial(
        mask, dx=2, dy=2, fill_value=False
    )

    im3_0 = im2_0 & _shift_spatial(mask, dx=-3, dy=0, fill_value=False)
    ip3_0 = im2_0 & _shift_spatial(mask, dx=3, dy=0, fill_value=False)
    i_0m3 = i_0m2 & _shift_spatial(mask, dx=0, dy=-3, fill_value=False)
    i_0p3 = i_0m2 & _shift_spatial(mask, dx=0, dy=3, fill_value=False)

    return [
        (im1_0, _W1, -1, 0),
        (ip1_0, _W1, 1, 0),
        (i_0m1, _W1, 0, -1),
        (i_0p1, _W1, 0, 1),
        (im2_0, _W2, -2, 0),
        (ip2_0, _W2, 2, 0),
        (i_0m2, _W2, 0, -2),
        (i_0p2, _W2, 0, 2),
        (im1_m1, _W11, -1, -1),
        (im1_p1, _W11, -1, 1),
        (ip1_m1, _W11, 1, -1),
        (ip1_p1, _W11, 1, 1),
        (im2_p1, _W21, -2, 1),
        (im2_m1, _W21, -2, -1),
        (ip2_p1, _W21, 2, 1),
        (ip2_m1, _W21, 2, -1),
        (im1_p2, _W21, -1, 2),
        (ip1_p2, _W21, 1, 2),
        (im1_m2, _W21, -1, -2),
        (ip1_m2, _W21, 1, -2),
        (im2_m2, _W22, -2, -2),
        (im2_p2, _W22, -2, 2),
        (ip2_m2, _W22, 2, -2),
        (ip2_p2, _W22, 2, 2),
        (im3_0, _W3, -3, 0),
        (ip3_0, _W3, 3, 0),
        (i_0m3, _W3, 0, -3),
        (i_0p3, _W3, 0, 3),
    ]


def _shift_spatial(
    array: np.ndarray,
    *,
    dx: int,
    dy: int,
    fill_value: bool | np.float32,
) -> np.ndarray:
    del fill_value

    pad = [
        (abs(dy), abs(dy)),
        (abs(dx), abs(dx)),
    ] + [(0, 0)] * (array.ndim - 2)
    padded = np.pad(array, pad_width=pad, mode='edge')

    y_start = abs(dy) + dy
    y_end = y_start + array.shape[0]
    x_start = abs(dx) + dx
    x_end = x_start + array.shape[1]

    return padded[y_start:y_end, x_start:x_end, ...]


def _get_missing_value(da: xr.DataArray) -> np.float32:
    for mapping in (da.encoding, da.attrs):
        for name in ('_FillValue', 'missing_value'):
            if name in mapping:
                value = np.asarray(mapping[name]).reshape(-1)[0]
                return np.float32(value)
    return np.float32(-1.0e6)


def _get_ocean_frac_name(ds_topo: xr.Dataset) -> str:
    for name in ('ocean_frac', 'mask_ocean', 'ocean_mask'):
        if name in ds_topo:
            return name
    raise KeyError(
        'Topography file is missing ocean fraction variable '
        "('ocean_frac', 'mask_ocean', or 'ocean_mask')."
    )


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
