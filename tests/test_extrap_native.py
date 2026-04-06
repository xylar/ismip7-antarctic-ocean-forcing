import numpy as np

from i7aof.extrap.horizontal import horizontal_extrapolation
from i7aof.extrap.vertical import vertical_extrapolation


def test_horizontal_extrapolation_matches_reference() -> None:
    miss = np.float32(-9999.0)
    z = np.array([-10.0, -30.0, -60.0], dtype=np.float32)

    data = np.full((3, 6, 6, 2), miss, dtype=np.float32)
    for z_index, base in enumerate((0.0, 100.0, 200.0)):
        data[z_index, 2, 1, :] = np.array([base + 1.0, base + 2.0])
        data[z_index, 2, 4, :] = np.array([base + 10.0, base + 20.0])
    data[0, 1, 1, :] = np.array([4.0, 5.0], dtype=np.float32)
    data[1, 4, 4, :] = np.array([115.0, 125.0], dtype=np.float32)

    basin = np.zeros((6, 6), dtype=np.int64)
    basin[:, 3:] = 1

    bed = np.full((6, 6), -100.0, dtype=np.float32)
    bed[1, 2] = -5.0
    bed[4, 3] = -15.0

    ocean_frac = np.ones((6, 6), dtype=np.float32)
    ocean_frac[2, 2] = 0.0
    ocean_frac[3, 3] = 0.0

    expected = _horizontal_reference(
        data=data,
        z=z,
        basin_number=basin,
        bed=bed,
        ocean_frac=ocean_frac,
        miss=miss,
    )
    actual = horizontal_extrapolation(
        data=data,
        z=z,
        basin_number=basin,
        bed=bed,
        ocean_frac=ocean_frac,
        miss=miss,
    )

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1.0e-6)


def test_vertical_extrapolation_fills_downward() -> None:
    miss = np.float32(-9999.0)
    data = np.array(
        [
            [[[miss], [1.0]], [[miss], [miss]]],
            [[[3.0], [2.0]], [[miss], [4.0]]],
            [[[miss], [miss]], [[5.0], [miss]]],
        ],
        dtype=np.float32,
    )

    expected = np.array(
        [
            [[[3.0], [1.0]], [[miss], [4.0]]],
            [[[3.0], [2.0]], [[miss], [4.0]]],
            [[[3.0], [2.0]], [[5.0], [4.0]]],
        ],
        dtype=np.float32,
    )

    actual = vertical_extrapolation(data=data, miss=miss)
    np.testing.assert_array_equal(actual, expected)


def test_native_file_pipeline_round_trip(tmp_path) -> None:
    try:
        import cftime
        import xarray as xr

        from i7aof.extrap.horizontal import extrapolate_horizontal
        from i7aof.extrap.vertical import extrapolate_vertical
        from i7aof.io import read_dataset, write_netcdf
    except ModuleNotFoundError:
        return

    miss = np.float32(-9999.0)
    z = np.array([-10.0, -30.0], dtype=np.float32)
    y = np.array([0.0, 1.0, 2.0], dtype=np.float32)
    x = np.array([0.0, 1.0, 2.0], dtype=np.float32)
    times = [
        cftime.DatetimeNoLeap(1850, 1, 1),
        cftime.DatetimeNoLeap(1850, 2, 1),
    ]
    time_bnds = np.array(
        [
            [
                cftime.DatetimeNoLeap(1850, 1, 1),
                cftime.DatetimeNoLeap(1850, 2, 1),
            ],
            [
                cftime.DatetimeNoLeap(1850, 2, 1),
                cftime.DatetimeNoLeap(1850, 3, 1),
            ],
        ],
        dtype=object,
    )

    values = np.full((2, 2, 3, 3), miss, dtype=np.float32)
    values[:, :, 1, 1] = np.array([[1.0, 2.0], [10.0, 11.0]], dtype=np.float32)

    ds_in = xr.Dataset(
        data_vars={
            'ct': (('time', 'z_extrap', 'y', 'x'), values),
            'time_bnds': (('time', 'bnds'), time_bnds),
        },
        coords={
            'time': ('time', times),
            'z_extrap': ('z_extrap', z),
            'y': ('y', y),
            'x': ('x', x),
        },
        attrs={'history': 'synthetic extrapolation test'},
    )
    ds_in['time'].attrs.update({'bounds': 'time_bnds', 'calendar': 'noleap'})
    ds_in['ct'].attrs['missing_value'] = miss

    basin = xr.Dataset(
        {'basinNumber': (('y', 'x'), np.zeros((3, 3), dtype=np.int64))}
    )
    topo = xr.Dataset(
        {
            'bed': (('y', 'x'), np.full((3, 3), -100.0, dtype=np.float32)),
            'ocean_frac': (
                ('y', 'x'),
                np.ones((3, 3), dtype=np.float32),
            ),
        }
    )

    in_path = tmp_path / 'input.nc'
    horiz_path = tmp_path / 'horizontal.nc'
    vert_path = tmp_path / 'vertical.nc'
    basin_path = tmp_path / 'basin.nc'
    topo_path = tmp_path / 'topo.nc'

    write_netcdf(ds_in, in_path, has_fill_values=['ct'], progress_bar=False)
    basin.to_netcdf(basin_path)
    topo.to_netcdf(topo_path)

    extrapolate_horizontal(
        in_path=str(in_path),
        out_path=str(horiz_path),
        basin_path=str(basin_path),
        topo_path=str(topo_path),
        variable='ct',
    )
    extrapolate_vertical(
        in_path=str(horiz_path),
        out_path=str(vert_path),
        variable='ct',
    )

    with read_dataset(vert_path) as ds_out:
        assert ds_out['ct'].dims == ('time', 'z_extrap', 'y', 'x')
        assert np.isfinite(ds_out['ct'].isel(time=0, z_extrap=0)).any()


def _horizontal_reference(
    *,
    data: np.ndarray,
    z: np.ndarray,
    basin_number: np.ndarray,
    bed: np.ndarray,
    ocean_frac: np.ndarray,
    miss: np.float32,
) -> np.ndarray:
    result = np.array(data, copy=True, dtype=np.float32)
    basin_ids = range(int(np.max(basin_number)) + 1)

    for z_index, z_level in enumerate(z):
        plane = result[z_index]
        below_level = bed < z_level

        for basin_id in basin_ids:
            basin_mask = basin_number == basin_id
            plane = _iterate_reference(
                plane=plane,
                source_allowed=below_level & basin_mask,
                target_allowed=below_level
                & basin_mask
                & (ocean_frac > np.float32(0.01)),
                miss=miss,
            )

        for basin_id in basin_ids:
            basin_mask = basin_number == basin_id
            plane = _iterate_reference(
                plane=plane,
                source_allowed=below_level & basin_mask,
                target_allowed=below_level & basin_mask,
                miss=miss,
            )

        plane = _iterate_reference(
            plane=plane,
            source_allowed=np.ones(plane.shape[:2], dtype=bool),
            target_allowed=np.ones(plane.shape[:2], dtype=bool),
            miss=miss,
            max_iterations=None if z_index == 1 else 5,
        )
        result[z_index] = plane

    return result


def _iterate_reference(
    *,
    plane: np.ndarray,
    source_allowed: np.ndarray,
    target_allowed: np.ndarray,
    miss: np.float32,
    max_iterations: int | None = None,
) -> np.ndarray:
    updated = np.array(plane, copy=True, dtype=np.float32)
    iterations = 0

    while True:
        valid = updated[..., 0] != miss
        source_mask = source_allowed & valid
        target_mask = target_allowed & ~valid
        if not np.any(target_mask):
            break

        next_plane, filled = _fill_once_reference(
            plane=updated,
            source_mask=source_mask,
            target_mask=target_mask,
            miss=miss,
        )
        if filled == 0:
            break

        updated = next_plane
        iterations += 1
        if max_iterations is not None and iterations >= max_iterations:
            break

    return updated


def _fill_once_reference(
    *,
    plane: np.ndarray,
    source_mask: np.ndarray,
    target_mask: np.ndarray,
    miss: np.float32,
) -> tuple[np.ndarray, int]:
    ny, nx, nt = plane.shape
    updated = np.array(plane, copy=True, dtype=np.float32)
    filled = 0

    for y_index in range(ny):
        for x_index in range(nx):
            if not target_mask[y_index, x_index]:
                continue

            denom = np.float32(0.0)
            numer = np.zeros(nt, dtype=np.float32)
            for flag, weight, src_y, src_x in _reference_terms(
                source_mask, y_index, x_index
            ):
                if not flag:
                    continue
                denom += weight
                numer += weight * plane[src_y, src_x, :]

            if denom > np.float32(0.1):
                updated[y_index, x_index, :] = numer / denom
                filled += 1

    return updated, filled


def _reference_terms(
    mask: np.ndarray, y_index: int, x_index: int
) -> list[tuple[bool, np.float32, int, int]]:
    ny, nx = mask.shape
    xm1 = max(x_index - 1, 0)
    xm2 = max(x_index - 2, 0)
    xm3 = max(x_index - 3, 0)
    xp1 = min(x_index + 1, nx - 1)
    xp2 = min(x_index + 2, nx - 1)
    xp3 = min(x_index + 3, nx - 1)
    ym1 = max(y_index - 1, 0)
    ym2 = max(y_index - 2, 0)
    ym3 = max(y_index - 3, 0)
    yp1 = min(y_index + 1, ny - 1)
    yp2 = min(y_index + 2, ny - 1)
    yp3 = min(y_index + 3, ny - 1)

    xm1_0 = bool(mask[y_index, xm1])
    xp1_0 = bool(mask[y_index, xp1])
    y_0m1 = bool(mask[ym1, x_index])
    y_0p1 = bool(mask[yp1, x_index])

    xm2_0 = xm1_0 and bool(mask[y_index, xm2])
    xp2_0 = xm1_0 and bool(mask[y_index, xp2])
    y_0m2 = y_0m1 and bool(mask[ym2, x_index])
    y_0p2 = y_0m1 and bool(mask[yp2, x_index])

    xm1_m1 = (xm1_0 or y_0m1) and bool(mask[ym1, xm1])
    xm1_p1 = (xm1_0 or y_0p1) and bool(mask[yp1, xm1])
    xp1_m1 = (xp1_0 or y_0m1) and bool(mask[ym1, xp1])
    xp1_p1 = (xp1_0 or y_0p1) and bool(mask[yp1, xp1])

    xm2_p1 = (xm2_0 or xm1_p1) and bool(mask[yp1, xm2])
    xm2_m1 = (xm2_0 or xm1_m1) and bool(mask[ym1, xm2])
    xp2_p1 = (xp2_0 or xp1_p1) and bool(mask[yp1, xp2])
    xp2_m1 = (xp2_0 or xp1_m1) and bool(mask[ym1, xp2])
    xm1_p2 = (xm1_p1 or y_0p2) and bool(mask[yp2, xm1])
    xp1_p2 = (xp1_p1 or y_0p2) and bool(mask[yp2, xp1])
    xm1_m2 = (xm1_m1 or y_0m2) and bool(mask[ym2, xm1])
    xp1_m2 = (xp1_m1 or y_0m2) and bool(mask[ym2, xp1])

    xm2_m2 = (xm2_m1 or xm1_m2) and bool(mask[ym2, xm2])
    xm2_p2 = (xm2_p1 or xm1_p2) and bool(mask[yp2, xm2])
    xp2_m2 = (xp2_m1 or xp1_m2) and bool(mask[ym2, xp2])
    xp2_p2 = (xp2_p1 or xp1_p2) and bool(mask[yp2, xp2])

    xm3_0 = xm2_0 and bool(mask[y_index, xm3])
    xp3_0 = xm2_0 and bool(mask[y_index, xp3])
    y_0m3 = y_0m2 and bool(mask[ym3, x_index])
    y_0p3 = y_0m2 and bool(mask[yp3, x_index])

    return [
        (xm1_0, np.float32(0.946), y_index, xm1),
        (xp1_0, np.float32(0.946), y_index, xp1),
        (y_0m1, np.float32(0.946), ym1, x_index),
        (y_0p1, np.float32(0.946), yp1, x_index),
        (xm2_0, np.float32(0.801), y_index, xm2),
        (xp2_0, np.float32(0.801), y_index, xp2),
        (y_0m2, np.float32(0.801), ym2, x_index),
        (y_0p2, np.float32(0.801), yp2, x_index),
        (xm1_m1, np.float32(0.895), ym1, xm1),
        (xm1_p1, np.float32(0.895), yp1, xm1),
        (xp1_m1, np.float32(0.895), ym1, xp1),
        (xp1_p1, np.float32(0.895), yp1, xp1),
        (xm2_p1, np.float32(0.757), yp1, xm2),
        (xm2_m1, np.float32(0.757), ym1, xm2),
        (xp2_p1, np.float32(0.757), yp1, xp2),
        (xp2_m1, np.float32(0.757), ym1, xp2),
        (xm1_p2, np.float32(0.757), yp2, xm1),
        (xp1_p2, np.float32(0.757), yp2, xp1),
        (xm1_m2, np.float32(0.757), ym2, xm1),
        (xp1_m2, np.float32(0.757), ym2, xp1),
        (xm2_m2, np.float32(0.641), ym2, xm2),
        (xm2_p2, np.float32(0.641), yp2, xm2),
        (xp2_m2, np.float32(0.641), ym2, xp2),
        (xp2_p2, np.float32(0.641), yp2, xp2),
        (xm3_0, np.float32(0.606), y_index, xm3),
        (xp3_0, np.float32(0.606), y_index, xp3),
        (y_0m3, np.float32(0.606), ym3, x_index),
        (y_0p3, np.float32(0.606), yp3, x_index),
    ]
