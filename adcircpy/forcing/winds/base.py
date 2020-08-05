from pyproj import CRS, Geod, Transformer


class WindForcing:
    pass


def ellipsoidal_distance(point_a: (float, float), point_b: (float, float), crs_a: CRS,
                         crs_b: CRS = None) -> float:
    if crs_b is not None:
        transformer = Transformer.from_crs(crs_b, crs_a)
        point_b = transformer.transform(*point_b)
    datum_json = crs_a.datum.to_json_dict()
    ellipsoid = Geod(a=datum_json['ellipsoid']['semi_major_axis'],
                     rf=datum_json['ellipsoid']['inverse_flattening'])
    return ellipsoid.inv(point_a[0], point_a[1], point_b[0], point_b[1])[2]
