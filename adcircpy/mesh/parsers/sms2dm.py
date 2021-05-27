import pathlib


def read(path):
    sms2dm = dict()
    with open(pathlib.Path(path), "r") as f:
        f.readline()
        while 1:
            line = f.readline().split()
            if len(line) == 0:
                break
            if line[0] in ["E3T", "E4Q"]:
                if line[0] not in sms2dm:
                    sms2dm[line[0]] = {}
                sms2dm[line[0]].update({line[1]: line[2:]})
            if line[0] == "ND":
                if line[0] not in sms2dm:
                    sms2dm[line[0]] = {}
                sms2dm[line[0]].update(
                    {line[1]: (list(map(float, line[2:-1])), float(line[-1]))}
                )
    return sms2dm


def write(sms2dm, path, overwrite=False):
    path = pathlib.Path(path)
    if path.is_file() and not overwrite:
        msg = "File exists, pass overwrite=True to allow overwrite."
        raise Exception(msg)
    with open(path, "w") as f:
        f.write(string(sms2dm))
    return 0  # for unittests


def string(sms2dm):
    f = graph(sms2dm)
    f += boundaries(sms2dm)
    return f


def graph(sms2dm):
    f = "MESH2D\n"
    # TODO: Make faster using np.array2string
    f += triangular_elements(sms2dm)
    f += quadrilateral_elements(sms2dm)
    f += nodes(sms2dm)
    return f


def nodes(sms2dm):
    assert all(int(id) > 0 for id in sms2dm["ND"])
    f = ""
    for id, (coords, value) in sms2dm["ND"].items():
        f += f"ND {int(id):d} "
        f += f"{coords[0]:<.16E} "
        f += f"{coords[1]:<.16E} "
        f += f"{value:<.16E}\n"
    return f


def boundaries(sms2dm):
    f = ""
    if "boundaries" in sms2dm.keys():
        for ibtype, bnds in sms2dm["boundaries"].items():
            for id, bnd in bnds.items():
                f += nodestring(bnd["indexes"])
    return f


def geom_string(geom_type, sms2dm):
    assert geom_type in ["E3T", "E4Q", "E6T", "E8Q", "E9Q"]
    assert all(int(id) > 0 for id in sms2dm[geom_type])
    f = ""
    for id, geom in sms2dm[geom_type].items():
        f += f"{geom_type} {id} "
        for j in range(len(geom)):
            f += f"{geom[j]} "
        f += "\n"
    return f


def nodestring(geom):
    f = "NS "
    for i in range(len(geom) - 1):
        f += f"{geom[i]} "
    f += f"-{geom[-1]}\n"
    return f


def nodestrings(geom):
    pass


def triangular_elements(geom):
    f = ""
    if geom is not None:
        f += geom_string("E3T", geom)
    return f


def quadrilateral_elements(geom):
    f = ""
    if geom is not None:
        f += geom_string("E4Q", geom)
    return f
