import pathlib
PARENT = pathlib.Path(__file__).parent.absolute()
FORT14 = PARENT / "data/NetCDF_Shinnecock_Inlet/fort.14"
# fetch shinnecock inlet test data
if not FORT14.is_file():
    import urllib
    import tempfile
    import tarfile
    url = "https://www.dropbox.com/s/1wk91r67cacf132/"
    url += "NetCDF_shinnecock_inlet.tar.bz2?dl=1"
    g = urllib.request.urlopen(url)
    tmpfile = tempfile.NamedTemporaryFile()
    with open(tmpfile.name, 'b+w') as f:
        f.write(g.read())
    with tarfile.open(tmpfile.name, "r:bz2") as tar:
        tar.extractall(PARENT / "data/NetCDF_Shinnecock_Inlet/")
