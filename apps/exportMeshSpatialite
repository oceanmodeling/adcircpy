def surface_to_spatialite(self, path, table, levels, overwrite=False, epsg=4326):

    conn = sqlite3.connect(path)
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')
    cur = conn.cursor()
    
    # Check if SpatialMetadata exists, if not, create it.
    cur.execute("PRAGMA table_info('spatial_ref_sys');")
    if len(cur.fetchall()) == 0:
        cur.execute("SELECT InitSpatialMetadata()")
    

    cur.execute("PRAGMA table_info({});".format(table))
    table_info = cur.fetchall()
    
    if overwrite==True or len(table_info)==0:
        cur.execute("DROP TABLE IF EXISTS {};".format(table))
        cur.execute("CREATE TABLE {} (id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT);".format(table))
    
    elif overwrite==False and len(table_info)>0:
        raise IOError("Table with name '{}' already exists in database. Set keyword argument overwrite=True if you wish to overwrite the existing table.".format(table))
    
    # cur.execute("PRAGMA table_info({});".format(table))
    # print(cur.fetchall())
    sql = "SELECT DiscardGeometryColumn('{}', 'geom')".format(table)
    cur.execute(sql)
    print(cur.fetchall())
    
    sql = "SELECT AddGeometryColumn('{}', 'geom', {}, 'MULTIPOLYGON', 'XY', 1)".format(table, epsg)
    # print(sql)
    cur.execute(sql)
    
    print(cur.fetchall())
    