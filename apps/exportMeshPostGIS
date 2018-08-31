#! /usr/bin/env python

def deploy_to_PostGIS(self, dbname, **kwargs):
    """
    Writes an ADCIRCgrid object to a PostGIS database.
    Database must be created manually on the server prior to deployment.
    """
    user        = kwargs.pop('user', 'postgres')
    password    = kwargs.pop('password', True)
    host        = kwargs.pop('host', 'localhost')
    port        = kwargs.pop('port', 5432)
    progressbar = kwargs.pop('progressbar', False)
    overwrite   = kwargs.pop('overwrite', False)
    schema      = kwargs.pop('schema', 'ADCIRC')
    boundaries_only = kwargs.pop('boundaries_only', False)


    if password==True:
        password = getpass.getpass('Password: ')
        

    con = psycopg2.connect(dbname=dbname, user=user, password=password, host=host, port=port)
    cur  = con.cursor()
    
    # Create schema to store grid
    cur.execute('CREATE SCHEMA IF NOT EXISTS {};'.format(schema))
    con.commit()
 
    # Ocean boundaries injection
    cur.execute("SELECT to_regclass('{}.ocean_boundaries');".format(schema))
    table_existance = cur.fetchall()[0][0]
    if table_existance is None or overwrite==True:
        if table_existance is not None and overwrite==True:
            cur.execute('DROP TABLE {}.ocean_boundaries;'.format(schema))
            con.commit()
        cur.execute('CREATE TABLE IF NOT EXISTS {}.ocean_boundaries (id BIGSERIAL PRIMARY KEY, geom geometry(Linestring, {}));'.format(schema, self.epsg))
        if self.oceanBoundaries is not None:
            for boundary in self.oceanBoundaries:
                geom = "ST_GeomFromText('LINESTRING("
                for index in boundary:
                    geom += '{:f} {:f}, '.format(self.x[index], self.y[index])
                geom = geom[:-2]
                geom += ")', {})".format(self.epsg)
                cur.execute("INSERT INTO {}.ocean_boundaries (geom) VALUES ({});".format(schema, geom))
            con.commit()
            cur.execute("CREATE INDEX sidx_ocean_boundaries_geom ON {}.ocean_boundaries USING GIST (geom);".format(schema))
            con.commit()
        
    # land boundaries injection
    cur.execute("SELECT to_regclass('{}.land_boundaries');".format(schema))
    table_existance = cur.fetchall()[0][0]
    if table_existance is None or overwrite==True:
        if table_existance is not None and overwrite==True:
            cur.execute('DROP TABLE {}.land_boundaries;'.format(schema))
            con.commit()
        cur.execute('CREATE TABLE IF NOT EXISTS {}.land_boundaries (id BIGSERIAL PRIMARY KEY, geom geometry(Linestring, {}), btype SMALLINT NOT NULL);'.format(schema, self.epsg))
        if self.landBoundaries is not None:
            for indices, btype in self.landBoundaries:
                geom = "ST_GeomFromText('LINESTRING("
                for index in indices:
                    geom += '{:f} {:f}, '.format(self.x[index], self.y[index])
                geom = geom[:-2]
                geom += ")', {})".format(self.epsg)
                cur.execute("INSERT INTO {}.land_boundaries (geom, btype) VALUES ({}, {});".format(schema, geom, btype))
            con.commit()
            cur.execute("CREATE INDEX sidx_land_boundaries_geom ON {}.land_boundaries USING GIST (geom);".format(schema))
            con.commit()
        
    # inner boundaries injection
    cur.execute("SELECT to_regclass('{}.inner_boundaries');")
    table_existance = cur.fetchall()[0][0]
    if table_existance is None or overwrite==True:
        if table_existance is not None and overwrite==True:
            cur.execute('DROP TABLE {}.inner_boundaries')
            con.commit()
        cur.execute('CREATE TABLE IF NOT EXISTS {}.inner_boundaries (id BIGSERIAL PRIMARY KEY, geom geometry(Linestring, {}), btype SMALLINT NOT NULL);'.format(schema, self.epsg))        
        for indices, btype in self.innerBoundaries:
            geom = "ST_GeomFromText('LINESTRING("
            for index in indices:
                geom += '{:f} {:f}, '.format(self.x[index], self.y[index])
            geom += '{:f} {:f}'.format(self.x[indices[0]], self.y[indices[0]])
            geom += ")', {})".format(self.epsg)
            cur.execute("INSERT INTO {}.inner_boundaries (geom, btype) VALUES ({}, {});".format(schema, geom, btype))
        con.commit()
        
    if self.weirBoundaries is not None:
        cur.execute("SELECT to_regclass('{}.weir_boundaries');".format(schema))
        table_existance = cur.fetchall()[0][0]
        if table_existance is None or overwrite==True:
            if table_existance is not None and overwrite==True:
                cur.execute('DROP TABLE {}.weir_boundaries;'.format(schema))
                con.commit()
        sql = "CREATE TABLE IF NOT EXISTS {}.weir_boundaries (id SMALLSERIAL PRIMARY KEY, faces geometry(MultiLinestring, {}), ".format(schema, self.epsg)
        sql+= "heights geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "supercritical_flow_coefficient geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "subcritical_flow_coefficient geometry(MultipointM, {}), ".format(self.epsg)
        sql+= " btype smallint);"
        cur.execute(sql)
        con.commit()
        for boundary in self.weirBoundaries:
            faces = "ST_GeomFromText('MULTILINESTRING(("
            heights = "ST_GeomFromText('MULTIpointM("
            supercritical_flow_coefficient = "ST_GeomFromText('MULTIpointM("
            subcritical_flow_coefficient = "ST_GeomFromText('MULTIpointM("
            for i, index in enumerate(boundary['front_face']):
                faces+= "{:f} {:f}, ".format(self.x[index], self.y[index])
                heights += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['height'][i])
                supercritical_flow_coefficient += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['supercritical_flow_coefficient'][i])
                subcritical_flow_coefficient += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['subcritical_flow_coefficient'][i])
            faces = faces[:-2]+"), ("
            for i, index in enumerate(boundary['back_face']):
                faces+= "{:f} {:f}, ".format(self.x[index], self.y[index])
            faces = faces[:-2]+"))', {})".format(self.epsg)
            heights = heights[:-2]+")', {})".format(self.epsg)
            supercritical_flow_coefficient = supercritical_flow_coefficient[:-2]+")', {})".format(self.epsg)
            subcritical_flow_coefficient = subcritical_flow_coefficient[:-2]+")', {})".format(self.epsg)
            sql = "INSERT INTO {}.weir_boundaries ".format(schema)
            sql+= "(faces, heights, supercritical_flow_coefficient, subcritical_flow_coefficient, btype) "
            sql+= "VALUES ({}, {}, {}, {}, {});".format(faces, heights, supercritical_flow_coefficient, subcritical_flow_coefficient, boundary['btype'])
            cur.execute(sql)
            con.commit()
        cur.execute("CREATE INDEX sidx_weir_boundaries_faces ON {}.weir_boundaries USING GIST (faces);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_heights ON {}.weir_boundaries USING GIST (heights);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_supercritical_flow_coefficient ON {}.weir_boundaries USING GIST (supercritical_flow_coefficient);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_subcritical_flow_coefficient ON {}.weir_boundaries USING GIST (subcritical_flow_coefficient);".format(schema))
        con.commit()
    
    if self.culvertBoundaries is not None:
        cur.execute("SELECT to_regclass('{}.culvert_boundaries');".format(schema))
        table_existance = cur.fetchall()[0][0]
        if table_existance is None or overwrite==True:
            if table_existance is not None and overwrite==True:
                cur.execute('DROP TABLE {}.culvert_boundaries;')
                con.commit()
        sql = "CREATE TABLE IF NOT EXISTS {}.culvert_boundaries (id SMALLSERIAL PRIMARY KEY, faces geometry(MultiLinestring, {}), ".format(schema, self.epsg)
        sql+= "heights geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "supercritical_flow_coefficient geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "subcritical_flow_coefficient geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "cross_barrier_pipe_height geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "friction_factor geometry(MultipointM, {}), ".format(self.epsg)
        sql+= "pipe_diameter geometry(MultipointM, {}), ".format(self.epsg)
        sql+= " btype smallint);"
        cur.execute(sql)
        con.commit()
        for boundary in self.weirBoundaries:
            faces = "ST_GeomFromText('MULTILINESTRING(("
            heights = "ST_GeomFromText('MULTIpointM("
            supercritical_flow_coefficient = "ST_GeomFromText('MULTIpointM("
            subcritical_flow_coefficient = "ST_GeomFromText('MULTIpointM("
            cross_barrier_pipe_height = "ST_GeomFromText('MULTIpointM("
            friction_factor = "ST_GeomFromText('MULTIpointM("
            pipe_diameter = "ST_GeomFromText('MULTIpointM("
            for i, index in enumerate(boundary['front_face']):
                faces+= "{:f} {:f}, ".format(self.x[index], self.y[index])
                heights += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['height'][i])
                supercritical_flow_coefficient += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['supercritical_flow_coefficient'][i])
                subcritical_flow_coefficient += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['subcritical_flow_coefficient'][i])
                cross_barrier_pipe_height += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['cross_barrier_pipe_height'][i])
                friction_factor += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['friction_factor'][i])
                pipe_diameter += "{:f} {:f} {:f}, ".format(self.x[index], self.y[index], boundary['pipe_diameter'][i])
            faces = faces[:-2]+"), ("
            for i, index in enumerate(boundary['back_face']):
                faces+= "{:f} {:f}, ".format(self.x[index], self.y[index])
            faces = faces[:-2]+"))', {})".format(self.epsg)
            heights = heights[:-2]+")', {})".format(self.epsg)
            supercritical_flow_coefficient = supercritical_flow_coefficient[:-2]+")', {})".format(self.epsg)
            subcritical_flow_coefficient = subcritical_flow_coefficient[:-2]+")', {})".format(self.epsg)
            cross_barrier_pipe_height = cross_barrier_pipe_height[:-2]+")', {})".format(self.epsg)
            friction_factor = friction_factor[:-2]+")', {})".format(self.epsg)
            pipe_diameter = pipe_diameter[:-2]+")', {})".format(self.epsg)
            sql = "INSERT INTO {}.weir_boundaries ".format(schema)
            sql+= "(faces, heights, supercritical_flow_coefficient, subcritical_flow_coefficient, cross_barrier_pipe_height, friction_factor, pipe_diameter,btype) "
            sql+= "VALUES ({}, {}, {}, {}, {}, {}, {}, {});".format(faces, heights, supercritical_flow_coefficient, subcritical_flow_coefficient, cross_barrier_pipe_height, friction_factor, pipe_diameter, boundary['btype'])
            cur.execute(sql)
            con.commit()
        cur.execute("CREATE INDEX sidx_weir_boundaries_faces ON {}.weir_boundaries USING GIST (faces);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_heights ON {}.weir_boundaries USING GIST (heights);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_supercritical_flow_coefficient ON {}.weir_boundaries USING GIST (supercritical_flow_coefficient);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_subcritical_flow_coefficient ON {}.weir_boundaries USING GIST (subcritical_flow_coefficient);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_cross_barrier_pipe_height ON {}.weir_boundaries USING GIST (cross_barrier_pipe_height);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_friction_factor ON {}.weir_boundaries USING GIST (friction_factor);".format(schema))
        cur.execute("CREATE INDEX sidx_weir_boundaries_pipe_diameter ON {}.weir_boundaries USING GIST (pipe_diameter);".format(schema))
        con.commit()
        
    # Triangle injection
    if boundaries_only==False:
        if progressbar==True:
            pbar = tqdm.tqdm(total=self.elements.shape[0])
        cur.execute("SELECT to_regclass('{}.elements');".format(schema))
        table_existance = cur.fetchall()[0][0]
        if table_existance is None or overwrite==True:
            if table_existance is not None and overwrite==True:
                cur.execute('DROP TABLE {}.elements;'.format(schema))
                con.commit()
            cur.execute('CREATE TABLE IF NOT EXISTS {}.elements (id BIGSERIAL PRIMARY KEY, geom geometry(LinestringZ, {}));'.format(schema, self.epsg))
            for row in self.elements:
                geom = "ST_GeomFromText('LinestringZ("
                geom += "%f %f %f, " % (self.x[row][0], self.y[row][0], self.values[row][0])
                geom += "%f %f %f, " % (self.x[row][1], self.y[row][1], self.values[row][1])
                geom += "%f %f %f, " % (self.x[row][2], self.y[row][2], self.values[row][2])
                geom += "%f %f %f" % (self.x[row][0], self.y[row][0], self.values[row][0])
                geom += ")', {})".format(self.epsg)
                cur.execute("INSERT INTO {}.elements (geom) VALUES ({})".format(schema, geom))
                if progressbar==True:
                    pbar.update(1)
            con.commit()
            cur.exectue("CREATE INDEX sidx_elements_geom ON {}.elements USING GIST (geom);".format(schema))
            con.commit()

    # inject depth nodal attribute        
    # if self.fort13 is not None:
        # cur.execute("SELECT to_regclass('grid.nodal_attributes');")
        # table_existance = cur.fetchall()[0][0]
        # if table_existance is None or overwrite==True:
            # if table_existance is not None and overwrite==True:
                # cur.execute('DROP TABLE grid.nodal_attributes;')
                # con.commit()
            # cur.execute("CREATE TABLE IF NOT EXISTS grid.nodal_attributes (id integer PRIMARY KEY, geom geometry(Point, {}), depth DOUBLE PRECISION);".format(self.epsg))
            # con.commit()
            # for index in self.nodeID:
                # geom = "ST_GeomFromText('POINT("
                # geom += "{:f} {:f})', {})".format(self.x[index], self.y[index], self.epsg)
                # cur.execute("INSERT INTO grid.nodal_attributes (id, geom, depth) VALUES ({}, {}, {});".format(index, geom, self.values[index]))
                # if progressbar ==True:
                    # pbar.update(1)
            # con.commit()
        

    
    # cur.execute("SELECT exists(select schema_name FROM information_schema.schemata WHERE schema_name = 'topoelements');")
    # schema_existance = cur.fetchone()[0]
    # if schema_existance==True and overwrite==True:
        # cur.execute("SELECT topology.DropTopology('topoelements')")
        # con.commit()
    # if schema_existance==False or overwrite==True:
        # cur.execute("SELECT topology.CreateTopology('topoelements', {}, 0, 't');".format(self.epsg))
        # con.commit()
        # for row in self.elements:       
            # geom = "ST_GeomFromText('LinestringZ("
            # geom += "%f %f %f, " % (self.x[row][0], self.y[row][0], self.values[row][0])
            # geom += "%f %f %f, " % (self.x[row][1], self.y[row][1], self.values[row][1])
            # geom += "%f %f %f, " % (self.x[row][2], self.y[row][2], self.values[row][2])
            # geom += "%f %f %f" % (self.x[row][0], self.y[row][0], self.values[row][0])
            # geom += ")', {})".format(self.epsg)
            # cur.execute("SELECT topology.TopoGeo_AddLinestring('topoelements', {});".format(geom))
            # pbar.update(1)
        # con.commit()   

        
            
    
    con.commit()
    con.close()
