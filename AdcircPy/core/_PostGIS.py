
class _PostGIS(object):

  def __init__(self):
    super(_PostGIS, self).__init__()

  def _parse_conn_string(self, conn_string):
    conn = "PG: "
    conn_string = conn_string.strip('PG: ').split()
    for item in conn_string:
      string = item.split("=")
      if string[0]=="schema":
        schema=string[1].strip("'")
      elif string[0]=="table":
        table=string[1].strip("'")
      else:
        conn+="{}={} ".format(string[0], string[1])
    if 'schema' not in locals():
      raise RuntimeError("No schema provided in connection string.")
    if 'table' not in locals():
      raise RuntimeError("No table provided in connection string.")
    return conn, schema, table

  def _get_PG_conn_string(cls, dbname, host, port, user, password):
    self._conn_string ="PG:"
    self._conn_string+="\""
    self._conn_string+="dbname='{}' ".format(dbname)
    self._conn_string+="host='{}' ".format(host)
    self._conn_string+="port='{}' ".format(port)
    self._conn_string+="user='{}' ".format(user)
    if len(password)>0:
      self._conn_string+="password='{}' ".format(password)
    self._conn_string+="\""