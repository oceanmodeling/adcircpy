#! /usr/bin/env python
import unittest
from datetime import datetime, timedelta
from AdcircPyTests import AdcircPyEnvironment
from AdcircPy import AdcircPy
from AdcircPy.core import PBS, ServerConfiguration


class TestServerConfiguration(AdcircPyEnvironment, unittest.TestCase):
  
  def setUp(self):
    super(TestServerConfiguration, self).__init__()
    self.read_environment_variables()
    AdcircMesh = AdcircPy.read_mesh(fort14=self.os.getenv('FORT14_PATH'),
                                    fort13=self.os.getenv('FORT13_PATH'))
    start_time = datetime.now()
    end_time = start_time + timedelta(days=5)
    self.TidalRun = AdcircMesh.TidalRun(start_time, end_time)

  def test_ServerConfiguration(self):
    walltime = timedelta(hours=2, minutes=38)
    numprocs = 666
    qsub = PBS('account_name', walltime, numprocs, email='test@email.test')
    module_list = ['intel/15.6.233', 'impi/5.0.3.048', 'hdf5parallel', 'netcdf-hdf5parallel']
    serverConfig = self.TidalRun.ServerConfiguration('./test_run_staging', self.os.getenv('ADCIRC_BINARIES_PATH'),
                                                     PBS=qsub,
                                                     module_list=module_list,
                                                     # environment_file='~/opt/environment.sh'
                                                     )
    serverConfig.printf()

if __name__ == '__main__':
  unittest.main()



