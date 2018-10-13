import os
import matplotlib
class AdcircPyEnvironment(object):
  if os.getenv('CIRCLECI') == 'true':
    matplotlib.use('Agg')
  def read_environment_variables(self):
    with open(os.getenv('HOME')+'/.adcpyrc') as lines:
      for line in lines:
        if 'export' in line and '#' not in line:
          line = line.replace('export', '').strip(' \n').split('=')
          os.environ[line[0]] = line[1]
    self.os=os
