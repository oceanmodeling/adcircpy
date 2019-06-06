# global imports
import numpy as np

# local imports
from AdcircPy.Mesh import UnstructuredMesh
from AdcircPy.Mesh.UnstructuredMesh import NodalAttributes
from AdcircPy.Model._TidalRun import _TidalRun
from AdcircPy.Model._BestTrackRun import _BestTrackRun

# unittest imports
import os
import unittest


class AdcircMesh(UnstructuredMesh):

    def __init__(self, fort14, SpatialReference, vertical_datum, fort13=None,
                 datum_mesh=None):
        fort14 = self.parse_fort14(fort14)
        x = np.asarray(fort14.pop('x'))
        y = np.asarray(fort14.pop('y'))
        xy = np.vstack([x, y]).T
        super(AdcircMesh, self).__init__(xy,
                                         SpatialReference=SpatialReference,
                                         vertical_datum=vertical_datum,
                                         **fort14)
        if fort13 is not None:
            self.set_fort13(fort13)
        if datum_mesh is not None:
            self.convert_datum(datum_mesh)

    def get_TidalRun(self, start_date, end_date, spinup_days=7,
                     constituents='all', **fort15):
        return _TidalRun(self, start_date, end_date, spinup_days=spinup_days,
                         constituents=constituents, **fort15)

    def get_BestTrackRun(self, storm_id, start_date=None, end_date=None,
                         spinup_days=7, constituents='all', **fort15):
        return _BestTrackRun(self, storm_id, start_date=start_date,
                             end_date=end_date, spinup_days=spinup_days,
                             constituents=constituents, **fort15)

    def TidalRun(self, start_date, end_date, spinup_days=7.,
                 constituents='all', **fort15):
        return self.get_TidalRun(start_date, end_date, spinup_days,
                                 constituents=constituents, **fort15)

    def BestTrackRun(self, start_date=None, end_date=None, spinup_days=7.,
                     constituents='all', **fort15):
        return self.get_BestTrackRun(start_date, end_date, spinup_days,
                                     constituents=constituents, **fort15)

    def add_nodal_attribute(self, attribute_name, full_values, defaults=None):
        self.NodalAttributes._add_attribute(self, attribute_name, full_values,
                                            defaults=None)

    def set_fort13(self, path, coldstart_attributes=[],
                   hotstart_attributes=[], use_all=True):
        self.__NodalAttributes = NodalAttributes.from_fort13(
                                            self, path, coldstart_attributes,
                                            hotstart_attributes, use_all)
        self.__fort13 = path

    def set_nodal_attribute_state(self, *attributes, coldstart=True,
                                  hotstart=True):
        self.NodalAttributes.set_attribute_state(*attributes,
                                                 coldstart=coldstart,
                                                 hotstart=hotstart)

    def auto_generate_TAU0(self, **kwargs):
        self.set_TAU0(3, **kwargs)

    def set_TAU0(self, TAU0, **kwargs):
        self.NodalAttributes.set_TAU0(self, TAU0, **kwargs)

    def make_plot(self, cmap='topobathy', **kwargs):
        super(AdcircMesh, self).make_plot(cmap='topobathy', **kwargs)

    def write_fort14(self, path, filename='fort.14'):
        raise NotImplementedError
        fort14 = list()
        fort14.append("{}\n".format(self.description))
        fort14.append("{}  {}\n".format(self.elements.shape[0], self.z.size))
        for i in range(self.z.size):
            fort14.append("{:>10} {:15.10f} {:15.10f} {:15.10f}\n".format(
                        self.node_id[i], self.x[i], self.y[i], self.z[i]))
        for i in range(self.elements.shape[0]):
            fort14.append("{:5d}{:5d} {:5d} {:5d} {:5d}\n".format(
                                    self.element_id[i], 3,
                                    self.elements[i, 0]+1,
                                    self.elements[i, 1]+1,
                                    self.elements[i, 2]+1))
        fort14.append("{:d}".format(len(self.OceanBoundaries))
                      + " = Number of open boundaries\n")
        for i, boundary in enumerate(self.OceanBoundaries):
            Geometry = self.OceanBoundaries[boundary]['Geometry']
            print(Geometry)
        with open(path + '/' + filename, 'w') as f:
            f.write(self.fort14)
            #     _sum += ocean_boundaries[i].size
            # f.write("{:d} = Number of open boundaries\n".format(len(ocean_boundaries)))
            # f.write("{:d} = Total number of open boundary nodes\n".format(_sum))
            # if self.ocean_boundaries is not None:
            #     for i in range(len(self.ocean_boundaries)):
            #         f.write("{:d} = Number of nodes for open boundary {:d}\n".format(ocean_boundaries[i].size, i+1))
            #         for j in range(ocean_boundaries[i].size):
            #             f.write("{:d}\n".format(ocean_boundaries[i][j]+1))
            # remainingBoundaries = 0
            # _sum=0
            # if self.land_boundaries is not None:
            #     remainingBoundaries+=len(self.land_boundaries)
            #     for i in range(len(self.land_boundaries)):
            #         _sum +=  self.land_boundaries[i][0].size

            # if self.inner_boundaries is not None:
            #     remainingBoundaries+=len(self.inner_boundaries)
            #     for i in range(len(self.inner_boundaries)):
            #         _sum +=  self.inner_boundaries[i][0].size

            # if self.inflow_boundaries is not None:
            #     remainingBoundaries+=len(self.inflow_boundaries)
            #     for i in range(len(self.inflow_boundaries)):
            #         _sum +=  self.inflow_boundaries[i][0].size

            # if self.outflow_boundaries is not None:
            #     remainingBoundaries+=len(self.outflow_boundaries)
            #     for i in range(len(self.outflow_boundaries)):
            #         _sum +=  self.outflow_boundaries[i][0].size

            # if self.weir_boundaries is not None:
            #     remainingBoundaries+=len(self.weir_boundaries)
            #     for i in range(len(self.weir_boundaries)):
            #         _sum +=  self.weir_boundaries[i]['front_face'].size
            #         _sum +=  self.weir_boundaries[i]['back_face'].size

            # if self.culvert_boundaries is not None:
            #     remainingBoundaries+=len(self.culvert_boundaries)
            #     for i in range(len(self.culvert_boundaries)):
            #         _sum +=  self.culvert_boundaries[i]['front_face'].size
            #         _sum +=  self.culvert_boundaries[i]['back_face'].size

            # f.write("{:d} = Number of land boundaries\n".format(remainingBoundaries))
            # f.write("{:d} = Total number of land boundary nodes\n".format(_sum))

            # _cnt = 1
            # if self.land_boundaries is not None:
            #     for i in range(len(self.land_boundaries)):
            #         f.write("{:d} {:d} = Number of nodes for land boundary {:d}\n".format(
            #                         self.land_boundaries[i][0].size, self.land_boundaries[i][-1], _cnt))
            #         _cnt+=1
            #         for j in range(self.land_boundaries[i][0].size):
            #             f.write("{:d}\n".format(self.land_boundaries[i][0][j]+1))

            # if self.inner_boundaries is not None:
            #     for i in range(len(self.inner_boundaries)):
            #         f.write("{:d} {:d} = Number of nodes for closed (\"island\")  boundary (land boundary {:d})\n".format(
            #                                 self.inner_boundaries[i][0].size, self.inner_boundaries[i][-1], _cnt))
            #         _cnt+=1
            #         for j in range(self.inner_boundaries[i][0].size):
            #             f.write("{:d}\n".format(self.inner_boundaries[i][0][j]+1))

            # if self.inflow_boundaries is not None:
            #     for i in range(len(self.inflow_boundaries)):
            #         f.write("{:d} {:d} = Number of nodes for inflow boundary (land boundary {:d})\n".format(
            #                                 self.inflow_boundaries[i][0].size, self.inflow_boundaries[i][-1], _cnt))
            #         _cnt+=1
            #         for j in range(self.inflow_boundaries[i][0].size):
            #             f.write("{:d}\n".format(self.inflow_boundaries[i][0][j]+1))

            # if self.outflow_boundaries is not None:
            #     for i in range(len(self.outflow_boundaries)):
            #         f.write("{:d} {:d} = Number of nodes for outflow boundary (land boundary {:d})\n".format(
            #                                 self.outflow_boundaries[i][0].size, self.outflow_boundaries[i][-1], _cnt))
            #         _cnt+=1
            #         for j in range(self.outflow_boundaries[i][0].size):
            #             f.write("{:d} {:.3f} {:.3f}\n".format(self.outflow_boundaries[i][0][j]+1, self.outflow_boundaries[i][1][j],
            #             self.outflow_boundaries[i][2][j]))

            # if self.weir_boundaries is not None:
            #     for i in range(len(self.weir_boundaries)):
            #         f.write("{:d} {:d} = Number of node pairs for weir (land boundary {:d})\n".format(
            #                                 self.weir_boundaries[i]['front_face'].size, self.weir_boundaries[i]['btype'], _cnt))
            #         _cnt+=1
            #         for j in range(self.weir_boundaries[i]['front_face'].size):
            #             f.write("{:d} {:d} {:.3f} {:.3f} {:.3f}\n".format(
            #             self.weir_boundaries[i]['front_face'][j]+1, self.weir_boundaries[i]['back_face'][j]+1,
            #             self.weir_boundaries[i]['barrier_heights'][j], self.weir_boundaries[i]['subcritical_flow_coefficient'][j],
            #             self.weir_boundaries[i]['supercritical_flow_coefficient'][j]))

            # if self.culvert_boundaries is not None:
            #     for i in range(len(self.culvert_boundaries)):
            #         f.write("{:d} {:d} = Number of nodes pairs for culvert boundary (land boundary {:d})\n".format(
            #                                 len(self.culvert_boundaries[i]['front_face']), self.culvert_boundaries[i]['btype'], _cnt))
            #         _cnt+=1
            #         for j in range(self.culvert_boundaries[i]['front_face'].size):
            #             f.write("{:d} {:d} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}\n".format(
            #             self.culvert_boundaries[i]['front_face'][j]+1, self.culvert_boundaries[i]['back_face'][j]+1,
            #             self.culvert_boundaries[i]['barrier_heights'][j],     self.culvert_boundaries[i]['subcritical_flow_coefficient'][j],
            #             self.culvert_boundaries[i]['supercritical_flow_coefficient'][j],     self.culvert_boundaries[i]['cross_barrier_pipe_height'][j],
            #             self.culvert_boundaries[i]['friction_factor'][j],     self.culvert_boundaries[i]['pipe_diameter'][j]))
            pass

    @staticmethod
    def parse_fort14(path):
        fort14 = dict()
        fort14['x'] = list()
        fort14['y'] = list()
        fort14['z'] = list()
        fort14['node_id'] = list()
        fort14['elements'] = list()
        fort14['element_id'] = list()
        # fort14['numberOfConnectedElements'] = list()
        fort14['OceanBoundaries'] = list()
        fort14['LandBoundaries'] = list()
        fort14['InnerBoundaries'] = list()
        fort14['InflowBoundaries'] = list()
        fort14['OutflowBoundaries'] = list()
        fort14['WeirBoundaries'] = list()
        fort14['CulvertBoundaries'] = list()
        with open(path, 'r') as f:
            fort14['description'] = "{}".format(f.readline())
            NE, NP = map(int, f.readline().split())
            _NP = len([])
            while _NP < NP:
                node_id, x, y, z = f.readline().split()
                fort14['node_id'].append(int(node_id)-1)
                fort14['x'].append(float(x))
                fort14['y'].append(float(y))
                fort14['z'].append(-float(z))
                _NP += 1
            _NE = len([])
            while _NE < NE:
                line = f.readline().split()
                fort14['element_id'].append(float(line[0]))
                # fort14['numberOfConnectedElements'].append(int(line[1]))
                if int(line[1]) != 3:
                    raise NotImplementedError('Package only supports '
                                              + 'triangular meshes for the '
                                              + 'time being.')
                fort14['elements'].append([int(x)-1 for x in line[2:]])
                _NE += 1
            NOPE = int(f.readline().split()[0])
            _NOPE = len([])
            f.readline()  # Number of total open ocean nodes. Not used.
            while _NOPE < NOPE:
                fort14['OceanBoundaries'].append({'indexes': list()})
                NETA = int(f.readline().split()[0])
                _NETA = len([])
                while _NETA < NETA:
                    fort14['OceanBoundaries'][_NOPE]['indexes'].append(
                                                int(f.readline().split()[0])-1)
                    _NETA += 1
                _NOPE += 1
            # Assume EOF if NBOU is empty.
            try:
                NBOU = int(f.readline().split()[0])
            except IndexError:
                return fort14
            # For now, let -1 mean a self closing mesh
            # reassigning NBOU to 0 until further implementation is applied.
            if NBOU == -1:
                NBOU = 0
            _NBOU = len([])
            f.readline()
            while _NBOU < NBOU:
                NVELL, IBTYPE = map(int, f.readline().split()[:2])
                _NVELL = 0
                if IBTYPE in [0, 10, 20]:
                    fort14['LandBoundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list()})
                elif IBTYPE in [1, 11, 21]:
                    fort14['InnerBoundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list()})
                elif IBTYPE in [2, 12, 22, 102, 122]:
                    fort14['InflowBoundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list()})
                elif IBTYPE in [3, 13, 23]:
                    fort14['OutflowBoundaries'].append({
                                    'ibtype': IBTYPE,
                                    'indexes': list(),
                                    'barrier_heights': list(),
                                    'supercritical_flow_coefficients': list()})

                elif IBTYPE in [4, 24]:
                    fort14['WeirBoundaries'].append({
                                    'ibtype': IBTYPE,
                                    'front_face_indexes': list(),
                                    'back_face_indexes': list(),
                                    'barrier_heights': list(),
                                    'subcritical_flow_coefficients': list(),
                                    'supercritical_flow_coefficients': list()})
                elif IBTYPE in [5, 25]:
                    fort14['CulvertBoundaries'].append({
                                    'ibtype': IBTYPE,
                                    'front_face_indexes': list(),
                                    'back_face_indexes': list(),
                                    'barrier_heights': list(),
                                    'subcritical_flow_coefficients': list(),
                                    'supercritical_flow_coefficients': list(),
                                    'cross_barrier_pipe_heights': list(),
                                    'friction_factors': list(),
                                    'pipe_diameters': list()})
                else:
                    raise Exception('IBTYPE={} '.format(IBTYPE)
                                    + 'found in fort.14 not recongnized. ')
                while _NVELL < NVELL:
                    line = f.readline().split()
                    if IBTYPE in [0, 10, 20]:
                        fort14['LandBoundaries'][-1]['indexes'] \
                            .append(int(line[0])-1)
                    elif IBTYPE in [1, 11, 21]:
                        fort14['InnerBoundaries'][-1]['indexes'] \
                            .append(int(line[0])-1)
                    elif IBTYPE in [3, 13, 23]:
                        fort14['OutflowBoundaries'][-1]['indexes'] \
                            .append(int(line[0])-1)
                        fort14['OutflowBoundaries'][-1]['external_'
                                                        + 'barrier_heights'] \
                            .append(float(line[1]))
                        fort14['OutflowBoundaries'][-1]['supercritical_flow'
                                                        + '_coefficients'] \
                            .append(float(line[2]))
                    elif IBTYPE in [2, 12, 22, 102, 122]:
                        fort14['iflowBoundaries'][-1]['indexes'] \
                            .append(int(line[0])-1)
                    elif IBTYPE in [4, 24]:
                        fort14['WeirBoundaries'][-1]['front_face_indexes'] \
                            .append(int(line[0])-1)
                        fort14['WeirBoundaries'][-1]['back_face_indexes'] \
                            .append(int(line[1])-1)
                        fort14['WeirBoundaries'][-1]['barrier_heights'] \
                            .append(float(line[2]))
                        fort14['WeirBoundaries'][-1]['subcritical_'
                                                     + 'flow_coefficients'] \
                            .append(float(line[3]))
                        fort14['WeirBoundaries'][-1]['supercritical_'
                                                     + 'flow_coefficients'] \
                            .append(float(line[4]))
                    elif IBTYPE in [5, 25]:
                        fort14['CulvertBoundaries'][-1]['front_face_indexes'] \
                            .append(int(line[0])-1)
                        fort14['CulvertBoundaries'][-1]['back_face_indexes'] \
                            .append(int(line[1])-1)
                        fort14['CulvertBoundaries'][-1]['barrier_heights'] \
                            .append(float(line[2]))
                        fort14['CulvertBoundaries'][-1]['subcritical_flow_'
                                                        + 'coefficients'] \
                            .append(float(line[3]))
                        fort14['CulvertBoundaries'][-1]['supercritical_flow_'
                                                        + 'coefficients'] \
                            .append(float(line[4]))
                        fort14['CulvertBoundaries'][-1]['friction'
                                                        + 'Factor'] \
                            .append(float(line[5]))
                        fort14['CulvertBoundaries'][-1]['pipe_diameter'] \
                            .append(float(line[6]))
                    else:
                        Exception("Duck-typing error. "
                                  + "This exception should be unreachable.")
                    _NVELL += 1
                _NBOU += 1
        return fort14

    @property
    def node_id(self):
        return self._node_id

    @property
    def element_id(self):
        return self._element_id

    @property
    def description(self):
        return self._description

    @property
    def fort14(self):
        return self.__fort14

    @property
    def _node_id(self):
        return self.__node_id

    @property
    def _element_id(self):
        return self.__element_id

    @description.setter
    def description(self, description):
        self._description = description

    @_node_id.setter
    def _node_id(self, node_id):
        if node_id is None:
            self.__node_id = np.arange(1, self.z.shape[0]+1)
        else:
            node_id = np.asarray(node_id).flatten()
            assert node_id.shape[0] == self.z.shape[0]
            self.__node_id = node_id

    @_element_id.setter
    def _element_id(self, element_id):
        if element_id is None:
            self.__element_id = np.arange(1, self.elements.shape[0]+1)
        else:
            element_id = np.asarray(element_id)
            assert element_id.shape[0] == self.elements.shape[0]
            self.__element_id = element_id


class AdcircMeshTestCase(unittest.TestCase):

    def setUp(self):
        self.HSOFS = AdcircMesh.from_fort14(os.getenv("HSOFS"), 4326, 'LMSL')

    def test_set_fort13(self):
        self.HSOFS.set_fort13(os.getenv("HSOFS_FORT13"))

    def _test_make_plot(self):
        self.HSOFS.make_plot(show=False)
