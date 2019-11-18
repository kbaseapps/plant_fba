# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from plant_fba.plant_fbaImpl import plant_fba
from plant_fba.plant_fbaServer import MethodContext
from plant_fba.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class plant_fbaTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('plant_fba'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'plant_fba',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = plant_fba(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_plant_fba_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_integrate_abundances_with_metabolism(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})

        ret = self.getImpl().reconstruct_plant_metabolism(self.getContext(), {'input_ws': 'plant_fba_testing',
                                                                              'input_genome': 'PlantSEED_Arabidopsis',

                                                                              #'output_ws': 'plant_fba_testing',
                                                                              'output_fbamodel': 'PlantSEED_Arabidopsis_FBAModel',

                                                                              #'template_ws': 'NewKBaseModelTemplates',
                                                                              #'template': 'PlantModelTemplate'
                                                                              })

        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
