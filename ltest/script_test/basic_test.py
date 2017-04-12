# Test script for genome_util package - it should be launched from
# the root of the genome_util module, ideally just with 'make test', as
# it looks for a hardcoded relative path to find the 'test.cfg' file
import unittest
import json
import ConfigParser
from os import environ
from string import Template
from pprint import pprint

from subprocess import call

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint
import sys


from biokbase.CoExpression.authclient import KBaseAuth as _KBaseAuth

from biokbase.auth import Token
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3
from biokbase.CoExpression.authclient import KBaseAuth as _KBaseAuth
from biokbase.workspace.client import Workspace as Workspace

INPUT_META_DATA_DIR = 'ltest/script_test/input_meta_data'
INPUT_DATA_DIR = 'ltest/script_test/input_data'

# Before all the tests, read the config file and get a user token and
# save it to a file used by the main service script
class TestCoExpressionMethodsSetup(unittest.TestCase):

    @classmethod
    def setUp(cls):

      token = environ.get('KB_AUTH_TOKEN', None)

      if token is None:
          sys.stderr.write("Error: Unable to run tests without authentication token!\n")
          sys.exit(1)

      token_file = open('ltest/script_test/token.txt', 'w')
      token_file.write(token)

      config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
      cls.cfg = {}
      config = ConfigParser()
      config.read(config_file)
      for nameval in config.items('CoExpression'):
          cls.cfg[nameval[0]] = nameval[1]
      auth_service_url = cls.cfg.get('auth-service-url',
                                     "https://kbase.us/services/authorization/Sessions/Login")
      ws_url = cls.cfg['ws_url']
      auth_service_url_allow_insecure = cls.cfg['auth-service-url-allow-insecure']
      auth_client = _KBaseAuth(auth_service_url)
      user_id = auth_client.get_user(token)

      ws = Workspace(url=ws_url, token=token, auth_svc=auth_service_url,
                             trust_all_ssl_certificates=auth_service_url_allow_insecure)

      # update input data in reverse order of references
      ordered_file_list = [INPUT_META_DATA_DIR+'/test_diff_p_distribution_input_ref2.json',
                      INPUT_META_DATA_DIR+'/test_diff_p_distribution_input_ref1.json',
                      INPUT_META_DATA_DIR+'/test_diff_p_distribution_input.json',
                      INPUT_META_DATA_DIR+'/test_view_heatmap_input_ref1.json',
                      INPUT_META_DATA_DIR+'/test_view_heatmap_input.json',
                      INPUT_META_DATA_DIR+'/test_coex_clust_input.json',
                      INPUT_META_DATA_DIR+'/test_filter_genes_input.json']

      for filename in ordered_file_list:
          with open(filename, 'r') as infile:
            input_meta_data = json.load(infile)

          # create workspace that is local to the user if it does not exist
          workspace_name_t = Template(str(input_meta_data['params'][0]['workspace_name']))
          workspace_name = workspace_name_t.substitute(user_id=user_id)
          print('workspace_name: ' + workspace_name)

          try:
              ws_info = ws.get_workspace_info({'workspace': workspace_name})
              print("workspace already exists: " + str(ws_info))
          except:
              ws_info = ws.create_workspace(
                {'workspace': workspace_name, 'description': 'Workspace for ' + str(input_meta_data['method'])})
              print("Created new workspace: " + str(ws_info))

          print('reading input file: '+filename)
          object_name = str(input_meta_data['params'][0]['object_name'])
          print('object_name: '+object_name)

          input_data_filename = INPUT_DATA_DIR + '/' + object_name + '.json'
          print('input data filename: ' + input_data_filename)

          with open(input_data_filename, 'r') as infile:
            input_data = json.load(infile)

          # update workspace name in input data
          input_data_str = json.dumps(input_data)
          input_data_t = Template(input_data_str)
          input_data_str = input_data_t.substitute(workspace_name=workspace_name)
          input_data = json.loads(input_data_str)

          print('type: '+input_data[0]['info'][2])

          #upload data (no effect if data already exists in workspace)
          print('uploading input data to workspace')
          ws.save_objects(
                  {'workspace': workspace_name, 'objects': [{'type': input_data[0]['info'][2],
                                                                  'data': input_data[0]['data'],
                                                                  'name': object_name}]})
      print('ws objects: ' + str(ws.list_objects({'workspaces': [workspace_name]})))

  # Define all our other test cases here
class TestCoExpressionMethods(TestCoExpressionMethodsSetup):

   def test_diff_p_distribution(self):
          print("\n\n----------- test diff_p_distribution  ----------")

          out =call(["run_CoExpression.sh",
          "ltest/script_test/input_meta_data/test_diff_p_distribution_input.json",
          "ltest/script_test/test_diff_p_distribution_output.json",
          "ltest/script_test/token.txt"])

          # print error code of Implementation
          print(out);

          with open('ltest/script_test/test_diff_p_distribution_output.json') as o:
                  output =json.load(o)
          pprint(output)

   def test_view_heatmap(self):
          print("\n\n----------- test view_heatmap  ----------")

          out =call(["run_CoExpression.sh",
          "ltest/script_test/input_meta_data/test_view_heatmap_input.json",
          "ltest/script_test/test_view_heatmap_output.json",
          "ltest/script_test/token.txt"])

          # print error code of Implementation
          print(out);

          with open('ltest/script_test/test_view_heatmap_output.json') as o:
                  output =json.load(o)
          pprint(output)

   def test_filter_genes(self):
          print("\n\n----------- test filter genes ----------")

          out =call(["run_CoExpression.sh",
          "ltest/script_test/input_meta_data/test_filter_genes_input.json",
          "ltest/script_test/test_filter_genes_output.json",
          "ltest/script_test/token.txt"])

          # print error code of Implementation
          print(out);


          with open('ltest/script_test/test_filter_genes_output.json') as o:
                  output =json.load(o)
          pprint(output)

   def test_coex_cluster(self):
          print("\n\n----------- test constcoex_net_clust ----------")

          out =call(["run_CoExpression.sh",
          "ltest/script_test/input_meta_data/test_coex_clust_input.json",
          "ltest/script_test/test_coex_clust_output.json",
          "ltest/script_test/token.txt"])

          # print error code of Implementation
          print(out);


          with open('ltest/script_test/test_coex_clust_output.json') as o:
                  output =json.load(o)
          pprint(output)

#start the tests if run as a script
if __name__ == '__main__':
    unittest.main()
