#BEGIN_HEADER
import os
import sys
import traceback
import argparse
import json
import logging
import time
from  pprint import pprint
import string
import subprocess
from os import environ
from ConfigParser import ConfigParser
import re
from collections import OrderedDict
import uuid
from string import Template

# 3rd party imports
import requests
import pandas as pd
import numpy as np 

# KBase imports
import biokbase.workspace.client
from biokbase.CoExpression.authclient import KBaseAuth as _KBaseAuth
import biokbase.Transform.script_utils as script_utils 

def error_report(err_msg, expr, workspace_service_url, workspace_name, provenance, ws):

    ## Create report object:
    report = "Failed : {0}".format(err_msg)
    reportObj = {
                    'objects_created':[],
                    'text_message':report
                }

    # generate a unique name for the Method report
    reportName = 'ErrorReport_'+str(hex(uuid.getnode()))
    report_info = ws.save_objects({
                                    'workspace':workspace_name,
                                    'objects':[
                                    {
                                    'type':'KBaseReport.Report',
                                    'data':reportObj,
                                    'name':reportName,
                                    'meta':{},
                                    'hidden':1, 
                                    'provenance':provenance
                                    }
                                    ]
                                    })[0]

    return { "report_name" : reportName,"report_ref" : "{0}/{1}/{2}".format(report_info[6],report_info[0],report_info[4]) }


#END_HEADER


class CoExpression:
    '''
    Module Name:
    CoExpression

    Module Description:
    Co-Expression Service APIs 

 This module provides services in support of the coexpression network. 
 The modules supports retrieval of the following information:
 1. Identify differentially expressed genes
 2. WGCNA clustering
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "1.0.4"
    GIT_URL = "https://github.com/sjyoo/coexpression"
    GIT_COMMIT_HASH = "cd2d4e4351c6e80cf4b642f72f4fd1a5c0576617"
    
    #BEGIN_CLASS_HEADER
    FVE_2_TSV = 'trns_transform_KBaseFeatureValues_ExpressionMatrix_to_TSV'
    TSV_2_FVE = 'trns_transform_TSV_Exspression_to_KBaseFeatureValues_ExpressionMatrix'
    RAWEXPR_DIR = 'raw_dir'
    FLTRD_DIR = 'fltr_dir'
    CLSTR_DIR = 'clstr_dir'
    FINAL_DIR = 'final_dir'
    EXPRESS_FN = 'expression.tsv'
    SAMPLE_FN = 'sample.tsv'
    COEX_FILTER = 'coex_filter'
    COEX_CLUSTER = 'coex_cluster2'
    FLTRD_FN = 'filtered.tsv'
    CLSTR_FN = 'clusters.tsv'
    CSTAT_FN = 'cluster_stat.tsv'
    FINAL_FN = 'filtered.json'
    PVFDT_FN = 'pv_distribution.json'
    GENELST_FN = 'selected.tsv'
    __WS_URL = 'https://ci.kbase.us/services/ws'
    __HS_URL = 'https://ci.kbase.us/services/handle_service'
    __SHOCK_URL = 'https://ci.kbase.us/services/shock-api'
    __PUBLIC_SHOCK_NODE = 'true'
    logger = None


    def _dumpExp2File(self, oexpr, odir, ofn):
        try:
            os.makedirs(odir)
        except:
            pass

        try:
            df = pd.DataFrame(oexpr['data']['values'], index=oexpr['data']['row_ids'], columns=oexpr['data']['col_ids'])
            mask = pd.Series(df.index == 'NA').add(pd.Series(df.index == '')).values == 0  # remove 'NA' or '' (missing gene name)
            df = df.iloc[mask,]
            df.to_csv(path_or_buf = odir + "/" + ofn, sep='\t', na_rep = 'NA' )
        except:
            self.logger.error("Failed to dump expression object into tsv file:" + traceback.format_exc());
            raise
        
    def _subselectExp(self, oexpr, gl):
        exp_idx = [oexpr['data']['row_ids'].index(x) for x in gl]
        oexpr['data']['row_ids'] = [oexpr['data']['row_ids'][x] for x in exp_idx]
        oexpr['data']['values'] = [oexpr['data']['values'][x] for x in exp_idx]
        return oexpr

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #pprint(config)
        if 'auth-service-url' in config:
              self.__AUTH_SERVICE_URL = config['auth-service-url']
        if 'ws_url' in config:
              self.__WS_URL = config['ws_url']
        if 'shock_url' in config:
              self.__SHOCK_URL = config['shock_url']
        if 'hs_url' in config:
          self.__HS_URL = config['hs_url']
        if 'cluster_dir' in config:
              self.__CLSTR_DIR = config['cluster_dir']
        if 'final_dir' in config:
              self.__FINAL_DIR = config['final_dir']
        if 'coex_filter' in config:
              self.__COEX_FILTER = config['coex_filter']
        if 'coex_cluster' in config:
              self.__COEX_CLUSTER = config['coex_cluster']
        if 'force_shock_node_2b_public' in config: # expect 'true' or 'false' string
          self.__PUBLIC_SHOCK_NODE = config['force_shock_node_2b_public']
    
        # logging
        self.logger = logging.getLogger('CoExpression')
        if 'log_level' in config:
              self.logger.setLevel(config['log_level'])
        else:
              self.logger.setLevel(logging.INFO)
        streamHandler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        streamHandler.setFormatter(formatter)
        self.logger.addHandler(streamHandler)
        self.logger.info("Logger was set")
        #END_CONSTRUCTOR
        pass
    

    def diff_p_distribution(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN diff_p_distribution
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FLTRD_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
        token = ctx['token']
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token

        param = args

        auth_client = _KBaseAuth(self.__AUTH_SERVICE_URL)
        user_id = auth_client.get_user(token)
        workspace_name_t = Template(param['workspace_name'])
        workspace_name = workspace_name_t.substitute(user_id=user_id)
 
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        expr = ws.get_objects([{'workspace': workspace_name, 'name' : param['object_name']}])[0]['data']
 
 
        self._dumpExp2File(expr, self.RAWEXPR_DIR, self.EXPRESS_FN)
 
        self.logger.info("Identifying differentially expressed genes")
 
        ## Prepare sample file
        # detect num of columns
        ncol = len(expr['data']['col_ids'])
        
        # force to use ANOVA if the number of sample is two
        if(ncol == 3): param['method'] = 'anova'
 
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN), 'wt') as s:
          s.write("0")
          for j in range(1,ncol):
            s.write("\t{0}".format(j))
          s.write("\n")
 
 
        ## Run coex_filter
        cmd_coex_filter = [self.COEX_FILTER, '-i', "{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), '-o', "{0}/{1}".format(self.FLTRD_DIR, self.FLTRD_FN),
           '-m', param['method'], '-n', '10', '-s', "{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN),
           '-x', "{0}/{1}".format(self.RAWEXPR_DIR, self.GENELST_FN), '-t', 'y', '-j', self.PVFDT_FN]
        if 'num_features' in param:
          cmd_coex_filter.append("-n")
          cmd_coex_filter.append(str(param['num_features']))
 
        if 'p_value' in param:
          cmd_coex_filter.append("-p")
          cmd_coex_filter.append(str(param['p_value']))
 
 
        tool_process = subprocess.Popen(cmd_coex_filter, stderr=subprocess.PIPE)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        ## loading pvalue distribution FDT
        pvfdt = {'row_labels' :[], 'column_labels' : [], "data" : [[]]};
        pvfdt = OrderedDict(pvfdt)
        with open(self.PVFDT_FN, 'r') as myfile:
           pvfdt = json.load(myfile)
        data_obj_name = "{0}.fdt".format(param['out_figure_object_name'])
        pvfdt['id'] = data_obj_name
 
 
        fig_properties = {"xlabel" : "-log2(p-value)", "ylabel" : "Number of features", "xlog_mode" : "-log2", "ylog_mode" : "none", "title" : "Histogram of P-values", "plot_type" : "histogram"}
        sstatus = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'MAK.FloatDataTable',
                                                                              'data' : pvfdt,
                                                                              'name' : data_obj_name}]})

        data_ref = "{0}/{1}/{2}".format(sstatus[0][6], sstatus[0][0], sstatus[0][4])
        fig_properties['data_ref'] = data_ref

        sstatus = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'CoExpression.FigureProperties',
                                                                              'data' : fig_properties,
                                                                              'name' : (param['out_figure_object_name'])}]})
        result = fig_properties
        #END diff_p_distribution

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method diff_p_distribution return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def filter_genes(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN filter_genes
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FLTRD_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
        token = ctx['token']
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token

        param = args

        auth_client = _KBaseAuth(self.__AUTH_SERVICE_URL)
        user_id = auth_client.get_user(token)
        workspace_name_t = Template(param['workspace_name'])
        workspace_name = workspace_name_t.substitute(user_id=user_id)
 
        provenance = [{}]
        if 'provenance' in ctx:
                provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[workspace_name+'/'+param['object_name']]
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        expr = ws.get_objects([{'workspace': workspace_name, 'name' : param['object_name']}])[0]['data']
 
        self._dumpExp2File(expr, self.RAWEXPR_DIR, self.EXPRESS_FN)
 
        self.logger.info("Identifying differentially expressed genes")
 
        ## Prepare sample file
        # detect num of columns
        ncol = len(expr['data']['col_ids'])
        
        # force to use ANOVA if the number of sample is two
        if(ncol == 3): param['method'] = 'anova'
 
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN), 'wt') as s:
          s.write("0")
          for j in range(1,ncol):
            s.write("\t{0}".format(j))
          s.write("\n")
 
 
        ## Run coex_filter
        cmd_coex_filter = [self.COEX_FILTER, '-i', "{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), '-o', "{0}/{1}".format(self.FLTRD_DIR, self.FLTRD_FN),
                           '-m', param['method'], '-s', "{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN),
                           '-x', "{0}/{1}".format(self.RAWEXPR_DIR, self.GENELST_FN), '-t', 'y']
        if 'num_features' in param:
          cmd_coex_filter.append("-n")
          cmd_coex_filter.append(str(param['num_features']))
 
        if 'p_value' in param:
          cmd_coex_filter.append("-p")
          cmd_coex_filter.append(str(param['p_value']))
 
        if 'p_value' not in param and 'num_features' not in param:
          self.logger.error("One of p_value or num_features must be defined");
          return error_report("One of p_value or num_features must be defined", expr,self.__WS_URL, workspace_name, provenance, ws)
          #sys.exit(2) #TODO: No error handling in narrative so we do graceful termination
 
        #if 'p_value' in param and 'num_features' in param:
        #  self.logger.error("Both of p_value and num_features cannot be defined together");
        #  sys.exit(3)
 
        tool_process = subprocess.Popen(cmd_coex_filter, stderr=subprocess.PIPE)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            self.logger.info(stderr)
 
        ## checking genelist
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.GENELST_FN),'r') as glh:
          gl = glh.readlines()
        gl = [x.strip('\n') for x in gl]
 
        if(len(gl) < 1) :
          self.logger.error("No genes are selected")
          return error_report("Increase p_value or specify num_features", expr,self.__WS_URL, workspace_name, provenance, ws)
          #sys.exit(4)
 
        ## Upload FVE
        if 'description' not in expr: 
            expr['description'] = "Filtered Expression Matrix"
        expr['description'] += " : Filtered by '{1}' method ".format(expr['description'], param['method'])
 
        expr = self._subselectExp(expr, gl)
 
        ex_info = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'KBaseFeatureValues.ExpressionMatrix',
                                                                              'data' : expr,
                                                                              'name' : (param['out_expr_object_name'])}]})[0]
 
        ## Upload FeatureSet
        fs ={'elements': {}}
        fs['description'] = "FeatureSet identified by filtering method '{0}' ".format(param['method'])
 
        fs['description'] += "from {0}/{1}".format(workspace_name, param['object_name'])
 
        for g in gl:
          if 'genome_ref' in expr:
            fs['elements'][g] = [expr['genome_ref']]
          else:
            fs['elements'][g] = []
 
        fs_info = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'KBaseCollections.FeatureSet',
                                                                              'data' : fs,
                                                                              'name' : (param['out_fs_object_name'])}]})[0]

        ## Create report object:
	report = "Filtering expression matrix using {0} on {1}".format(param['method'],param['object_name'])
        reportObj = {
                        'objects_created':[{
                                'ref':"{0}/{1}/{2}".format(fs_info[6], fs_info[0], fs_info[4]),
                                'description':'Filtered FeatureSet' },
                             {
                                'ref':"{0}/{1}/{2}".format(ex_info[6], ex_info[0], ex_info[4]),
                                'description':'Filetered ExpressionMatrix' 
                             }],
                        'text_message':report
                    }

        # generate a unique name for the Method report
        reportName = 'FilterExpression_'+str(hex(uuid.getnode()))
        report_info = ws.save_objects({
                                        'id':ex_info[6],
                                        'objects':[
                                        {
                                        'type':'KBaseReport.Report',
                                        'data':reportObj,
                                        'name':reportName,
                                        'meta':{},
                                        'hidden':1, 
                                        'provenance':provenance
                                        }
                                        ]
                                        })[0]

        result = { "report_name" : reportName,"report_ref" : "{0}/{1}/{2}".format(report_info[6],report_info[0],report_info[4]) }



        #result = {'workspace_name' : workspace_name, 'out_expr_object_name' : param['out_expr_object_name'], 'out_fs_object_name' : param['out_fs_object_name']}
        #END filter_genes

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method filter_genes return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def const_coex_net_clust(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN const_coex_net_clust
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.CLSTR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Starting conversion of KBaseFeatureValues.ExpressionMatrix to TSV")
        token = ctx['token']

        param = args

        auth_client = _KBaseAuth(self.__AUTH_SERVICE_URL)
        user_id = auth_client.get_user(token)
        workspace_name_t = Template(param['workspace_name'])
        workspace_name = workspace_name_t.substitute(user_id=user_id)

        provenance = [{}]
        if 'provenance' in ctx:
                provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[workspace_name+'/'+param['object_name']]
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        expr = ws.get_objects([{'workspace': workspace_name, 'name' : param['object_name']}])[0]['data']
 
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token
        self._dumpExp2File(expr, self.RAWEXPR_DIR, self.EXPRESS_FN)
 
        self.logger.info("Identifying differentially expressed genes")
 
        ## Prepare sample file
        # detect num of columns
        ncol = len(expr['data']['col_ids'])
        
        # grouping information 
        with open("{0}/{1}".format(self.RAWEXPR_DIR, self.SAMPLE_FN), 'wt') as s:
          s.write("0")
          for j in range(1,ncol):
            s.write("\t{0}".format(j))
          s.write("\n")
 
 
        ## Run coex_cluster
        cmd_coex_cluster = [self.COEX_CLUSTER, '-t', 'y',
                           '-i', "{0}/{1}".format(self.RAWEXPR_DIR, self.EXPRESS_FN), 
                           '-o', "{0}/{1}".format(self.CLSTR_DIR, self.CLSTR_FN), '-m', "{0}/{1}".format(self.CLSTR_DIR, self.CSTAT_FN) ]
 
        for p in ['net_method', 'minRsq', 'maxmediank', 'maxpower', 'clust_method', 'minModuleSize', 'detectCutHeight']:
           if p in param:
             cmd_coex_cluster.append("--{0}".format(p))
             cmd_coex_cluster.append(str(param[p]))
  
 
        #sys.exit(2) #TODO: No error handling in narrative so we do graceful termination
 
        #if 'p_value' in param and 'num_features' in param:
        #  self.logger.error("Both of p_value and num_features cannot be defined together");
        #  sys.exit(3)
 
        tool_process = subprocess.Popen(cmd_coex_cluster, stderr=subprocess.PIPE)
        stdout, stderr = tool_process.communicate()
        
        if stdout is not None and len(stdout) > 0:
            self.logger.info(stdout)
 
        if stderr is not None and len(stderr) > 0:
            if re.search(r'^There were \d+ warnings \(use warnings\(\) to see them\)', stderr):
              self.logger.info(stderr)
            else:
              self.logger.error(stderr)
              raise Exception(stderr)
 
        
        # build index for gene list
        pos_index ={expr['data']['row_ids'][i]: i for i in range(0, len(expr['data']['row_ids']))}
 
 
        # parse clustering results
        cid2genelist = {}
        cid2stat = {}
        with open("{0}/{1}".format(self.CLSTR_DIR, self.CSTAT_FN),'r') as glh:
            glh.readline() # skip header
            for line in glh:
                cluster, mcor, msec = line.rstrip().replace('"','').split("\t")
                cid2stat[cluster]= [mcor, msec]
        with open("{0}/{1}".format(self.CLSTR_DIR, self.CLSTR_FN),'r') as glh:
            glh.readline() # skip header
            for line in glh:
                gene, cluster = line.rstrip().replace('"','').split("\t")
                if cluster not in cid2genelist:
                    cid2genelist[cluster] = []
                cid2genelist[cluster].append(gene)
 
        if(len(cid2genelist) < 1) :
          self.logger.error("Clustering failed")
          return error_report("Error: No cluster output", expr,self.__WS_URL, workspace_name, provenance, ws)
          #sys.exit(4)
 
        self.logger.info("Uploading the results onto WS")
        feature_clusters = []
        for cluster in cid2genelist:
            feature_clusters.append( {"meancor": float(cid2stat[cluster][0]), "msec": float(cid2stat[cluster][0]), "id_to_pos" : { gene : pos_index[gene] for gene in cid2genelist[cluster]}})

        ## Upload Clusters
        feature_clusters ={"original_data": "{0}/{1}".format(workspace_name,param['object_name']),
                           "feature_clusters": feature_clusters}
 
        cl_info = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'KBaseFeatureValues.FeatureClusters',
                                                                          'data' : feature_clusters,
                                                                          'name' : (param['out_object_name'])}]})[0]
        ## Create report object:
	report = "Clustering expression matrix using WGCNA on {0}".format(param['object_name'])
        reportObj = {
                        'objects_created':[                             {
                                'ref':"{0}/{1}/{2}".format(cl_info[6], cl_info[0], cl_info[4]),
                                'description':'WGCNA FeatureClusters' 
                             }],
                        'text_message':report
                    }

        # generate a unique name for the Method report
        reportName = 'WGCNA_Clusters_'+str(hex(uuid.getnode()))
        report_info = ws.save_objects({
                                        'id':cl_info[6],
                                        'objects':[
                                        {
                                        'type':'KBaseReport.Report',
                                        'data':reportObj,
                                        'name':reportName,
                                        'meta':{},
                                        'hidden':1, 
                                        'provenance':provenance
                                        }
                                        ]
                                        })[0]

        #result = { "report_name" : reportName,"report_ref" : "{0}/{1}/{2}".format(report_info[6],report_info[0],report_info[4]) }
        #result = {'workspace_name' : workspace_name, 'out_object_name' : param['out_object_name']}
        result = {'workspace' : workspace_name, 'output' : param['out_object_name']}
        #END const_coex_net_clust

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method const_coex_net_clust return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def view_heatmap(self, ctx, args):
        # ctx is the context object
        # return variables are: result
        #BEGIN view_heatmap
        try:
            os.makedirs(self.RAWEXPR_DIR)
        except:
            pass
        try:
            os.makedirs(self.FLTRD_DIR)
        except:
            pass
        try:
            os.makedirs(self.FINAL_DIR)
        except:
            pass
 
        if self.logger is None:
            self.logger = script_utils.stderrlogger(__file__)
        
        result = {}
        self.logger.info("Loading data")
        token = ctx['token']
 
        eenv = os.environ.copy()
        eenv['KB_AUTH_TOKEN'] = token

        param = args

        auth_client = _KBaseAuth(self.__AUTH_SERVICE_URL)
        user_id = auth_client.get_user(token)
        workspace_name_t = Template(param['workspace_name'])
        workspace_name = workspace_name_t.substitute(user_id=user_id)
 
 
        from biokbase.workspace.client import Workspace
        ws = Workspace(url=self.__WS_URL, token=token)
        fc = ws.get_objects([{'workspace': workspace_name, 'name' : param['object_name']}])[0]['data']
        if 'original_data' not in fc:
            raise Exception("FeatureCluster object does not have information for the original ExpressionMatrix")
        oexpr = ws.get_objects([{ 'ref' : fc['original_data']}])[0]

        df2 = pd.DataFrame(oexpr['data']['data']['values'], index=oexpr['data']['data']['row_ids'], columns=oexpr['data']['data']['col_ids'])

        # L2 normalization
        df3 = df2.div(df2.pow(2).sum(axis=1).pow(0.5), axis=0)

        # type - ? level, ratio, log-ratio  <---> "untransformed"
        # scale - ? probably: raw, ln, log2, log10
        self.logger.info("Expression matrix type: {0}, scale: {1}".format(oexpr['data']['type'],oexpr['data']['scale'] ))
        # do default behavior
        factor = 0.125
        fc_df = df2 + df2[df2 !=0].abs().min().min() * factor
        if param['control_condition']  in fc_df.columns:
            fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[param['control_condition']]], axis=0)).apply(np.log2)
        else:
            fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[0]], axis=0)).apply(np.log2)
        # now fc_df will be reset
        if oexpr['data']['type'] == 'level' or oexpr['data']['type'] == 'untransformed': # need to compute fold changes
            if 'scale' not in oexpr['data'] or oexpr['data']['scale'] == 'raw' or oexpr['data']['scale'] == "1.0":
              factor = 0.125
              fc_df = df2 + df2[df2 !=0].abs().min().min() * factor
              if param['control_condition']  in fc_df.columns:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[param['control_condition']]], axis=0)).apply(np.log2)
              else:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[0]], axis=0)).apply(np.log2)
            else:
              fc_df = df2
              if param['control_condition']  in fc_df.columns:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[param['control_condition']]], axis=0))
              else:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[0]], axis=0))
              if oexpr['data']['scale'] == "log10":
                  fc_df = fc_df/np.log10(2)
              elif oexpr['data']['scale'] == "ln":
                  fc_df = fc_df/np.log(2)
              else:
                  pass
        elif oexpr['data']['type'] == 'ratio':
            fc_df = df2.apply(np.log2)
        elif oexpr['data']['type'] == 'log-ratio':
            fc_df = df2
            if oexpr['data']['scale'] == "log10":
                fc_df = fc_df/np.log10(2)
            elif oexpr['data']['scale'] == "ln":
                fc_df = fc_df/np.log(2)
            else:
                pass

        else: # do the same thing with simple level or untransformed
            if 'scale' not in oexpr['data'] or oexpr['data']['scale'] == 'raw' or oexpr['data']['scale'] == "1.0":
              factor = 0.125
              fc_df = df2 + df2[df2 !=0].abs().min().min() * factor
              if param['control_condition']  in fc_df.columns:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[param['control_condition']]], axis=0)).apply(np.log2)
              else:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[0]], axis=0)).apply(np.log2)
            else:
              fc_df = df2
              if param['control_condition']  in fc_df.columns:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[param['control_condition']]], axis=0))
              else:
                  fc_df = (fc_df.div(fc_df.loc[:,fc_df.columns[0]], axis=0))
              if oexpr['data']['scale'] == "log10":
                  fc_df = fc_df/np.log10(2)
              elif oexpr['data']['scale'] == "ln":
                  fc_df = fc_df/np.log(2)
              else:
                  pass
       
        self.logger.info("Compute cluster statistics")

        cl = {}
        afs = [];
        cid = 1;

        c_stat = pd.DataFrame()
        for cluster in fc['feature_clusters']:
         
          try: 
            fs  = cluster['id_to_pos'].keys()
          except:
            continue # couldn't find feature_set

          fsn = "Cluster_{0}".format(cid)
          cid +=1
          c_stat.loc[fsn,'size'] = len(fs)
          if 'meancor' in cluster:
              c_stat.loc[fsn,'mcor'] = cluster['meancor']
          else:
            pass
            # TODO: Add mean cor calculation later
            #raise Exception("Mean correlation is not included in FeatureCluster object") # now it is NaN

          if 'quantile' in param:
              # enforcing quantile to be in [0 .. 1] rnage
              qt = float(param['quantile'])
              if qt > 1.0: qt = 1.0
              if qt < 0.0: qt = 0.0
              c_stat.loc[fsn,'stdstat'] = fc_df.loc[fs,].std(axis=1).quantile(qt)
          else:
              c_stat.loc[fsn,'stdstat'] = fc_df.loc[fs,].std(axis=1).quantile(0.75)
         

          c1 = df3.loc[fs,].sum(axis=0)
          if df3.loc[fs,].shape[0] < 1: # empty
            continue
          cl[fsn] = fs
          #afs.extend(fs)

          #c1 = df3.loc[fs,].sum(axis=0)
          #c1 = c1 / np.sqrt(c1.pow(2).sum())
          #if(len(cl.keys()) == 1):
          #  centroids = c1.to_frame(fsn).T
          #else:
          #  centroids.loc[fsn] = c1

        # now we have centroids and statistics
        # let's subselect clusters
        min_features = 200
        if 'min_features' in param :
          min_features = param['min_features']
        
        c_stat.loc[:,'nmcor'] = c_stat.loc[:,'mcor'] / c_stat.loc[:,'mcor'].max()
        c_stat.loc[:,'nstdstat'] = c_stat.loc[:,'stdstat'] / c_stat.loc[:,'stdstat'].max()
        
        if 'use_norm_weight' in param and param['use_norm_weight'] != 0:
            if 'quantile_weight' in param:
                c_stat.loc[:,'weight'] = c_stat.loc[:,'nmcor'] + float(param['quantile_weight']) * c_stat.loc[:,'nstdstat']
            else:
                c_stat.loc[:,'weight'] = c_stat.loc[:,'nmcor'] + 1.0                             * c_stat.loc[:,'nstdstat']
        else:
            if 'quantile_weight' in param:
                c_stat.loc[:,'weight'] = c_stat.loc[:,'mcor'] + float(param['quantile_weight']) * c_stat.loc[:,'stdstat']
            else:
                c_stat.loc[:,'weight'] = c_stat.loc[:,'mcor'] + 0.1                             * c_stat.loc[:,'stdstat']

        c_stat.sort_values('weight', inplace=True, ascending=False)

        pprint(c_stat)

        centroids = pd.DataFrame()
        for i in range(c_stat.shape[0]):
            fsn = c_stat.index[i]
            fs = cl[fsn]
            if i != 0 and len(afs) + len(fs) > min_features :
                break;
           
            afs.extend(fs)

            c1 = df3.loc[fs,].sum(axis=0)
            c1 = c1 / np.sqrt(c1.pow(2).sum())
            if(centroids.shape[0] < 1):
              centroids = c1.to_frame(fsn).T
            else:
              centroids.loc[fsn] = c1
           
        pprint(centroids)
        
        if len(cl.keys()) == 0:
            raise Exception("No feature ids were mapped to dataset or no clusters were selected")
        
        # dataset centroid
        dc = df3.loc[afs,].sum(axis=0)
        dc = dc / np.sqrt(dc.pow(2).sum())
    
        
        self.logger.info("Ordering Centroids and Data")
        # the most far away cluster centroid from dataset centroid
        fc = (centroids * dc).sum(axis=1).idxmin()
        # the most far away centroid centroid from fc
        ffc = (centroids * centroids.loc[fc,]).sum(axis=1).idxmin()
        
        # major direction to order on unit ball space
        md = centroids.loc[ffc,] - centroids.loc[fc,]
        
        # unnormalized component of projection to the major direction (ignored md quantities because it is the same to all)
        corder = (centroids * md).sum(axis=1).sort_values() # cluster order
        coidx = corder.index
        
        dorder =(df3.loc[afs,] * md).sum(axis=1).sort_values() # data order
        
        # get first fs table    
        fig_properties = {"xlabel" : "Conditions", "ylabel" : "Features", "xlog_mode" : "none", "ylog_mode" : "none", "title" : "Log Fold Changes", "plot_type" : "heatmap", 'ygroup': []}
        fig_properties['ygtick_labels'] = coidx.tolist()

        if 'fold_change' in param and param['fold_change'] == 1:
            frange = 2
            if 'fold_change_range' in param:
                frange = float(param['fold_change_range'])
            final=fc_df.loc[dorder.loc[cl[coidx[0]],].index,]
            fig_properties['ygroup'].append(final.shape[0])
            
            for i in range(1,len(coidx)):
                tf = fc_df.loc[dorder.loc[cl[coidx[i]],].index,]
                fig_properties['ygroup'].append(tf.shape[0])
                final = final.append(tf)

            if 'fold_cutoff' in param and param['fold_cutoff'] == 1:
                final[final > frange] = frange
                final[final < - frange] = - frange
            else:
                fc_df0b = final.sub(final.min(axis=1), axis=0)
                final = (fc_df0b.div(fc_df0b.max(axis=1), axis=0) - 0.5) * 2 * frange
        else:
            final=df2.loc[dorder.loc[cl[coidx[0]],].index,]
            fig_properties['ygroup'].append(final.shape[0])
            
            for i in range(1,len(coidx)):
                tf = df2.loc[dorder.loc[cl[coidx[i]],].index,]
                fig_properties['ygroup'].append(tf.shape[0])
                final = final.append(tf)
        
        ## loading pvalue distribution FDT
        fdt = {'row_labels' :[], 'column_labels' : [], "data" : [[]]};
        #fdt = OrderedDict(fdt)
        # Nan to None
        final = final.where(pd.notnull(final),None)
        fdt['data'] = final.T.as_matrix().tolist() # make sure Transpose
        fdt['row_labels'] = final.columns.tolist()
        fdt['column_labels'] = final.index.tolist()
        # TODO: Add group label later
        fdt['id'] = "{0}.fdt".format(param['out_figure_object_name'])
 
        self.logger.info("Saving the results")
        sstatus = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'MAK.FloatDataTable',
                                                                              'data' : fdt,
                                                                              'hidden':1, 
                                                                              'name' : "{0}.fdt".format(param['out_figure_object_name'])}]})

        data_ref = "{0}/{1}/{2}".format(sstatus[0][6], sstatus[0][0], sstatus[0][4])
        fig_properties['data_ref'] = data_ref

        sstatus = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'CoExpression.FigureProperties',
                                                                              'data' : fig_properties,
                                                                              #'hidden':1, 
                                                                              'name' : "{0}".format(param['out_figure_object_name'])}]})
                                                                              #'name' : "{0}.fp".format(param['out_figure_object_name'])}]})

        #mchp = {}
        #mchp['figure_obj'] = "{0}/{1}/{2}".format(sstatus[0][6], sstatus[0][0], sstatus[0][4])
        #sstatus = ws.save_objects({'workspace' : workspace_name, 'objects' : [{'type' : 'CoExpression.MulticlusterHeatmapPlot',
        #                                                                      'data' : mchp,
        #                                                                      'name' : (param['out_figure_object_name'])}]})

        result = fig_properties
        #END view_heatmap

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method view_heatmap return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
