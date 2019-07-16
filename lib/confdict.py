#!/usr/bin/env python

import yaml
import argparse
import shutil
import pandas as pd
import csv
import re
import os
import ntpath
import codecs
import getpass

class ConfDict(object):

    def __init__(self, options, basepath):
        # Store all the properties to the dict "datastore"
        object.__setattr__(self, 'datastore', {})
        self.datastore['basepath'] = basepath

        self.options = options
        self.config_file = options.config_file
        self.mydict = yaml.load(open(self.config_file))
        self.ldelim = "/"
        self.inidict = dict()
        self.popu_inidict()

    def popu_inidict(self):
        inidict = self.inidict
        inidict['LC_method'] = 'LC_method'
        inidict['Read_pairing'] = 'Read_pairing'

    def __getattr__(self, key):
        if key in ['KeyTbl', 'options', 'mydict']:
            return super(ConfDict, self).__getattr__(key)
        else:
            return self.datastore[key]

    def __setattr__(self, key, value):
        if key in ['KeyTbl', 'options', 'mydict']:
            super(ConfDict, self).__setattr__(key, value)
        else:
            self.datastore[key] = value

    def get_project_set(self, ltab):
        project_id = self.project_id
        project_lst = ltab[project_id].tolist()
        project_set = list(set(project_lst))
        return project_set

    def loadKeyTable(self):
        Key_path = self.Key_path
        sample_id = self.sample_id
        lfile = codecs.open(Key_path, encoding='utf-8')
        tab0 = csv.reader(lfile, delimiter = '\t')
        tab1 = [row for row in tab0]
        tab2 = filter(lambda x: not x[0].startswith("###"), tab1)
        header = tab2.pop(0)	
        tab3 = pd.DataFrame(tab2, columns = header)
        for index, row in tab3.iterrows():
            lsample = row[sample_id]
            lsample2 = re.sub('\s+', '_', lsample)
            lsample3 = re.sub('/', '_', lsample2)
            tab3.loc[index, sample_id] = lsample3
        self.KeyTbl = tab3

    def storeConfigFromKeyTbl(self):
        KeyTbl = self.KeyTbl
        options = self.options
        mydict = self.mydict

        # Here we are modifying the code a little and still trying to keep it
        # backward compatible. If the LC_method_val obtained from the key file
        # ends with _alr1 ot _alr2, then set those flags on and keep that part
        # separate from the rest of LC_method_val

        LCmv_ori = KeyTbl[self.LC_method][0]
        self.align_readnum = ''
        LCmv_ori_u = LCmv_ori.lower()
        if LCmv_ori_u.endswith('_alr1') or LCmv_ori_u.endswith('_alr2'):
            lparts = LCmv_ori.rsplit('_', 1)
            self.LC_method_val = lparts[0]
            self.align_readnum = lparts[1].lower() 
        else:
            self.LC_method_val = LCmv_ori

        Read_pairing_val_ci = KeyTbl[self.Read_pairing][0]
        if not Read_pairing_val_ci:
            Read_pairing_val_ci = 'Paired'
            print("Rad_pairing_val empty. setting it: " + Read_pairing_val_ci)
        Read_pairing_val = Read_pairing_val_ci.upper()
        self.Read_pairing_val = Read_pairing_val.upper()

         

        self.soft_link_from_input = False
        lbc_splitter = ''
        lc_lower = self.LC_method_val.lower()
        if lc_lower == 'rts' or lc_lower == 'rts-ts':
            if Read_pairing_val == 'PAIRED':
                lbc_splitter = self.bc_splitter_rts
            elif Read_pairing_val == 'SINGLE':
                lbc_splitter = self.bc_splitter_rts_se
            self.Dict_infile = self.RtS_dict_file
        elif lc_lower == 'alignr2':
            self.Dict_infile = self.RtS_dict_file
        elif lc_lower == 'alignr1':
            self.Dict_infile = self.RtS_dict_file


        no_bc_split = self.no_bc_split
        if no_bc_split:
            # This may not be required: Disturbing the dropseq
            #self.soft_link_from_input = True
            lbc_splitter = 'None'
                

        do_R2_trim = False
        if lc_lower == 'rts-ts':
            do_R2_trim = True


        self.lbc_splitter = lbc_splitter
        self.do_R2_trim = do_R2_trim
        self.lc_lower = lc_lower


    def get_from_mydict(self, key):
        mydict = self.mydict
        inidict = self.inidict

        if key in mydict:
            return mydict[key]
        elif key in inidict:
            return inidict[key]
        else:
            return '__NA__'

    def storeConfigFromConfig(self):
        mydict = self.mydict

        self.project_id      = self.get_from_mydict('Proj')
        self.P7_index        = self.get_from_mydict('P7')
        self.P5_index        = self.get_from_mydict('P5')
        self.Lane_index      = self.get_from_mydict('Lane')
        self.sample_id       = self.get_from_mydict('ID')
        self.LC_method       = self.get_from_mydict('LC_method')
        self.Path_to_SeqFile = self.get_from_mydict('Path_to_SeqFile')
        self.Ref_accession   = self.get_from_mydict('Ref_accession')
        self.Host_reference  = self.get_from_mydict('Host_reference')
        self.lbc             = self.get_from_mydict('bc')         
        self.ldelim          = self.get_from_mydict('delim')
        self.dict_builder    = self.get_from_mydict('dict_builder')
        self.bc_splitter     = self.get_from_mydict('bc_splitter')
        self.bc_splitter_rts = self.get_from_mydict('bc_splitter_rts')
        self.bc_splitter_rts_se = self.get_from_mydict('bc_splitter_rts_se')
        self.sam_fragcount   = self.get_from_mydict('sam_fragcount')
        self.paired_only_script = self.get_from_mydict('paired_only_script')
        self.frag_to_gene_count = self.get_from_mydict('frag_to_gene_count')
        self.metrics_gen     =  self.get_from_mydict('metrics_gen')
        self.UGER_cbp        = self.get_from_mydict('UGER_cbp')
        self.bwa_path        = self.get_from_mydict('bwa')
        self.patho_dbpath    = self.get_from_mydict('patho_dbpath')
        self.Read_pairing    = self.get_from_mydict('Read_pairing')
        self.samtools       = self.get_from_mydict('samtools')
        self.featureCounts  = self.get_from_mydict('featureCounts')
        self.JLCounter      = self.get_from_mydict('JLCounter')
        self.strand_dir     = self.get_from_mydict('strand_dir')
        self.tdf_str        = self.get_from_mydict('tdf_str')
        self.patho_thread_count = self.get_from_mydict('patho_thread_count')
        self.patho_memory   = self.get_from_mydict('patho_memory')
        self.picard_bindir = self.get_from_mydict('picard_bindir')
        self.cutadapt       = self.get_from_mydict('cutadapt')   
        self.RtS_dict_file = self.get_from_mydict('RtS_dict_file')
        self.Suffix_s1 = self.get_from_mydict('Suffix_s1')
        self.Suffix_s2 = self.get_from_mydict('Suffix_s2')
        self.Suffix_ne = self.get_from_mydict('Suffix_ne')
        self.read_trimmer_path = self.get_from_mydict('read_trimmer')
        self.sam_splitter_path = self.get_from_mydict('sam_splitter')
        self.Results_path = self.get_from_mydict('Results_path')
        self.Data_dir = self.get_from_mydict('Data_dir')
        
        # Added as a fix for issue 41 (https://github.com/broadinstitute/PipelineII/issues/41)
        self.Pipeline_bcLog = self.get_from_mydict('Pipeline_bcLog')
        self.corr_script = self.get_from_mydict('corr_script')
        self.gff_parser = self.get_from_mydict('gff_parser')
        self.RPG_metrics_script = self.get_from_mydict('RPG_metrics_script')
        self.Data_finish_script = self.get_from_mydict('Data_finish_script')
        self.bcLog_metrics_script = self.get_from_mydict('bcLog_metrics_script')
        self.trim_script = self.get_from_mydict('trim_script')
        self.picard_metrics = self.get_from_mydict('picard_metrics')
        self.picard_metrics_parse = self.get_from_mydict('picard_metrics_parse')
        self.fpkm_script = self.get_from_mydict('fpkm_script')
        self.bestacc_script = self.get_from_mydict('bestacc_script')
        self.remove_dup_script = self.get_from_mydict('remove_dup_script')
        

    def storeDerivedPaths(self):
        ldelim = self.ldelim
        Script_dir = self.basepath 
        self.UniCore      = Script_dir + ldelim + "lib" + ldelim + "unicore.py" 
        self.UniMerger      = Script_dir + ldelim + "lib" + ldelim + "unimerger.py" 


    def storeConfigFromOptions(self):
        options = self.options
        self.Project_ids = options.project_ids
        self.do_replace_refname = options.do_replace_refname
        self.add5 = options.add5
        self.add3 = options.add3
        self.trim_rs_5p = options.trim_rs_5p
        self.trim_rs_3p = options.trim_rs_3p
        self.keep_rs_5p = options.keep_rs_5p
        self.keep_rs_3p = options.keep_rs_3p
        self.trim_r1_5p = options.trim_r1_5p
        self.trim_r1_3p = options.trim_r1_3p
        self.trim_r2_5p = options.trim_r2_5p
        self.trim_r2_3p = options.trim_r2_3p
        self.keep_r1_5p = options.keep_r1_5p
        self.keep_r1_3p = options.keep_r1_3p
        self.keep_r2_5p = options.keep_r2_5p
        self.keep_r2_3p = options.keep_r2_3p
        self.read_counter = options.read_counter
        self.use_p7 = options.use_p7
        self.use_p5 = options.use_p5
        self.use_lane = options.use_lane
        self.do_patho = options.do_patho
        self.remove_split = options.remove_split
        self.gzip_merged = options.gzip_merged
        self.use_qsub = options.use_qsub
        self.do_split = options.do_split
        self.do_merge = options.do_merge
        self.do_expand = options.do_expand
        self.no_bc_split = options.no_bc_split
        self.count_strand_rev = options.count_strand_rev
        self.use_seq_path = options.use_seq_path
        self.MOC_id = options.MOC_id
        self.MOC_id_ref = options.MOC_id_ref
        self.use_sample_id = options.use_sample_id
        self.do_rerun = options.do_rerun
        self.ucore_time = options.ucore_time
        self.paired_only_patho = options.paired_only_patho


    def storeConfigMixed(self):
        options = self.options
        mydict = self.mydict
        ldelim = self.ldelim

        suffix_s1 = ''
        if options.Suffix_s1 != 'none':
            suffix_s1 = options.Suffix_s1
        else:
            suffix_s1 = self.Suffix_s1
        self.suffix_s1 = suffix_s1
        
        suffix_s2 = ''
        if options.Suffix_s2 != 'none':
            suffix_s2 = options.Suffix_s2
        else:
            suffix_s2 = self.Suffix_s2
        self.suffix_s2 = suffix_s2

        suffix_ne = ''
        if options.Suffix_ne != 'none':
            suffix_ne = options.Suffix_ne
        else:
            suffix_ne = self.Suffix_ne
        self.suffix_ne = suffix_ne

        if self.use_sample_id:
            self.use_p7 = False
            self.no_bc_split = True
            self.soft_link_from_input = True

        regex = '\s*[,|;]\s*'
        raw_seq_path = list()
        if options.seq_id != 'none':
            Seq_base = mydict['Seq_base']
            lseq_ids = re.split(regex, options.seq_id)
            for lseq_id in lseq_ids:
                raw_seq_path_l = Seq_base + self.ldelim + lseq_id
                raw_seq_path.append(raw_seq_path_l)
        elif options.raw_seq_path != 'none':
            raw_seq_path = re.split(regex, options.raw_seq_path)
        else:
            raise argparse.ArgumentTypeError('Either seq_id or raw_seq_path need to be provided.')
        print(raw_seq_path)
        self.Input_dir = raw_seq_path 


        self.login_name = getpass.getuser()
        temp_path = ''
        if options.temp_path != 'none':
            temp_path = options.temp_path
        else:
            temp_path = mydict['Temp_path']

        if options.use_login_name:
            self.Temp_dir = temp_path + self.ldelim + self.login_name
        else:
            self.Temp_dir = temp_path
        
        bam_path = ''
        if options.bam_path != 'none':
            bam_path = options.bam_path
        else:
            bam_path = mydict['Bam_path']

        if options.use_login_name:
            self.Bam_path = bam_path + self.ldelim + self.login_name
        else:
            self.Bam_path = bam_path

        results_path = ''
        if options.results_path != 'none':
            results_path = options.results_path 
        else:
            results_path = mydict['Results_path'] 
       
        if options.use_login_name:
            self.Results_path = results_path + self.ldelim + self.login_name
        else:
            self.Results_path = results_path

        # Add MOC _id hierarchy
        if self.MOC_id != 'none':
            self.Temp_dir = self.Temp_dir + self.ldelim + self.MOC_id
            self.Bam_path = self.Bam_path + self.ldelim + self.MOC_id
            self.Results_path = self.Results_path + self.ldelim + self.MOC_id



    def createSubPaths(self):
        Temp_dir = self.Temp_dir
        ldelim = self.ldelim
        use_qsub = self.use_qsub
        do_patho = self.do_patho
        do_rerun = self.do_rerun

        project_set = ''
        regex = '\s*[,|;]\s*'
        if self.Project_ids != 'none':
             project_set = set(re.split(regex, self.Project_ids))
        else:
            project_set = self.get_project_set(self.KeyTbl)
        self.project_set = sorted(project_set)
        out_dir_name = '_'.join(self.project_set)
        Out_dir = Temp_dir + ldelim + out_dir_name
        self.Out_dir = Out_dir
        
        UGER_cbp_dir = Out_dir + ldelim + "UGER_cbp"
        Split_dir = Out_dir + ldelim + "splitdir"
        Merge_dir = Out_dir + ldelim + "mergedir"
        Merge_real_dir = Out_dir + ldelim + "mergedir/real"
        Merge_garbage_dir = Out_dir + ldelim + "mergedir/garbage"
        Log_dir = Out_dir + ldelim + "logdir"

        if use_qsub and os.path.exists(UGER_cbp_dir):
            if not do_rerun:
                shutil.rmtree(UGER_cbp_dir)
        if not os.path.exists(Out_dir):
            os.makedirs(Out_dir)
        if not os.path.exists(Split_dir):
            os.makedirs(Split_dir)
        if not os.path.exists(Merge_dir):
            os.makedirs(Merge_dir)
        if not os.path.exists(Log_dir):
            os.makedirs(Log_dir)

        self.UGER_cbp_dir = UGER_cbp_dir
        self.Log_dir = Log_dir
        self.Split_dir = Split_dir
        self.Merge_dir = Merge_dir

        # Create the host_ref_path 
        self.prefix_dict_out = self.Log_dir + ldelim + "prefix_dict_out.txt"
       
    def initMain(self):
        options = self.options
        mydict = self.mydict
        key_path = ''
        if options.key_id != 'none':
            Key_base = mydict['Key_base']
            key_path = Key_base + self.ldelim + options.key_id + "_key.txt"
        elif key_path != 'none':
            key_path = options.key_path
        else:
            raise argparse.ArgumentTypeError('Either key_id or key_path need to be provided.')    
        self.Key_path = key_path


    def createSubDirs(self):
        do_patho = self.do_patho
        project_set = self.project_set
        Bam_path = self.Bam_path
        ldelim = self.ldelim
        ldel = ldelim


    def doEndUpdate(self):
        if self.MOC_id_ref != 'none':
            ldel = self.ldelim
            self.patho_dbpath = self.patho_dbpath + ldel + self.MOC_id_ref
            
        
    def loadConfig(self):
        self.initMain()
        self.storeConfigFromConfig()
        self.storeConfigFromOptions()
        self.loadKeyTable()    
        self.storeConfigFromKeyTbl()
        self.storeDerivedPaths()
        self.storeConfigMixed()
        self.createSubPaths()
        self.createSubDirs()
        self.doEndUpdate()

    def dumpConfigLog(self):
        config_log_out = self.Log_dir + self.ldelim + "config_log_out.txt"
        with open(config_log_out, 'w') as outfile:
            yaml.dump(self.datastore, outfile, default_flow_style=False)
        return config_log_out

    def getKeyTblFinal(self):
        return self.KeyTbl

