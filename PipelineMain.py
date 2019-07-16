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
from itertools import chain
    
from subprocess import call
from collections import defaultdict
from lib.mapper import Mapper
from lib.confdict import ConfDict
from lib.splitter import Splitter
from lib.expander import Expander
    
class PipelineII:   
    def __init__(self, options):

        self.options = options
        basepath = os.path.dirname(os.path.realpath(__file__))
        confd = ConfDict(options, basepath)
        confd.loadConfig()
        self.config_log_out = confd.dumpConfigLog()
        self.confd = confd

        self.mapo  = Mapper(confd)

    def exe_bestacc(self):
        confd = self.confd
        ldel = confd.ldelim
        Key_path = confd.Key_path
        config_file = confd.config_file
        Merge_dir = confd.Merge_dir
        bestacc_script = confd.bestacc_script
        basepath = confd.basepath
        Script_dir = basepath + ldel + "shell_scripts"
        bestacc_cmd = "sh " + bestacc_script + " " + Key_path + " " + Merge_dir + " " + Merge_dir + " " + config_file + " " + Script_dir
        print("Starting bestacc: " + bestacc_cmd)
        call(bestacc_cmd.split())

    
    def build_dict(self):
        confd = self.confd
        ldelim = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        ldict_infile  = confd.Dict_infile
        UGER_cbp_dir = confd.UGER_cbp_dir
        Data_dir = confd.Data_dir
        Log_dir = confd.Log_dir
        use_qsub = confd.use_qsub
        dict_builder = confd.dict_builder

        ldict_infile2 = os.path.basename(os.path.normpath(ldict_infile))
        ldict_infile3 = ldict_infile2.replace(".txt", "")
        ldict_file1 = ldict_infile3 + ".dat"
        ldict_file = Log_dir + ldelim + ldict_file1 
        dict_builder_cmd = dict_builder + " -i " + ldict_infile + " -o " + \
            ldict_file 
        ldict_build_job_path = Log_dir + ldelim + "ldict_build_job.txt"
        ldict_build_jf = open(ldict_build_job_path, "w")
        ldict_build_jf.write(dict_builder_cmd + "\n")
        ldict_build_jf.close()
        ldict_build_cmd = UGER_cbp + " --cmds_file " + ldict_build_job_path + \
                                " --batch_size 1" + \
                                " --memory 1" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header" 
        if use_qsub:
            call(ldict_build_cmd.split())
        print(ldict_build_cmd)
        return ldict_file
   
    def create_soft_links(self, prefix_lst):
        print("Creating soft links from input directory.")
        confd = self.confd
        splito = Splitter(confd)
        for lprefix in prefix_lst:
            splito.create_soft_links(lprefix) 
    
    def create_split_job_files(self, prefix_lst, ldict_file):
        confd = self.confd
        ldelim = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        UGER_cbp_dir = confd.UGER_cbp_dir
        Log_dir = confd.Log_dir
        use_qsub = confd.use_qsub

        splito = Splitter(confd, ldict_file)
        joblist_path = Log_dir + ldelim + "bc_split_joblist.txt"
        jfile = open(joblist_path, "w")
        for lprefix in prefix_lst:
            cmd_str2 = splito.get_split_command(lprefix)
            print(cmd_str2)
            jfile.write(cmd_str2)

        jfile.close()
        joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --memory 8" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())
 
        return joblist_path


    def exe_merge(self, ltab, sample_prefix, sample_bc):
        
        confd = self.confd
        mapo = self.mapo
 
        config_log_out = self.config_log_out
        UniMerger = confd.UniMerger
        Log_dir = confd.Log_dir
        ldelim = confd.ldelim
        project_set = confd.project_set
        use_qsub = confd.use_qsub
        UGER_cbp_dir = confd.UGER_cbp_dir
        UGER_cbp = confd.UGER_cbp

        sample_project = self.sample_project

        unicore_job_path = Log_dir + ldelim + "unimerger_joblist.txt"
        jfile = open(unicore_job_path, "w")

        for lsample in sample_prefix:
            lproject_id = sample_project[lsample]
            prefix_set = sample_prefix[lsample]
            prefix_str = ';'.join(prefix_set)
            bc_set = sample_bc[lsample]
            bc_str = ';'.join(bc_set)
            if not lproject_id in project_set:
                continue
            
            cmd_str = UniMerger + " --config_log_file " + config_log_out + \
              " --sample_id \"" + lsample + "\" --project_id " + lproject_id + \
              ' --prefix_set "' + prefix_str + '" --bc_set "' +  bc_str + '"'

            loutfile = Log_dir + ldelim + lsample + "_out.txt"
            lerrfile = Log_dir + ldelim + lsample + "_err.txt"
            cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
            print(cmd_str2)
            jfile.write(cmd_str2)
        jfile.close()
        lmemory = "8"
        joblist_cmd = UGER_cbp + " --cmds_file " + unicore_job_path + \
                                " --batch_size 1" + \
                                " --num_cores 1" + \
                                " --memory " + lmemory + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())
        return unicore_job_path

  
    def create_temp_bamdir(self, lproject_id, subdir = ''):
        confd = self.confd
        ldelim = confd.ldelim
        Patho_dir = confd.Patho_dir
        patho_project_id = Patho_dir + ldelim + lproject_id
        if subdir:
            patho_project_id += ldelim + subdir
        patho_temp_bamdir = patho_project_id + ldelim + "temp_bamdir"
        if not os.path.exists(patho_temp_bamdir):
            print("Creating temp_bam directory:" + patho_temp_bamdir)
            os.makedirs(patho_temp_bamdir)


    def create_picard_dir(self, lproject_id, species_type):
        confd = self.confd
        ldel = confd.ldelim
        picard_outdir = None
        Results_path = confd.Results_path
        picard_project = Results_path + ldel + lproject_id
        picard_prefix = picard_project + ldel + "bam_metrics"
        if "patho" == species_type:
            picard_outdir = picard_prefix + "_patho"
        else:
            raise StandardError("Unknown species type: " + species_type)
            
        if not os.path.exists(picard_outdir):
            print("Creating picard metrics dir: " + picard_outdir)
            os.makedirs(picard_outdir)
        
          
    def exe_splitting(self, KeyTblFinal, prefix_lst):
        confd = self.confd
        do_split = confd.do_split
        use_qsub = confd.use_qsub
        soft_link_from_input = confd.soft_link_from_input

        if do_split:
            if soft_link_from_input:
                self.create_soft_links(prefix_lst)
            else:
                ldict_file = self.build_dict()
                joblist_path = self.create_split_job_files(prefix_lst, ldict_file)


    def copy_summary_to_result(self, summary_tag, lproject_id):
        confd = self.confd
        ldelim = confd.ldelim
        Results_path = confd.Results_path
        Summary_dir = confd.Summary_dir
        lpatho_summary_dir = Summary_dir + ldelim + lproject_id
        results_dir = Results_path + ldelim + lproject_id
        src_file = lpatho_summary_dir + ldelim + summary_tag
        dest_file = results_dir + ldelim + summary_tag
        print("Copying metrics file")
        print("src_file: " + src_file)
        print("dest_file: " + dest_file)
        copyLargeFile(src_file, dest_file) 


    def copy_moc_files(self):
        confd = self.confd
        Results_path = confd.Results_path
        ldelim = confd.ldelim
        # At this moment we shall copy a file to the results folder, this 
        # would be a hardcoded path.
        moc_dir = "/broad/IDP-Dx_storage/MOC/files"
        moc_file = "MOC_RNA-Seq_data_UserGuide_CURRENT.docx"
 
        moc_file_path = moc_dir + ldelim + moc_file
        result_dir = Results_path + ldelim + lproject_id
        moc_file_out_path = result_dir + ldelim + moc_file
        print("Copying: " + moc_file_path + " to " + moc_file_out_path)
        copyLargeFile(moc_file_path, moc_file_out_path)
         

    def getKeyTblFinal(self):
        confd = self.confd
        KeyTblFinal_ori = confd.getKeyTblFinal()
        print("KeyTblFinal is loaded!")
        do_expand = confd.do_expand
        if not do_expand:
            return KeyTblFinal_ori
        else: 
            # Expand the table if required
            expo = Expander(confd)
            KeyTblFinal = expo.expandKeyTbl(KeyTblFinal_ori)
            return KeyTblFinal
       

    def mainFunc(self):
        confd = self.confd
        mapo = self.mapo

        do_split = confd.do_split
        do_merge = confd.do_merge

        KeyTblFinal = self.getKeyTblFinal()
        prefix_lst = mapo.get_prefix_set(KeyTblFinal)
        prefix_map = mapo.build_prefix_map(prefix_lst, KeyTblFinal)
        print(prefix_lst)

        mapo = self.mapo
        sample_prefix = mapo.create_sample_map_prefix(KeyTblFinal, prefix_map)
        print("sample_prefix:")
        print(sample_prefix)
        bc_delim = '\s*;\s*'
        sample_bc = mapo.create_sample_map_bc(KeyTblFinal, bc_delim)

        sample_project = mapo.get_sample_project(KeyTblFinal)
        sample_refacc = mapo.create_sample_refacc(KeyTblFinal) 
        
        self.sample_project = sample_project
        self.sample_refacc = sample_refacc
        
        print(sample_prefix)
        if do_split:
            self.exe_splitting(KeyTblFinal, prefix_lst)
        if do_merge:
            self.exe_merge(KeyTblFinal, sample_prefix, sample_bc)
             

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process the options.')
    parser.add_argument('--config_file', '-c', type = str, required = True, help ='Path to main config file')
    parser.add_argument('--key_id', dest = 'key_id', required = False, default = 'none', help = 'Key file would be "key_id"_key.txt')
    parser.add_argument('--key_path', dest = 'key_path', required = False, default = 'none', help = 'Key file path (absolute)')
    parser.add_argument('--project_ids', required = False, default = 'none', help = 'Project ids')
    parser.add_argument('--seq_id', required = False, default = 'none', help = 'Sequencing id used to make raw_seq_path')
    parser.add_argument('--raw_seq_path', required = False, default = 'none', help = 'Directory for raw sequene files (absolute)')
    parser.add_argument('--temp_path', required = False, default = 'none', help = 'Will contain the temporary results')
    parser.add_argument('--bam_path', required = False, default = 'none', help = 'Will contain the aligned sorted bam files')
    parser.add_argument('--results_path', required = False, default = 'none', help = 'Will contain the path to the results')
    parser.add_argument("--do_patho", required=False, action="store_true", default=False, help="Do the patho")
    parser.add_argument("--remove_split", required=False, action="store_true", default=False, help="Remove the split files")
    parser.add_argument("--gzip_merged", required = False, action = "store_true", default = False, help = "Compress the merged files to gzip format")
    parser.add_argument('--read_counter', dest = "read_counter", type = str, \
        default = 'JL_counter', help = 'Tool for counting reads (default: Counter written by JLivny)')
    parser.add_argument('--no_p7', dest = 'use_p7', action = 'store_false', default = True, help = 'Use if P7 index is not used.' )
    parser.add_argument('--use_p5', dest = 'use_p5', action = 'store_true', default = False, help = 'Use if P5 index is used.' )
    parser.add_argument('--use_lane', dest = 'use_lane', action = 'store_true', default = False, help = 'Use if lane specific merging is required.' )
    parser.add_argument('--use_seq_path', dest = 'use_seq_path', action = 'store_true', default = False, help = 'Use if direct mapping of sample to raw seq path is required.' )
    parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )
    parser.add_argument('--no_split', dest = 'do_split', action = 'store_false', default = True, help = 'Does not split the fastq files.' )
    parser.add_argument('--no_merge', dest = 'do_merge', action = 'store_false', default = True, help = 'Does not merge the split fastq files.' )
    parser.add_argument('--no_replace_refname', dest = 'do_replace_refname', action = 'store_false', default = True, help = 'Does not replace reference name in fna')
    parser.add_argument('--no_expand', dest = 'do_expand', action = 'store_false', default = True, help = 'Does not expand p7 and barcode entries in the key file if required.')
    parser.add_argument('--no_bc_split', dest = 'no_bc_split', action = 'store_true', default = False, help = 'Does not run the barcode splitter, instead create softlinks of the raw seq files to the split directory.')
    parser.add_argument('--no_ref', dest = 'do_ref', action = 'store_false', default = True, help = 'Does not generate patho ref.')
    parser.add_argument('--Suffix_s1', default = 'none', required = False, help = 'Update the value of Suffix_s1')
    parser.add_argument('--Suffix_s2', default = 'none', required = False, help = 'Update the value of Suffix_s2')
    parser.add_argument('--Suffix_ne', default = 'none', required = False, help = 'Update the value of Suffix_ne')
    parser.add_argument('--ADD5', dest = 'add5', type = int, default = 0, help = 'ADD5 for gff parser') 
    parser.add_argument('--ADD3', dest = 'add3', type = int, default = 0, help = 'ADD3 for gff parser') 
    parser.add_argument('--MOC_id', dest = 'MOC_id', type = str, default = 'none', help = 'Provide MOC string for adding MOC hierarchy')
    parser.add_argument('--trim_rs_5p', dest = 'trim_rs_5p', type = int, default = 0, help = '5p trim count for read-single') 
    parser.add_argument('--trim_rs_3p', dest = 'trim_rs_3p', type = int, default = 0, help = '3p trim count for read-single') 
    parser.add_argument('--keep_rs_5p', dest = 'keep_rs_5p', type = int, default = -1, help = '5p keep count for read-single') 
    parser.add_argument('--keep_rs_3p', dest = 'keep_rs_3p', type = int, default = -1, help = '3p keep count for read-single') 
    parser.add_argument('--trim_r1_5p', dest = 'trim_r1_5p', type = int, default = 0, help = '5p trim count for read1') 
    parser.add_argument('--trim_r1_3p', dest = 'trim_r1_3p', type = int, default = 0, help = '3p trim count for read1') 
    parser.add_argument('--trim_r2_5p', dest = 'trim_r2_5p', type = int, default = 0, help = '5p trim count for read2') 
    parser.add_argument('--trim_r2_3p', dest = 'trim_r2_3p', type = int, default = 0, help = '3p trim count for read2') 
    parser.add_argument('--keep_r1_5p', dest = 'keep_r1_5p', type = int, default = -1, help = '5p keep count for read1') 
    parser.add_argument('--keep_r1_3p', dest = 'keep_r1_3p', type = int, default = -1, help = '3p keep count for read1') 
    parser.add_argument('--keep_r2_5p', dest = 'keep_r2_5p', type = int, default = -1, help = '5p keep count for read2') 
    parser.add_argument('--keep_r2_3p', dest = 'keep_r2_3p', type = int, default = -1, help = '3p keep count for read2') 
    parser.add_argument('--MOC_id_ref', dest = 'MOC_id_ref', type = str, default = 'none', help = 'Adds MOC_id_ref to Bacterial_Ref_path')
    parser.add_argument("--no_login_name", dest = "use_login_name", action = 'store_false', default = True, help = 'Generate results in a username specific directory')
    parser.add_argument("--count_strand_rev", dest = "count_strand_rev", type = str, default = "N", help = 'Reverse the counting mechanism by making it "forward"')
    parser.add_argument("--do_bestacc", dest = "do_bestacc", action = "store_true", default = False, help = "Run bestacc after splitting and merging")
    parser.add_argument("--use_sample_id", dest = "use_sample_id", action = "store_true", default = False, help = "Use sample id for sample id mapping")
    parser.add_argument("--bwa_mem", dest = "bwa_mem", action = "store_true", default = False, help = "Use bwa mem instread of bwa backtrack")
    parser.add_argument("--rm_rts_dup", dest = "rm_rts_dup", action = "store_true", default = False, help = "Remove PCR duplicates from RNA-tagseq")
    parser.add_argument("--do_rerun", dest = "do_rerun", action = "store_true", default = False, help = "Rerun pipeline after a previous incomplete run")
    parser.add_argument("--paired_only_patho", dest = "paired_only_patho", action = "store_true", default = False, help = "Align/count only the reads when both ends are mapped on pathogen side")
    parser.add_argument('--ucore_time', dest = 'ucore_time', type = int, default = 0, help = 'Timelimit of unicore for this run') 
    
    
    
    options = parser.parse_args()

    pipel = PipelineII(options) 
    pipel.mainFunc()

