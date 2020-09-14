#!/usr/bin/env python                                                                            

""" MultiQC module to parse output from Kallisto """

from __future__ import print_function
from collections import OrderedDict
import os
import logging
import re
import json

from multiqc import config
from multiqc.plots import bargraph

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger                                                                         
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """alignstats module"""

    def __init__(self):
        #Init parent obj
        super(MultiqcModule, self).__init__(name='alignstats', anchor='alignstats',href="https://github.com/jfarek/alignstats",info="Alignstats generates a ton of NGS stats for target capture and WGS data")

        self.dat_json = dict()
        for f in self.find_log_files('alignstats', filehandles=False):
            self.parse_alignstats_data(f)
            print('xxxx', f)

        # Write parsed report data to a file
        self.write_data_file(self.dat_json, 'multiqc_alignstats')
        self.alignstats_general_stats_table(f['s_name'])
        self.add_report_section()
        

    def parse_alignstats_data(self,f):

        self.dat_json[f['s_name']] =  {}

        fdat = ""
        file_loc = "{0}/{1}".format(f['root'], f['fn'])
        for i in open(file_loc):
            fdat = "{0}{1}".format(fdat, i.rstrip())
        temp_ds = json.loads(fdat)
        for k in temp_ds:
            self.dat_json[f['s_name']][k] = temp_ds[k]
        
    def alignstats_general_stats_table(self, samp_name):

        headers = OrderedDict()
        headers['MismatchedBasesPct'] = {
            'title':'Mismatch Base Pct',
            'description':'MIS',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['InsertedBasesPct'] = {
            'title':'Inserted Bases Pct',
            'description':'Ins Pct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['DeletedBasesPct'] = {
            'title': 'Deleted Bases Pct',
            'description':'Del Base Pct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['SoftClippedBasesPct'] = {
            'title': 'Soft Clip Pct',
            'description':'Soft Clip Pct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['PerfectReadsPct'] = {
            'title': 'Perfect Reads Pct',
            'description':'Perfect Reads Pct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['R1MappedReadsPct'] = {
            'title': 'R1 Mapped Reads Pct',
            'description':'R1 Mapped Reads Pct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['R1AlignedBasesPct'] = {
            'title': 'R1AlignedBasesPct',
            'description':'R1AlignedBasesPct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['R1PerfectReadsPct'] = {
            'title': 'R1PerfectReadsPct',
            'description':'R1PerfectReadsPct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['R2MappedReadsPct'] = {
            'title': 'R2 Mapped Reads Pct',
            'description':'R2 Mapped Reads Pct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['R2AlignedBasesPct'] = {
            'title': 'R2AlignedBasesPct',
            'description':'R2AlignedBasesPct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }
        headers['R2PerfectReadsPct'] = {
            'title': 'R2PerfectReadsPct',
            'description':'R2PerfectReadsPct',
            'min': 0,
            'suffix' : 'pct',
            'scale': 'RdYlGn'
        }

        data = {
            samp_name: {
                'MismatchedBasesPct' : self.dat_json[samp_name]['MismatchedBasesPct'],
                'InsertedBasesPct' : self.dat_json[samp_name]['InsertedBasesPct'],
                'DeletedBasesPct' : self.dat_json[samp_name]['DeletedBasesPct'],
                'SoftClippedBasesPct' : self.dat_json[samp_name]['SoftClippedBasesPct'],
                'PerfectReadsPct' : self.dat_json[samp_name]['PerfectReadsPct'],
                'R1MappedReadsPct' : self.dat_json[samp_name]['R1MappedReadsPct'],
                'R1AlignedBasesPct' : self.dat_json[samp_name]['R1AlignedBasesPct'],
                'R1PerfectReadsPct' : self.dat_json[samp_name]['R1PerfectReadsPct'],
                'R2PerfectReadsPct' : self.dat_json[samp_name]['R2PerfectReadsPct'],
                'R2AlignedBasesPct' : self.dat_json[samp_name]['R2AlignedBasesPct'],
                'UnmappedBasesPct' : self.dat_json[samp_name]['UnmappedBasesPct'],
                'DuplicateBasesPct' : self.dat_json[samp_name]['DuplicateBasesPct'],
                'MappedBasesPct' : self.dat_json[samp_name]['MappedBasesPct']
            }
        }
        print(headers)
        print(data)
        self.general_stats_addcols(data,headers)

        
    def add_report_section(self):
        #UnmappedBasesPct
        #DuplicateBasesPct
        #MappedBasesPct
        #MismatchedBasesPct
        #InsertedBasesPct
        #DeletedBasesPct
        #SoftClippedBasesPct
        #PerfectReadsPct
        #R1MappedReadsPct
        #R1AlignedBasesPct
        #R1PerfectReadsPct
        #R2PerfectReadsPct
        #R2AlignedBasesPct
        #R2MappedReadsPct

        keys = OrderedDict()
        keys['MismatchedBasesPct'] = { 'color': '#437bb1', 'name': 'MismatchedBasesPct' }
        keys['DeletedBasesPct'] =   { 'color': '#b1084c', 'name': 'DeletedBasesPct' }
        keys['InsertedBasesPct'] =   { 'color': '#b1084c', 'name': 'InsertedBasesPct' }

        # Config for the plot                                                                    
        config = {
            'id': 'alignstats_a',
            'title': 'alignstats Test Plot'
        }


        plot = bargraph.plot(self.dat_json, keys, config)
        self.add_section(name='alignstats', anchor='alignstats', plot= plot)
