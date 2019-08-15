#! /usr/bin/env python

from __future__ import print_function

import unittest
import sys, os, tempfile, shutil, filecmp, tarfile
import iravnet 
from .check_download import check_download
# from .make_savnet_input import *

class TestMain(unittest.TestCase):

    def setUp(self):

        def extract_tar_gz(input_tar_gz_file, out_path):
            tar = tarfile.open(input_tar_gz_file)
            tar.extractall(out_path)
            tar.close()

        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        check_download("https://storage.googleapis.com/friend1ws_package_data/iravnet/PC-7.19.bam", \
                       cur_dir + "/resource/bam/PC-7.19.bam")

        check_download("https://storage.googleapis.com/friend1ws_package_data/iravnet/PC-7.19.bam.bai", \
                       cur_dir + "/resource/bam/PC-7.19.bam.bai")
    
        self.parser = iravnet.parser.create_parser()


    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        
        input_bam = cur_dir + "/resource/bam/PC-7.19.bam"
        output_prefix = tmp_dir + "/PC-7.19.iravnet.result.txt"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"

        iravnet_get_args = ["get", input_bam, output_prefix, ref_genome]
        print("iravnet " + ' '.join(iravnet_get_args))

        args = self.parser.parse_args(iravnet_get_args)
        iravnet.run.get_main(args)

        with open(tmp_dir + "/PC-7.19.iravnet.result.txt", 'r') as hin: record_num = len(hin.readlines())
        self.assertTrue(record_num == 38)
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

