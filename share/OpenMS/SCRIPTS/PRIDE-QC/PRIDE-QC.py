#!/usr/bin/env python
__author__ = 'walzer'
import sys
import os
import subprocess
from datetime import datetime
#import pickle
import logging
import argparse

VERSION = "0.1"


def create_qcml(mzid, spectra, directory):
    """

    :param mzid:
    :param spectra:
    """
    ott = str(directory) + '/' + os.path.splitext(os.path.basename(spectra))[0]
    exe = "QCCalculator -in {inp} -out {out} -id {id}".format(inp=spectra, out=ott, id=mzid)
    logging.warning("Starting to create qcml for  " + str(mzid) + ', ' + str(spectra) + " ...")
    try:
        o = subprocess.check_output(exe, stderr=subprocess.STDOUT, shell=True)
    except Exception as e:
        logging.warning("Could not create qcml for  " + str(mzid) + ', ' + str(spectra) + " - skipping.")


def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-id', '--mzid', dest="mzid_file", help='<Required> full path to the input mzid', required=True)
    parser.add_argument('-sp', '--spectras', nargs='+', dest="spectra_list",
                        help='<Required> list of the spectra carrying files referenced in the given mzid - each full path', required=True)
    parser.add_argument('-od', "--outdir", dest="outfiles_directory", help="<Required> Outfile for qcml", required=True)
    parser.add_argument('-of', "--outfile", dest="outfile_path", help="Outfile for final qcml.")

    #parser.add_option("--outdir", "-od", action="store", dest="outfile_directory", help="Outfile for predictions")
    options = parser.parse_args()
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)

    if not (options.mzid_file or options.spectra_list or options.outfile_directory):
        parser.print_help()
        sys.exit(1)
    #~ else:
    #~ print dir(options)

    logging.basicConfig(filename=options.outfile_directory + str(datetime.now()) + '.log', filemode='w+',
                        level=logging.DEBUG)  #, format='%(levelname)s:%(message)s'
    logging.info("Starting PRIDE-QC for " + options.mzid_file + " at " + str(datetime.now()))
    args = parser.parse_args()
    logging.warning("verbosity turned on")

    #print options.mzid_file, options.spectra_list, options.outfile_directory
    # TODO go through all spectra files and create one by one the qcmls
    for mse in options.spectra_list:
        create_qcml(mzid=options.mzid_file, spectra=mse, dir=options.outfile_directory)


if __name__ == '__main__':
    __main__()
