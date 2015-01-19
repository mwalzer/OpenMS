#!/usr/bin/env python
__author__ = 'walzer'
import sys
import os
import subprocess
from datetime import datetime
import logging
import argparse
#import pickle

VERSION = "0.1"

IDRATIO_RPLOT = "ProduceQCFigures_idmap.R"
TIC_RPLOT = "ProduceQCFigures_idmap.R"
MASSERR_RPLOT = "ProduceQCFigures_idmap.R"


def create_qp(qcml, qp_cv, type_cv, temp, rscript, qpcv_export_list):
    export_files = list()
    try:
        for export_qccv in qpcv_export_list:
            export_path = str(temp) + '/' + str(export_qccv)
            eexe = "QCExporter -in {inpu} -out {out} -qp {id}".format(inpu=qcml, out=export_path, qp=str(export_qccv))
            o = subprocess.check_output(eexe, stderr=subprocess.STDOUT, shell=True)
            # TODO possible to check if output is legit?!
    except Exception as e:
        logging.warning("Could not find " + str(export_qccv) + ' in qcml ' + str(qcml) + " - no " + str(rscript) + " QC metric")
        return None

    plot_path = str(temp) + '/' + str(qp_cv) + ".png"
    try:
        rexe = "RScript {inpu} {out}".format(inpu=' '.join(export_files), out=plot_path)
        o = subprocess.check_output(rexe, stderr=subprocess.STDOUT, shell=True)
    except Exception as e:
        logging.warning("Could not execute R " + str(rscript) + ' for ' + str(qcml) + " - no " + str(rscript) + " QC metric")
        return None

    try:
        eexe = "QCEmbedder -in {inpu} -out {out} -plot {plot} -qp_att_acc {type} -cv_acc {cv}".format(inpu=qcml, out=qcml, plot=plot_path, type=type_cv, cv=qp_cv)
        o = subprocess.check_output(eexe, stderr=subprocess.STDOUT, shell=True)
    except Exception as e:
        logging.warning("Could not embed plot for " + str(rscript) + ' into ' + str(qcml) + " - no " + str(rscript) + " QC metric")
    return qcml


def create_qcml(mzid, spectra, directory):
    qcml = str(directory) + '/' + os.path.splitext(os.path.basename(spectra))[0]
    exe = "QCCalculator -in {inpu} -out {out} -id {id}".format(inpu=spectra, out=qcml, id=mzid)
    logging.warning("Starting to create qcml for  " + str(mzid) + ', ' + str(spectra) + " ...")
    try:
        o = subprocess.check_output(exe, stderr=subprocess.STDOUT, shell=True)
    except Exception as e:
        logging.warning("Could not create qcml for  " + str(mzid) + ', ' + str(spectra) + " - skipping.")
    return qcml


def __main__():
    parser = argparse.ArgumentParser(version=VERSION)
    parser.add_argument('-id', '--mzid', dest="mzid_file", help='<Required> full path to the input mzid', required=True)
    parser.add_argument('-sp', '--spectras', nargs='+', dest="spectra_list",
                        help='<Required> list of the spectra carrying files referenced in the given mzid - each full path', required=True)
    parser.add_argument('-od', "--outdir", dest="outfiles_directory", help="<Required> Outfile for qcml", required=True)
    parser.add_argument('-of', "--outfile", dest="outfile_path", help="Outfile for final qcml.")
    parser.add_argument('-temp', "--tempdir", dest="temp_path", help="Temporary directory for intermediate files.")

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
        qcml = create_qcml(mzid=options.mzid_file, spectra=mse, dir=options.outfile_directory)
        create_qp(qcml, options.temp_path, IDRATIO_RPLOT, {'precursors' : 'QC:0000044', 'ids' : 'QC:0000038'})


if __name__ == '__main__':
    __main__()
