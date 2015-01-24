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

#QC:0000044  -> ms2s.csv
#QC:0000038 -> ids.csv

EXPMAP_RPLOT = {'R': "../ProduceQCFigures_precursormap.R", 'qp_att_acc': 'QC:0000004', 'cv_acc': 'QC:0000055', 'qpcv_export_list': ['QC:0000044']}
IDMAP_RPLOT = {'R': "../ProduceQCFigures_idmap.R", 'qp_att_acc': 'QC:0000035', 'cv_acc': 'QC:0000052', 'qpcv_export_list': ['QC:0000038', 'QC:0000044']}
ERRORDISTR_RPLOT = {'R': "../ProduceQCFigures_acc.R", 'qp_att_acc': 'QC:0000041', 'cv_acc': 'QC:0000053', 'qpcv_export_list': ['QC:0000038']}
ERRORTIME_RPLOT = {'R': "../ProduceQCFigures_acc_time.R", 'qp_att_acc': 'QC:0000041', 'cv_acc': 'QC:0000054', 'qpcv_export_list': ['QC:0000038', 'QC:0000044']}
TIC_RPLOT = {'R': "../ProduceQCFigures_precursor_histogram.R", 'qp_att_acc': 'QC:0000023', 'cv_acc': 'MS:1000235', 'qpcv_export_list': ['QC:0000044']}


def create_qp(qcml, temp_dir, qc_metric):
    parent_cv = qc_metric['qp_att_acc']
    plot_cv = qc_metric['cv_acc']
    rscript = qc_metric['R']
    qpcv_export_list = qc_metric['qpcv_export_list']
    export_files = list()
    try:
        for export_qccv in qpcv_export_list:
            export_path = str(temp_dir) + '/' + str(export_qccv).replace(':', '_') + '.csv'
            eexe = "QCExtractor -in {inpu} -out_csv {out} -qp {qp}".format(inpu=qcml, out=export_path, qp=str(export_qccv))
            o = subprocess.check_output(eexe, stderr=subprocess.STDOUT, shell=True)
            # TODO possible to check if output is legit?!
            if o:
                export_files.append(export_path)
            else:
                raise('no ' + export_qccv)
    except Exception as e:
        logging.warning("Could not execute QCExtractor successfully for " + str(qcml) +
                        " - no QC metric " + str(export_qccv) + "(.." + str(e) + ")")
        return None

    plot_file = str(temp_dir) + '/' + str(plot_cv).replace(':', '_') + ".png"
    try:
        rexe = "Rscript {inpu} {out}".format(inpu=' '.join([rscript] + export_files), out=plot_file)
        o = subprocess.check_output(rexe, stderr=subprocess.STDOUT, shell=True)
        if not o:
            raise('no ' + plot_file)

    except Exception as e:
        logging.warning("Could not execute R " + str(rscript) + ' successfully for ' + str(qcml) +
                        " - no " + str(rscript) + " QC metric" + "(.." + str(e) + ")")
        return None

    try:
        eexe = "QCEmbedder -in {inpu} -out {out} -plot {plot} -qp_att_acc {parent} -cv_acc {type}".format(
            inpu=qcml, out=qcml, plot=plot_file, parent=parent_cv, type=plot_cv)
        o = subprocess.check_output(eexe, stderr=subprocess.STDOUT, shell=True)
        if not o:
            raise('no ' + qcml)

    except Exception as e:
        logging.warning("Could not execute QCEembedder successfully for " + str(plot_file) +
                        " - no " + str(rscript) + " QC metric" + "(.." + str(e) + ")")
        return None
    return qcml


def create_qcml(mzid, spectra, directory):
    qcml = str(directory) + '/' + os.path.splitext(os.path.basename(spectra))[0] + '.qcML'
    exe = "QCCalculator -in {inpu} -out {out} -id {id}".format(inpu=spectra, out=qcml, id=mzid)
    logging.warning("Starting to create qcml for  " + str(mzid) + ', ' + str(spectra) + " ...")
    try:
        o = subprocess.check_output(exe, stderr=subprocess.STDOUT, shell=True)
    except Exception as e:
        logging.warning("Could not create qcml for  " + str(mzid) + ', ' + str(spectra) +
                        " - skipping. (" + str(e) + ")")
    logging.warning("Created qcml for  " + str(mzid) + ', ' + str(spectra) + " - " + qcml)
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

    if options.temp_path:
        temp_dir = options.temp_path
    else:
        options.temp_dir = "/tmp/"

    logging.basicConfig(filename=options.outfiles_directory + str(datetime.now()) + '.log', filemode='w+',
                        level=logging.DEBUG)  #, format='%(levelname)s:%(message)s'
    logging.info("Starting PRIDE-QC for " + options.mzid_file + " at " + str(datetime.now()))
    args = parser.parse_args()
    logging.warning("verbosity turned on")

    #print options.mzid_file, options.spectra_list, options.outfile_directory
    # TODO go through all spectra files and create one by one the qcmls
    metrics_list = [
       IDMAP_RPLOT,
       ERRORDISTR_RPLOT,
       ERRORTIME_RPLOT,
       TIC_RPLOT,
       EXPMAP_RPLOT
    ]
    for mse in options.spectra_list:
        qcml = create_qcml(mzid=options.mzid_file, spectra=mse, directory=options.outfiles_directory)
        for metric in metrics_list:
            create_qp(qcml, temp_dir, metric)



if __name__ == '__main__':
    __main__()
