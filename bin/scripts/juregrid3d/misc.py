#
# Copyright 2014 by Forschungszentrum Juelich GmbH
# Author: J. Ungermann
#

import logging


def setup_logging(logfile="juregrid3d.log"):
    logger = logging.getLogger('juregrid3d')
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    cf = logging.Formatter('%(message)s')
    ch.setFormatter(cf)
    logger.addHandler(ch)
    if logfile:
        fh = logging.FileHandler(logfile)
        ff = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(ff)
        fh.setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
        logger.addHandler(fh)
