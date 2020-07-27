#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/25 4:07 PM
import os
import logging
from logging.handlers import RotatingFileHandler


def set_logging(log_name):
    u"""
    set logging handler
    modified from Zhang yiming at 2018.11.13
    :return:
    """
    sh = logging.StreamHandler()
    fh = RotatingFileHandler(
        os.path.join("{}_error.log".format(log_name)),
        mode='a+',
        encoding="utf-8",
        maxBytes=5 * 1024 * 1024,
        backupCount=2,
        delay=0
    )

    formatter = logging.Formatter(
        fmt="[%(asctime)s] - [%(levelname)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    sh.setFormatter(formatter)
    fh.setFormatter(formatter)
    sh.setLevel(logging.INFO)
    fh.setLevel(logging.DEBUG)

    log = logging.getLogger(log_name)
    log.setLevel(logging.DEBUG)
    log.addHandler(sh)
    log.addHandler(fh)
    return log

# logger = set_logging("pySplice")
