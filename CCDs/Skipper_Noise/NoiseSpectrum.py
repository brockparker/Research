#!/usr/bin/python

import xlrd
import numpy as np


def getSpectrum():

    wb = xlrd.open_workbook("STANoiseSpectrum.xls");
    sheet = wb.sheet_by_index(0);
    
    freq = [];
    noise = [];

    for row in range(sheet.nrows):
        freq.append(sheet.cell_value(row,0));
        noise.append(sheet.cell_value(row,1));

    freq = np.array(freq);
    noise = np.array(noise);

    return freq,noise
