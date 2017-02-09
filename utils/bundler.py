#!/usr/bin/python
# #### BEGIN LICENSE BLOCK ####
#
# bundler.py - Python convenience module for running Bundler.
# Copyright (C) 2013 Isaac Lenton (aka ilent2)
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# #### END LICENSE BLOCK ####

import argparse
import gzip
import os
import sys
import glob
import subprocess
import tempfile
from multiprocessing import Pool
from PIL import Image, ExifTags

VERSION = "Bundler 0.4"
DESCRIPTION = """\
Python convenience module to process a series of images and reconstruct
the scene using Bundler.

Bundler is a structure-from-motion system for unordered image
collections (for instance, images from the Internet). Bundler takes a
set of images, image features, and image matches as input, and
produces a 3D reconstruction of the camera and (sparse) scene geometry
as output."""

# This module replaces the existing RunBundler.sh script with a more
# cross platform implementation.  Additional elements replaced:
#   - RunBundler.sh             2008-2013 Noah Snavely
#   - ToSift.sh
#   - extract_focal.pl          2005-2009 Noah Snavely
#   - jhead

MOD_PATH = os.path.dirname(__file__)
BIN_PATH = os.path.join(MOD_PATH, "../bin")
LIB_PATH = os.path.join(MOD_PATH, "../lib")
BIN_SIFT = None
BIN_BUNDLER = None
BIN_MATCHKEYS = None

CCD_WIDTHS = {
     "Asahi Optical Co.,Ltd.  PENTAX Optio330RS" : 7.176, # 1/1.8"
     "Canon Canon DIGITAL IXUS 400"              : 7.176,  # 1/1.8"
     "Canon Canon DIGITAL IXUS 40"               : 5.76,   # 1/2.5"
     "Canon Canon DIGITAL IXUS 430"              : 7.176,  # 1/1.8"
     "Canon Canon DIGITAL IXUS 500"              : 7.176,  # 1/1.8"
     "Canon Canon DIGITAL IXUS 50"               : 5.76,   # 1/2.5"
     "Canon Canon DIGITAL IXUS 55"               : 5.76,   # 1/2.5"
     "Canon Canon DIGITAL IXUS 60"               : 5.76,   # 1/2.5"
     "Canon Canon DIGITAL IXUS 65"               : 5.76,   # 1/2.5"
     "Canon Canon DIGITAL IXUS 700"              : 7.176,  # 1/1.8"
     "Canon Canon DIGITAL IXUS 750"              : 7.176,  # 1/1.8"
     "Canon Canon DIGITAL IXUS 800 IS"           : 5.76,   # 1/2.5"
     "Canon Canon DIGITAL IXUS II"               : 5.27,   # 1/2.7"
     "Canon Canon EOS 10D"                       : 22.7,
     "Canon Canon EOS-1D Mark II"                : 28.7,
     "Canon Canon EOS-1Ds Mark II"               : 35.95,   
     "Canon Canon EOS  20D"                      : 22.5,
     "Canon Canon EOS 20D"                       : 22.5,
     "Canon Canon EOS 300D DIGITAL"              : 22.66,
     "Canon Canon EOS 30D"                       : 22.5,
     "Canon Canon EOS 350D DIGITAL"              : 22.2,
     "Canon Canon EOS 400D DIGITAL"              : 22.2,
     "Canon Canon EOS 40D"                       : 22.2,
     "Canon Canon EOS 5D"                        : 35.8,
     "Canon Canon EOS 5D Mark II"                : 36.0,
     "Canon Canon EOS 5D Mark III"               : 36.0,
     "Canon EOS DIGITAL REBEL"                   : 22.66,
     "Canon Canon EOS DIGITAL REBEL"             : 22.66,
     "Canon Canon EOS DIGITAL REBEL XT"          : 22.2,
     "Canon Canon EOS DIGITAL REBEL XTi"         : 22.2,
     "Canon Canon EOS REBEL T5i"                 : 22.3,
     "Canon Canon EOS Kiss Digital"              : 22.66,
     "Canon Canon IXY DIGITAL 600"               : 7.176,  # 1/1.8"
     "Canon Canon PowerShot A20"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot A400"                : 4.54,   # 1/3.2"
     "Canon Canon PowerShot A40"                 : 5.27,   # 1/2.7"
     "Canon Canon PowerShot A510"                : 5.76,   # 1/2.5"
     "Canon Canon PowerShot A520"                : 5.76,   # 1/2.5"
     "Canon Canon PowerShot A530"                : 5.76,   # 1/2.5"
     "Canon Canon PowerShot A60"                 : 5.27,   # 1/2.7"
     "Canon Canon PowerShot A620"                : 7.176,  # 1/1.8"
     "Canon Canon PowerShot A630"                : 7.176,  # 1/1.8"
     "Canon Canon PowerShot A640"                : 7.176,  # 1/1.8"
     "Canon Canon PowerShot A700"                : 5.76,   # 1/2.5"
     "Canon Canon PowerShot A70"                 : 5.27,   # 1/2.7"
     "Canon Canon PowerShot A710 IS"             : 5.76,   # 1/2.5"
     "Canon Canon PowerShot A75"                 : 5.27,   # 1/2.7"
     "Canon Canon PowerShot A80"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot A85"                 : 5.27,   # 1/2.7"
     "Canon Canon PowerShot A95"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G1"                  : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G2"                  : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G3"                  : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G5"                  : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G6"                  : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G7"                  : 7.176,  # 1/1.8"
     "Canon Canon PowerShot G9"                  : 7.600,  # 1/1.7"
     "Canon Canon PowerShot Pro1"                : 8.8,    # 2/3"
     "Canon Canon PowerShot S110"                : 5.27,   # 1/2.7"
     "Canon Canon PowerShot S1 IS"               : 5.27,   # 1/2.7"
     "Canon Canon PowerShot S200"                : 5.27,   # 1/2.7"
     "Canon Canon PowerShot S2 IS"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot S30"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S3 IS"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot S400"                : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S40"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S410"                : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S45"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S500"                : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S50"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S60"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S70"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot S80"                 : 7.176,  # 1/1.8"
     "Canon Canon PowerShot SD1000"              : 5.75,   # 1/2.5"
     "Canon Canon PowerShot SD100"               : 5.27,   # 1/2.7"
     "Canon Canon PowerShot SD10"                : 5.75,   # 1/2.5"
     "Canon Canon PowerShot SD110"               : 5.27,   # 1/2.7"
     "Canon Canon PowerShot SD200"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD300"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD400"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD450"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD500"               : 7.176,  # 1/1.8"
     "Canon Canon PowerShot SD550"               : 7.176,  # 1/1.8"
     "Canon Canon PowerShot SD600"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD630"               : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD700 IS"            : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SD750"               : 5.75,   # 1/2.5"
     "Canon Canon PowerShot SD800 IS"            : 5.76,   # 1/2.5"
     "Canon Canon PowerShot SX500 IS"            : 6.17,   # 1/2.3"
     "Canon EOS 300D DIGITAL"                    : 22.66,
     "Canon PowerShot A510"                      : 5.76,   # 1/2.5" ???
     "Canon PowerShot S30"                       : 7.176,  # 1/1.8"
     "CASIO COMPUTER CO.,LTD. EX-S500"           : 5.76,   # 1/2.5"
     "CASIO COMPUTER CO.,LTD. EX-Z1000"          : 7.716, # 1/1.8"
     "CASIO COMPUTER CO.,LTD  EX-Z30"            : 5.76,   # 1/2.5 "
     "CASIO COMPUTER CO.,LTD. EX-Z600"           : 5.76,   # 1/2.5"
     "CASIO COMPUTER CO.,LTD. EX-Z60"            : 7.176, # 1/1.8"
     "CASIO COMPUTER CO.,LTD  EX-Z750"           : 7.176, # 1/1.8"
     "CASIO COMPUTER CO.,LTD. EX-Z850"           : 7.176,
     "EASTMAN KODAK COMPANY KODAK CX7330 ZOOM DIGITAL CAMERA" : 5.27, # 1/2.7"
     "EASTMAN KODAK COMPANY KODAK CX7530 ZOOM DIGITAL CAMERA" : 5.76, # 1/2.5"
     "EASTMAN KODAK COMPANY KODAK DX3900 ZOOM DIGITAL CAMERA" : 7.176, # 1/1.8"
     "EASTMAN KODAK COMPANY KODAK DX4900 ZOOM DIGITAL CAMERA" : 7.176, # 1/1.8"
     "EASTMAN KODAK COMPANY KODAK DX6340 ZOOM DIGITAL CAMERA" : 5.27, # 1/2.7"
     "EASTMAN KODAK COMPANY KODAK DX6490 ZOOM DIGITAL CAMERA" : 5.76, # 1/2.5"
     "EASTMAN KODAK COMPANY KODAK DX7630 ZOOM DIGITAL CAMERA" : 7.176, # 1/1.8"
     "EASTMAN KODAK COMPANY KODAK Z650 ZOOM DIGITAL CAMERA"   : 5.76, # 1/2.5"
     "EASTMAN KODAK COMPANY KODAK Z700 ZOOM DIGITAL CAMERA"   : 5.76, # 1/2.5"
     "EASTMAN KODAK COMPANY KODAK Z740 ZOOM DIGITAL CAMERA"   : 5.76, # 1/2.5"
     "EASTMAN KODAK COMPANY KODAK Z740 ZOOM DIGITAL CAMERA"   : 5.76, # 1/2.5"?
     "FUJIFILM FinePix2600Zoom"                  : 5.27,   # 1/2.7"
     "FUJIFILM FinePix40i"                       : 7.600,  # 1/1.7" 
     "FUJIFILM FinePix A310"                     : 5.27,   # 1/2.7"
     "FUJIFILM FinePix A330"                     : 5.27,   # 1/2.7"
     "FUJIFILM FinePix A600"                     : 7.600,  # 1/1.7"
     "FUJIFILM FinePix E500"                     : 5.76,   # 1/2.5" 
     "FUJIFILM FinePix E510"                     : 5.76,   # 1/2.5"
     "FUJIFILM FinePix E550"                     : 7.600,  # 1/1.7" 
     "FUJIFILM FinePix E900"                     : 7.78,   # 1/1.6"
     "FUJIFILM FinePix F10"                      : 7.600, # 1/1.7"
     "FUJIFILM FinePix F30"                      : 7.600,  # 1/1.7"
     "FUJIFILM FinePix F450"                     : 5.76,   # 1/2.5"
     "FUJIFILM FinePix F601 ZOOM"                : 7.600,  # 1/1.7"
     "FUJIFILM FinePix S3Pro"                    : 23.0,
     "FUJIFILM FinePix S5000"                    : 5.27,   # 1/2.7"
     "FUJIFILM FinePix S5200"                    : 5.76,   # 1/2.5"
     "FUJIFILM FinePix S5500"                    : 5.27,   # 1/2.7"
     "FUJIFILM FinePix S6500fd"                  : 7.600,  # 1/1.7"
     "FUJIFILM FinePix S7000"                    : 7.600,  # 1/1.7"
     "FUJIFILM FinePix Z2"                       : 5.76,   # 1/2.5"
     "Hewlett-Packard hp 635 Digital Camera"     : 4.54, # 1/3.2"
     "Hewlett-Packard hp PhotoSmart 43x series"  : 5.27,  # 1/2.7"
     "Hewlett-Packard HP PhotoSmart 618 (V1.1)"  : 5.27,  # 1/2.7"
     "Hewlett-Packard HP PhotoSmart C945 (V01.61)" : 7.176, # 1/1.8"
     "Hewlett-Packard HP PhotoSmart R707 (V01.00)" : 7.176, # 1/1.8"
     "KONICA MILOLTA  DYNAX 5D"                  : 23.5,
     "Konica Minolta Camera, Inc. DiMAGE A2"     : 8.80, # 2/3"
     "KONICA MINOLTA CAMERA, Inc. DiMAGE G400"   : 5.76, # 1/2.5"
     "Konica Minolta Camera, Inc. DiMAGE Z2"     : 5.76, # 1/2.5"
     "KONICA MINOLTA DiMAGE A200"                : 8.80,   # 2/3"
     "KONICA MINOLTA DiMAGE X1"                  : 7.176,  # 1/1.8"
     "KONICA MINOLTA  DYNAX 5D"                  : 23.5,
     "Minolta Co., Ltd. DiMAGE F100"             : 7.176,  # 1/2.7"
     "Minolta Co., Ltd. DiMAGE Xi"               : 5.27,   # 1/2.7"
     "Minolta Co., Ltd. DiMAGE Xt"               : 5.27,   # 1/2.7"
     "Minolta Co., Ltd. DiMAGE Z1"               : 5.27, # 1/2.7" 
     "NIKON COOLPIX L3"                          : 5.76,   # 1/2.5"
     "NIKON COOLPIX P2"                          : 7.176,  # 1/1.8"
     "NIKON COOLPIX S4"                          : 5.76,   # 1/2.5"
     "NIKON COOLPIX S7c"                         : 5.76,   # 1/2.5"
     "NIKON CORPORATION NIKON D100"              : 23.7,
     "NIKON CORPORATION NIKON D1"                : 23.7,
     "NIKON CORPORATION NIKON D1H"               : 23.7,
     "NIKON CORPORATION NIKON D200"              : 23.6,
     "NIKON CORPORATION NIKON D2H"               : 23.3,
     "NIKON CORPORATION NIKON D2X"               : 23.7,
     "NIKON CORPORATION NIKON D40"               : 23.7,
     "NIKON CORPORATION NIKON D50"               : 23.7,
     "NIKON CORPORATION NIKON D60"               : 23.6,
     "NIKON CORPORATION NIKON D70"               : 23.7,
     "NIKON CORPORATION NIKON D70s"              : 23.7,
     "NIKON CORPORATION NIKON D80"               : 23.6,
     "NIKON E2500"                               : 5.27,   # 1/2.7"
     "NIKON E2500"                               : 5.27,   # 1/2.7"
     "NIKON E3100"                               : 5.27,   # 1/2.7"
     "NIKON E3200"                               : 5.27,
     "NIKON E3700"                               : 5.27,   # 1/2.7"
     "NIKON E4200"                               : 7.176,  # 1/1.8"
     "NIKON E4300"                               : 7.18,
     "NIKON E4500"                               : 7.176,  # 1/1.8"
     "NIKON E4600"                               : 5.76,   # 1/2.5"
     "NIKON E5000"                               : 8.80,   # 2/3"
     "NIKON E5200"                               : 7.176,  # 1/1.8"
     "NIKON E5400"                               : 7.176,  # 1/1.8"
     "NIKON E5600"                               : 5.76,   # 1/2.5"
     "NIKON E5700"                               : 8.80,   # 2/3"
     "NIKON E5900"                               : 7.176,  # 1/1.8"
     "NIKON E7600"                               : 7.176,  # 1/1.8"
     "NIKON E775"                                : 5.27,   # 1/2.7"
     "NIKON E7900"                               : 7.176,  # 1/1.8"
     "NIKON E7900"                               : 7.176,  # 1/1.8"
     "NIKON E8800"                               : 8.80,   # 2/3"
     "NIKON E990"                                : 7.176,  # 1/1.8"
     "NIKON E995"                                : 7.176,  # 1/1.8"
     "NIKON S1"                                  : 5.76,   # 1/2.5"
     "Nokia N80"                                 : 5.27,   # 1/2.7"
     "Nokia N80"                                 : 5.27,   # 1/2.7"
     "Nokia N93"                                 : 4.536,  # 1/3.1"
     "Nokia N95"                                 : 5.7,    # 1/2.7"
     "OLYMPUS CORPORATION     C-5000Z"           : 7.176,  # 1/1.8"
     "OLYMPUS CORPORATION C5060WZ"               : 7.176, # 1/1.8"
     "OLYMPUS CORPORATION C750UZ"                : 5.27,   # 1/2.7"
     "OLYMPUS CORPORATION C765UZ"                : 5.76,   # 1//2.5"
     "OLYMPUS CORPORATION C8080WZ"               : 8.80,   # 2/3"
     "OLYMPUS CORPORATION X250,D560Z,C350Z"      : 5.76, # 1/2.5" 
     "OLYMPUS CORPORATION     X-3,C-60Z"         : 7.176, # 1.8"
     "OLYMPUS CORPORATION X400,D580Z,C460Z"      : 5.27,  # 1/2.7"
     "OLYMPUS IMAGING CORP.   E-500"             : 17.3,  # 4/3?
     "OLYMPUS IMAGING CORP.   FE115,X715"        : 5.76, # 1/2.5"
     "OLYMPUS IMAGING CORP. SP310"               : 7.176, # 1/1.8"
     "OLYMPUS IMAGING CORP.   SP510UZ"           : 5.75,   # 1/2.5"
     "OLYMPUS IMAGING CORP.   SP550UZ"           : 5.76, # 1/2.5"
     "OLYMPUS IMAGING CORP.   uD600,S600"        : 5.75, # 1/2.5" 
     "OLYMPUS_IMAGING_CORP.   X450,D535Z,C370Z"  : 5.27, # 1/2.7" 
     "OLYMPUS IMAGING CORP. X550,D545Z,C480Z"    : 5.76, # 1/2.5" 
     "OLYMPUS OPTICAL CO.,LTD C2040Z"            : 6.40,  # 1/2"
     "OLYMPUS OPTICAL CO.,LTD C211Z"             : 5.27,   # 1/2.7"
     "OLYMPUS OPTICAL CO.,LTD C2Z,D520Z,C220Z"   : 4.54, # 1/3.2"
     "OLYMPUS OPTICAL CO.,LTD C3000Z"            : 7.176, # 1/1.8"
     "OLYMPUS OPTICAL CO.,LTD C300Z,D550Z"       : 5.4,
     "OLYMPUS OPTICAL CO.,LTD C4100Z,C4000Z"     : 7.176,  # 1/1.8" 
     "OLYMPUS OPTICAL CO.,LTD C750UZ"            : 5.27,  # 1/2.7"
     "OLYMPUS OPTICAL CO.,LTD X-2,C-50Z"         : 7.176, # 1/1.8"
     "OLYMPUS SP550UZ"                           : 5.76,  # 1/2.5"
     "OLYMPUS X100,D540Z,C310Z"                  : 5.27,   # 1/2.7"
     "Panasonic DMC-FX01"                        : 5.76,   # 1/2.5"
     "Panasonic DMC-FX07"                        : 5.75,   # 1/2.5"
     "Panasonic DMC-FX9"                         : 5.76,   # 1/2.5"
     "Panasonic DMC-FZ20"                        : 5.760,  # 1/2.5"
     "Panasonic DMC-FZ2"                         : 4.54,   # 1/3.2"
     "Panasonic DMC-FZ30"                        : 7.176,  # 1/1.8"
     "Panasonic DMC-FZ50"                        : 7.176,  # 1/1.8"
     "Panasonic DMC-FZ5"                         : 5.760,  # 1/2.5"
     "Panasonic DMC-FZ7"                         : 5.76,   # 1/2.5"
     "Panasonic DMC-LC1"                         : 8.80,   # 2/3"
     "Panasonic DMC-LC33"                        : 5.760,  # 1/2.5"
     "Panasonic DMC-LX1"                         : 8.50,   # 1/6.5"
     "Panasonic DMC-LZ2"                         : 5.76,   # 1/2.5"
     "Panasonic DMC-TZ1"                         : 5.75,   # 1/2.5"
     "Panasonic DMC-TZ3"                         : 5.68,   # 1/2.35"
     "PENTAX Corporation  PENTAX *ist DL"        : 23.5,
     "PENTAX Corporation  PENTAX *ist DS2"       : 23.5,
     "PENTAX Corporation  PENTAX *ist DS"        : 23.5,
     "PENTAX Corporation  PENTAX K100D"          : 23.5,
     "PENTAX Corporation PENTAX Optio 450"       : 7.176, # 1/1.8"
     "PENTAX Corporation PENTAX Optio 550"       : 7.176, # 1/1.8"
     "PENTAX Corporation PENTAX Optio E10"       : 5.76, # 1/2.5"
     "PENTAX Corporation PENTAX Optio S40"       : 5.76, # 1/2.5"
     "PENTAX Corporation  PENTAX Optio S4"       : 5.76, # 1/2.5" 
     "PENTAX Corporation PENTAX Optio S50"       : 5.76, # 1/2.5"
     "PENTAX Corporation  PENTAX Optio S5i"      : 5.76, # 1/2.5" 
     "PENTAX Corporation  PENTAX Optio S5z"      : 5.76, # 1/2.5" 
     "PENTAX Corporation  PENTAX Optio SV"       : 5.76, # 1/2.5" 
     "PENTAX Corporation PENTAX Optio WP"        : 5.75, # 1/2.5" 
     "RICOH CaplioG3 modelM"                     : 5.27,   # 1/2.7"
     "RICOH       Caplio GX"                     : 7.176,  # 1/1.8"
     "RICOH       Caplio R30"                    : 5.75,   # 1/2.5"
     "Samsung  Digimax 301"                      : 5.27,   # 1/2.7"
     "Samsung Techwin <Digimax i5, Samsung #1>"  : 5.76,   # 1/2.5"
     "SAMSUNG TECHWIN Pro 815"                   : 8.80,   # 2/3"
     "SONY DSC-F828"                             : 8.80,   # 2/3"
     "SONY DSC-N12"                              : 7.176,  # 1/1.8"
     "SONY DSC-P100"                             : 7.176,  # 1/1.8"
     "SONY DSC-P10"                              : 7.176,  # 1/1.8"
     "SONY DSC-P12"                              : 7.176,  # 1/1.8"
     "SONY DSC-P150"                             : 7.176,  # 1/1.8"
     "SONY DSC-P200"                             : 7.176,   # 1/1.8");
     "SONY DSC-P52"                              : 5.27,   # 1/2.7"
     "SONY DSC-P72"                              : 5.27,   # 1/2.7"
     "SONY DSC-P73"                              : 5.27,
     "SONY DSC-P8"                               : 5.27,   # 1/2.7"
     "SONY DSC-R1"                               : 21.5,
     "SONY DSC-S40"                              : 5.27,   # 1/2.7"
     "SONY DSC-S600"                             : 5.760,  # 1/2.5"
     "SONY DSC-T9"                               : 7.18,
     "SONY DSC-V1"                               : 7.176,  # 1/1.8"
     "SONY DSC-W1"                               : 7.176,  # 1/1.8"
     "SONY DSC-W30"                              : 5.760,  # 1/2.5"
     "SONY DSC-W50"                              : 5.75,   # 1/2.5"
     "SONY DSC-W5"                               : 7.176,  # 1/1.8"
     "SONY DSC-W7"                               : 7.176,  # 1/1.8"
     "SONY DSC-W80"                              : 5.75,   # 1/2.5"
     "SONY DSLR-A700"                            : 23.5,   # 23.5 x 15.6mm
     "SONY SLT-A99V"                             : 35.8,   # 35.8 x 23.9mm
}

def get_images():
    """Searches the present directory for JPEG images."""
    images = glob.glob("./*.[jJ][pP][gG]")
    if len(images) == 0:
        error_str = ("Error: No images supplied!  "
                     "No JPEG files found in directory!")
        raise Exception(error_str)
    return images

def extract_focal_length(images=[], scale=1.0, verbose=False):
    """Extracts (pixel) focal length from images where available.
    The functions returns a dictionary of image, focal length pairs.
    If no focal length is extracted for an image, the second pair is None.
    """
    if len(images) == 0:
        if verbose: print "[- Creating list of images -]"
        images = get_images()

    ret = {}
    for image in images:
        if verbose: print "[Extracting EXIF tags from image {0}]".format(image)

        tags = {}
        with open(image, 'rb') as fp:
            img = Image.open(fp)
            if hasattr(img, '_getexif'):
                exifinfo = img._getexif()
                if exifinfo is not None:
                    for tag, value in exifinfo.items():
                        tags[ExifTags.TAGS.get(tag, tag)] = value

        ret[image] = None

        # Extract Focal Length
        focalN, focalD = tags.get('FocalLength', (0, 1))
        focal_length = float(focalN)/float(focalD)

        # Extract Resolution
        img_width = tags.get('ExifImageWidth', 0)
        img_height = tags.get('ExifImageHeight', 0)
        if img_width < img_height:
            img_width,img_height = img_height,img_width



        # Extract CCD Width (Prefer Lookup Table)
        ccd_width = 1.0
        make_model = tags.get('Make', '') + ' ' + tags.get('Model', '')
        if CCD_WIDTHS.has_key(make_model.strip()):
            ccd_width = CCD_WIDTHS[make_model.strip()]
        else:
            fplaneN, fplaneD = tags.get('FocalPlaneXResolution', (0, 1))
            if fplaneN != 0:
                ccd_width = 25.4*float(img_width)*float(fplaneD)/float(fplaneN)
                if verbose: print "  [Using CCD width from EXIF tags]"
            else:
                ccd_width = 0

        if verbose:
            print "  [EXIF focal length = {0}mm]".format(focal_length)
            print "  [EXIF CCD width = {0}mm]".format(ccd_width)
            print "  [EXIF resolution = {0} x {1}]".format(
                img_width, img_height)
            if ccd_width == 0:
                print "  [No CCD width available for camera {0}]".format(
                    make_model)

        if (img_width==0 or img_height==0 or focalN==0 or ccd_width==0):
            if verbose: print "  [Could not determine pixel focal length]"
            continue

        # Compute Focal Length in Pixels
        ret[image] = img_width * (focal_length / ccd_width) * scale
        if verbose: print "  [Focal length (pixels) = {0}]".format(ret[image])

    return ret

def sift_image(image, verbose=False):
    """Extracts SIFT features from a single image.  See sift_images."""
    global BIN_SIFT, BIN_PATH

    if BIN_SIFT is None:
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            BIN_SIFT = os.path.join(BIN_PATH, "siftWin32.exe")
        else:
            BIN_SIFT = os.path.join(BIN_PATH, "sift")

    pgm_filename = image.rsplit('.', 1)[0] + ".pgm"
    key_filename = image.rsplit('.', 1)[0] + ".key"

    # Convert image to PGM format (grayscale)
    with open(image, 'rb') as fp_img:
        image = Image.open(fp_img)
        image.convert('L').save(pgm_filename)

    # Extract SIFT data
    if verbose:
        with open(pgm_filename, 'rb') as fp_in:
            with open(key_filename, 'wb') as fp_out:
                subprocess.call(BIN_SIFT, stdin=fp_in, stdout=fp_out)
    else:
        with open(pgm_filename, 'rb') as fp_in:
            with open(key_filename, 'wb') as fp_out:
                with open(os.devnull, 'w') as fp_err:
                    subprocess.call(BIN_SIFT, stdin=fp_in, stdout=fp_out,
                                    stderr=fp_err)

    # Remove pgm file
    os.remove(pgm_filename)

    # GZIP compress key file (and remove)
    with open(key_filename, 'rb') as fp_in:
        with gzip.open(key_filename + ".gz", 'wb') as fp_out:
            fp_out.writelines(fp_in)
    os.remove(key_filename)

    return key_filename

def sift_images(images, verbose=False, parallel=True):
    """Extracts SIFT features from images in 'images'.

    'images' should be a list of file names.  The function creates a
    SIFT compressed key file for each image in 'images' with a '.key.gz'
    extension.  A list of the uncompressed key file names is returned.

    If 'parallel' is True, the function executes SIFT in parallel.
    """
    global BIN_SIFT, BIN_PATH

    key_filenames = []

    if BIN_SIFT is None:
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            BIN_SIFT = os.path.join(BIN_PATH, "siftWin32.exe")
        else:
            BIN_SIFT = os.path.join(BIN_PATH, "sift")
        
    if parallel:
        pool = Pool()
        key_filenames = pool.map(sift_image, images)
    else:
        for image in images:
            key_filenames.append(sift_image(image, verbose=verbose))

    return key_filenames

def match_images(key_files, matches_file, verbose=False):
    "Executes KeyMatchFull to match key points in images."""
    global BIN_MATCHKEYS, BIN_PATH

    if BIN_MATCHKEYS is None:
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            BIN_MATCHKEYS = os.path.join(BIN_PATH, "KeyMatchFull.exe")
        else:
            BIN_MATCHKEYS = os.path.join(BIN_PATH, "KeyMatchFull")

    keys_file = ""
    with tempfile.NamedTemporaryFile(delete=False) as fp:
        for key in key_files:
            fp.write(key + '\n')
        keys_file = fp.name

    # Add lib folder to LD_LIBRARY_PATH
    env = dict(os.environ)
    if env.has_key('LD_LIBRARY_PATH'):
        env['LD_LIBRARY_PATH'] = env['LD_LIBRARY_PATH'] + ':' + LIB_PATH
    else:
        env['LD_LIBRARY_PATH'] = LIB_PATH

    if verbose:
        subprocess.call([BIN_MATCHKEYS, keys_file, matches_file], env=env)
    else:
        with open(os.devnull, 'w') as fp_out:
            subprocess.call([BIN_MATCHKEYS, keys_file, matches_file],
                            stdout=fp_out, env=env)
            
    os.remove(keys_file)

def bundler(image_list=None, options_file=None, shell=False, *args, **kwargs):
    """Run bundler, parsing arguments from args and kwargs through.
    For Bundler usage run bundler("--help").

    image_list : File containing list of images.
    options_file : Specify an options file for bundler (optional).
    shell : Enable full shell support for parsing args (default: False).
    """
    global BIN_BUNDLER, BIN_PATH

    if BIN_BUNDLER is None:
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            BIN_BUNDLER = os.path.join(BIN_PATH, "Bundler.exe")
        else:
            BIN_BUNDLER = os.path.join(BIN_PATH, "bundler")

    def kwargs_bool(b, r):
        if b: return r
        else: return []

    kwargs_dict = {
        'match_table'            : lambda k,v: ['--'+k,v],
        'output'                 : lambda k,v: ['--'+k,v],
        'output_all'             : lambda k,v: ['--'+k,v],
        'output_dir'             : lambda k,v: ['--'+k,v],
        'variable_focal_length'  : lambda k,v: kwargs_bool(v, ['--'+k]),
        'use_focal_estimate'     : lambda k,v: kwargs_bool(v, ['--'+k]),
        'constrain_focal'        : lambda k,v: kwargs_bool(v, ['--'+k]),
        'constrain_focal_weight' : lambda k,v: ['--'+k,str(v)],
        'estimate_distortion'    : lambda k,v: kwargs_bool(v, ['--'+k]),
        'run_bundle'             : lambda k,v: kwargs_bool(v, ['--'+k]),
    }

    str_args = [a for a in args if type(a) == str]
    for k,v in kwargs.items():
        if not kwargs_dict.has_key(k): continue
        str_args.extend(kwargs_dict[k](k,v))

    if len(str_args) != 0 and options_file is not None:
        with open(options_file, 'wb') as fp:
            for o in str_args:
                if o.startswith('--'): fp.write('\n')
                else: fp.write(' ')
                fp.write(o)

    image_list_file = ""
    if type(image_list) == dict:
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            for image,value in image_list.items():
                if value == None: fp.write(image + '\n')
                else: fp.write(' '.join([image, '0', str(value), '\n']))
            image_list_file = fp.name
    elif type(image_list) == str:
        image_list_file = image_list
    else:
        raise Exception("Error: Not a valid list or filename for image_list!")

    # Add lib folder to LD_LIBRARY_PATH
    env = dict(os.environ)
    if env.has_key('LD_LIBRARY_PATH'):
        env['LD_LIBRARY_PATH'] = env['LD_LIBRARY_PATH'] + ':' + LIB_PATH
    else:
        env['LD_LIBRARY_PATH'] = LIB_PATH

    try:    os.mkdir("bundle")
    except: pass

    with open(os.path.join("bundle", "out"), 'wb') as fp_out:
        if options_file is not None:
            subprocess.call([BIN_BUNDLER, image_list_file, "--options_file",
                options_file], shell=shell, env=env, stdout=fp_out)
        else:
            subprocess.call([BIN_BUNDLER, image_list_file] + str_args,
                shell=shell, env=env, stdout=fp_out)

    if type(image_list) == dict:
        os.remove(image_list_file)

def run_bundler(images=[], verbose=False, parallel=True):
    """Prepare images and run bundler with default options."""
    # Create list of images
    if len(images) == 0:
        if verbose: print "[- Creating list of images -]"
        images = get_images()

    # Extract focal length
    if type(images) == list:
        if verbose: print "[- Extracting EXIF tags from images -]"
        images = extract_focal_length(images, verbose=verbose)

    # Extract SIFT features from images
    if verbose: print "[- Extracting keypoints -]"
    key_files = sift_images(images, parallel=parallel, verbose=verbose)

    # Match images
    if verbose: print "[- Matching keypoints (this can take a while) -]"
    matches_file = "matches.init.txt"
    match_images(key_files, matches_file, verbose=verbose)

    # Run Bundler
    if verbose: print "[- Running Bundler -]"
    bundler(image_list=images,
            options_file="options.txt",
            verbose=verbose,
            match_table=matches_file,
            output="bundle.out",
            output_all="bundle_",
            output_dir="bundle",
            variable_focal_length=True,
            use_focal_estimate=True,
            constrain_focal=True,
            constrain_focal_weight=0.0001,
            estimate_distortion=True,
            run_bundle=True)

    if verbose: print "[- Done -]"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('-v', '--verbose', action='store_true',
        help="generate additional output in execution", default=False)
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--no-parallel', action='store_true',
        help="disable parallelisation", default=False)
    parser.add_argument('--extract-focal', action='store_true',
        help="only create list of images to be reconstructed", default=False)
    args = parser.parse_args()

    if args.extract_focal:
        images = extract_focal_length(verbose=args.verbose)
        with open("list.txt", 'w') as fp:
            for image,value in images.items():
                if value == None: fp.write(image + '\n')
                else: fp.write(' '.join([image, '0', str(value), '\n']))
    else:
        run_bundler(verbose=args.verbose, parallel=not args.no_parallel)

