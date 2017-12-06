#!/usr/bin/env python

import string, re, sys
from subprocess import *
#import Bio.Emboss.Primer3 as Primer3
#from Bio.Emboss.Applications import Primer3Commandline


class Primer:
	#PRIMER_TM_FORMULA=1 - refers to Santa Lucia JR (1998)
	#PRIMER_TM_FORMULA=0 - refers to Berlauer KJ (1986)
	#PRIMER_SALT_CORRECTIONS = 1 - refers to Santa Lucia JR (1998)
	#PRIMER_SALT_CORRECTIONS = 0 - refers to Schildkraut C (1965)
	#PRIMER_SALT_CORRECTIONS = 2 - refers to Owczarzy R (2008)
	def __init__(self,SEQUENCE_ID):
		self.SEQUENCE_ID="SEQUENCE_ID="+SEQUENCE_ID
		self.SEQUENCE_TEMPLATE="SEQUENCE_TEMPLATE="
		self.GENERIC_PARAMS="""PRIMER_TASK=pick_pcr_primers
PRIMER_TM_FORMULA=1
PRIMER_SALT_CORRECTIONS=1
PRIMER_OPT_SIZE=25
PRIMER_MAX_SIZE=28
P3_FILE_FLAG=1
PRIMER_GC_CLAMP=1
PRIMER_EXPLAIN_FLAG=1
PRIMER_OPT_GC=50
PRIMER_OPT_TM=60
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_MISPRIMING_LIBRARY=/home/clarkg/INSTALLS/primer3-2.4.0/src/humrep_and_simple.txt
"""
		self.ENDFILE="="
		self.PRIMER_PRODUCT_SIZE_RANGE="PRIMER_PRODUCT_SIZE_RANGE=150-250 251-500 501-750 751-1000 1001-1200\n"
		self.PRIMER_MIN_GC="PRIMER_MIN_GC=45\n"
		self.PRIMER_MAX_GC="PRIMER_MAX_GC=54\n"
		self.PRIMER_MIN_SIZE="PRIMER_MIN_SIZE=22\n"
		self.HAIRPIN="PRIMER_MAX_HAIRPIN_TH=10\n"
		self.SELF_TH="PRIMER_MAX_SELF_ANY_TH=10\n"
		self.END_TH="PRIMER_MAX_SELF_END_TH=10\n"
		self.PRIMER_MIN_TM="PRIMER_MIN_TM=58\n"
		self.PRIMER_MAX_TM="PRIMER_MAX_TM=64\n"
		self.MASKING="PRIMER_LOWERCASE_MASKING=1\n"

	def updateParams(self,attemptNum):
		if attemptNum == 1:
			self.MASKING="PRIMER_LOWERCASE_MASKING=0\n"
		elif attemptNum == 2:
			##Have already unset masking
			self.PRIMER_MIN_SIZE="PRIMER_MIN_SIZE=20\n"
		elif attemptNum == 3:
			##Have already unset masking
			##Have already changed min primer size
			self.PRIMER_MIN_GC="PRIMER_MIN_GC=40\n"
		elif attemptNum == 4:
			##Have already unset masking
			##Have already changed min primer size
			##Have already changed min gc
			self.PRIMER_MIN_GC="PRIMER_MAX_GC=58\n"
		elif attemptNum == 5:
			##Now we're getting desperate
			self.PRIMER_MIN_TM="PRIMER_MIN_TM=56\n"
			self.PRIMER_MAX_TM="PRIMER_MAX_TM=66\n"

