#!/usr/bin/env python

import string, re, sys
import glob,os
import subprocess
from primerDesign import Primer
from subprocess import Popen, PIPE
import operator
from collections import defaultdict
import getpass
import datetime
import numpy as np


class PrimerPairing:

	def __init__(self):
		self.For_ENDTH=""	
		self.For_ANYTH=""	
		self.For_SEQUENCE=""	
		self.For_GC=""
		self.For_TM=""	
		self.For_LENGTH=""	
		self.For_START=""	
		self.For_HAIRPINTH=""
	
		self.Rev_ENDTH=""	
		self.Rev_ANYTH=""	
		self.Rev_SEQUENCE=""	
		self.Rev_GC=""	
		self.Rev_TM=""	
		self.Rev_LENGTH=""	
		self.Rev_START=""	
		self.Rev_HAIRPINTH=""	

		self.PairANYTH=""
		self.PairENDTH=""
		self.Complete=False

	def CheckComplete(self):
		varia=[getattr(self,attr) for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__") and attr != "Complete"]
		filled_variables=filter(lambda p: len(str(p)),varia)
		##We want to check if all the class variables are populated 
		if len(filled_variables) >= 18:		##could just be == 18
			self.Complete=True
class APEfile:
	if not glob.glob("/mnt/cmmr/"):
		print "Must mount /mnt/cmmr/ using 'rootmount.sh' in home directory. Ensure that passwords for user have no changed"
		sys.exit()
	misc=re.compile("[\W]+misc_feature[\W]+")
	locus=re.compile("^[\W]+/locus_tag=")
	label=re.compile("^[\W]+/label=")
	exonlabel=re.compile("^[\W]+exon[\W]+[(]{0,1}[0-9]{1,7}..[0-9]{1,7}")
	exontag=re.compile("^[\W]+exon[\W]+[(]{0,1}[0-9]{1,7}..[0-9]{1,7}")
	exonlabel=re.compile("^[\W]+/label=ENSMUSE")
	def __init__(self,symbol):
		self.symbol=symbol
		self.mountpoint="/mnt/cmmr/Cas9"
		self.fileDIR=[]
		self.possibleDIR=[]
		self.misc=re.compile("[\W]+misc_feature[\W]+")
		self.coordinates=re.compile("[0-9]{1,7}..[0-9]{1,7}[)]{0,1}")
		self.fileFound=False
		self.feature_desc=""
		self.Attempt=0
		self.AttemptType={'EM':0,'WT':0}
		self.ValidPairs=False
		self.pairs={}
		self.FOR={}
		self.REV={}
		self.ChromFiles=glob.glob("/home/clarkg/Chromosomes/*.2bit")
		self.ChromFiles.sort()
		self.pcrdir="/home/clarkg/PrimerDesign/PCRdir"
		self.PcrExec='/home/clarkg/PrimerDesign/isPcr'
		self.RunDir='/home/clarkg/PrimerDesign/RunData'
		self.PrimerPairs={}
		self.PrimerTypedPairs={}
		self.VerifiedSequencePairs={}
		self.RankVerifiedPairs={}
		self.Primer3Output=""
		self.Primer3Error=""

	def findFile(self):
		self.possibleDIR+=glob.glob(os.path.join(self.mountpoint,"_IMPC_Mice",self.symbol+"*"))
		self.possibleDIR+=glob.glob(os.path.join(self.mountpoint,"_External",self.symbol+"*"))
		if len(self.possibleDIR) > 1:
			print "We have found %s possible directories for %s. They are listed below:" % (len(self.possibleDIR),self.symbol)
			for x,k in enumerate(self.possibleDIR):
				print "\t"+str(x+1)+". "+k
			print "\n"
		elif len(self.possibleDIR) ==0:
			print "No directory found for %s. Please check spelling." % (self.symbol)
			sys.exit()
		else:
			self.fileDIR=self.possibleDIR[0]
			self.APEtarget=glob.glob(os.path.join(self.fileDIR,"*target*.ape")) 
			self.APEgDNA=glob.glob(os.path.join(self.fileDIR,"*gDNA_masked.ape")) 
			if len(self.APEtarget) != 1:
				pass
				#print "Can't find target file, not fatal"
			else:
				self.APEtarget=self.APEtarget[0]	

			if len(self.APEgDNA) != 1:
				print "Can't find gDNA_masked file"
				print len(self.APEgDNA)
				sys.exit()
			else:
				self.APEgDNA=self.APEgDNA[0]
				self.fullfile=open(self.APEgDNA,'r').readlines()
				self.fileFound=True
				
	def readAPESequence(self):
		self.OriginStart=False
		fullfile=self.fullfile
		self.Sequence=""
		for x,line in enumerate(fullfile):
			if line.startswith("ORIGIN"):
				self.OriginStart=True
			if self.OriginStart:
				seqline=line.split()[1:]
				if len(seqline):
					self.Sequence+="".join(seqline)


	def RevSeq(self,sequence):
		sequence.upper()
		#print sequence
		pairs={'A':'T' , 'C':'G', 'G':'C' , 'T':'A', 'N':'N'}
		lclS=sequence[::-1]
		#lclS.upper()
		lclU=lclS.upper()
		#print lclU
		rev=[]
		for s in lclU:
			rev.append(pairs[s])
		return "".join(rev)

	def readAPEFeatures(self):
		#self.OriginStart=False
		fullfile=open(self.APEgDNA,'r').readlines()
		self.features={}
		#self.Sequence=""
		for x,line in enumerate(fullfile):
			if self.misc.search(line):
				new_feature=line.strip().split()
				if len(self.coordinates.findall(line)) == 1:
					coordinfo=self.coordinates.findall(line)[0]
					if coordinfo.endswith(")"):
						feature_name=coordinfo[:-1]+"_rev"
					else:
						feature_name=coordinfo
					lcl_feature=APEfeatures(feature_name)
					lcl_feature.grabFeatures(fullfile[x+1:])
					self.features[lcl_feature.feature_desc]=lcl_feature
			elif self.exontag.match(line):
				new_feature=line.strip().split()
				if len(self.coordinates.findall(line)) == 1:
					coordinfo=self.coordinates.findall(line)[0]
					if coordinfo.endswith(")"):
						feature_name=coordinfo[:-1]+"_rev"
					else:
						feature_name=coordinfo
					lcl_feature=APEfeatures(feature_name)
					lcl_feature.grabFeatures(fullfile[x+1:])
					self.features[lcl_feature.feature_desc]=lcl_feature
				else:
					print "COORDINATES In the APE file are too confusing\nExamine APE file and restart"
					print line
					sys.exit()
#	def Primer3Select(self):
#		if not len(self.Primer3Output):
#			print "NO OUTPUT"
#		else:
#			for line in self.Primer3Output.split("\n"):
#
#				print line
#		sys.exit()

	def PrimerSelect(self):

		##We want to capture the largest potential cutsite
		##	APE files have a variety of ways of defining gRNA sites (via the user)
		##	So we are necessarily handcuffed here at defining these sites
		try:
			self.upstream2=self.features['gRNA_U3']
			self.upstream=self.features['gRNA_U5']
			self.upstreamcut=self.features['gRNA_U5'].cutsite
			self.multiUP=True
		except KeyError:
			try:
				self.upstream=self.features['gRNA_U']
				self.upstreamcut=self.features['gRNA_U'].cutsite
				self.multiUP=False
			except KeyError:
				pass
		try:
			self.downstream=self.features['gRNA_D5']
			self.multiDN=True

			self.downstream2=self.features['gRNA_D3']
			self.downstreamcut=self.features['gRNA_D3'].cutsite

		except KeyError:
			try:
				self.downstream=self.features['gRNA_D']
				self.downstreamcut=self.features['gRNA_D'].cutsite

				self.multiDN=False
			except KeyError:
				pass

		try:
			self.upstream=self.features['gRNA_E2U']
			self.multiUP=False
		except KeyError:
			pass
		try:
			self.downstream=self.features['gRNA_E2D']
			self.multiDN=False
		except KeyError:
			pass

		try:
			self.upstream
			self.downstream
		except AttributeError:
			print "SOMETHING IS WRONG HERE"
			print self.features.keys()
			sys.exit()

		self.targetopen=35
		if self.multiUP and self.multiDN:
			self.wtStart=self.upstream.start
			self.wtEnd=self.downstream2.end
		elif not self.multiUP and self.multiDN:
			self.wtStart=self.upstream.start
			self.wtEnd=self.downstream2.end
			self.wtEnd=self.downstream2.end
		elif not self.multiUP and not self.multiDN:
			self.wtStart=self.upstream.start
			self.wtEnd=self.downstream.end

		try:
	 		firstcut=filein.upstream.cutsite
	 		secondcut=filein.downstream2.cutsite
	 	except AttributeError:
	 		secondcut=filein.downstream.cutsite
	 	self.delSize=secondcut-firstcut		


	def emDelSizes(self):
		initial=self.delSize-1
		prRange=np.arange(initial,initial+1200,200)
		stringadd=[]
		for j in range(len(prRange)-1):
			stringadd.append(str(prRange[j]+1)+"-"+str(prRange[j+1]))
		self.EMranges=" ".join(stringadd)


	def CreatePrimer3File(self,pr):


		if self.PrimerType == "WT":
			self.target_length=self.targetopen	##g
		elif self.PrimerType == "EM":
			self.target_length=self.delSize+self.targetopen
			self.emDelSizes()
			pr.PRIMER_PRODUCT_SIZE_RANGE="PRIMER_PRODUCT_SIZE_RANGE="+self.EMranges+"\n"

		
		pr.updateParams(self.Attempt)
		if filein.wtStart > 2000:
			self.initCUT=2000
			self.upstreamCUT=filein.wtStart-2000
			self.targetStart=self.wtStart - self.upstreamCUT

		else:
			self.upstreamCUT=0
			self.initCUT=filein.wtStart
			self.targetStart=self.wtStart

		if len(filein.Sequence) - filein.wtEnd > 3000:
			downstreamCUT=filein.wtEnd+3000
		else:
			downstreamCUT=len(filein.Sequence)
		self.cullSeq=filein.Sequence[(self.upstreamCUT):downstreamCUT]
		secondY=self.targetStart-150+self.target_length+150+1
		output=self.cullSeq[:(self.targetStart-150)]+"[" +self.cullSeq[(self.targetStart-150):]	
		Nput=output[:secondY]+"]" +self.cullSeq[secondY:]	
#		print "\n\n\n"
#		print self.PrimerType
#		print Nput
#		print "\n\n\n"
#		print pr.PRIMER_PRODUCT_SIZE_RANGE
		pr.SEQUENCE_TEMPLATE="SEQUENCE_TEMPLATE="+self.cullSeq+"\n"
		self.P3filename=os.path.join(self.RunDir,symbol+"_"+self.PrimerType+".txt")
		io=open(self.P3filename,'w')
		io.write(pr.SEQUENCE_ID+"\n")
		io.write(pr.SEQUENCE_TEMPLATE)
		io.write("SEQUENCE_TARGET="+str(self.targetStart-150)+","+str(self.target_length+150)+"\n")
		io.write(pr.GENERIC_PARAMS_A)
		io.write(pr.MASKING)
		io.write(pr.GENERIC_PARAMS_Aa)
		io.write(pr.PRIMER_MIN_SIZE)
		io.write(pr.GENERIC_PARAMS_B)
		io.write(pr.PRIMER_MIN_TM)
		io.write(pr.GENERIC_PARAMS_C)	
		io.write(pr.PRIMER_MIN_GC)
		io.write(pr.PRIMER_OPT_GC_PERCENT)
		io.write(pr.PRIMER_MAX_GC)
		io.write(pr.PRIMER_PRODUCT_SIZE_RANGE)
		io.write(pr.GENERIC_PARAMS_D)
		io.write(pr.ENDFILE)
		io.close()


	def Primer3exe(self):
	#	print self.filename
		cmd='/usr/bin/primer3_core < '+self.P3filename
		p = Popen(cmd, shell=True,stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
		self.Primer3Output=stdout
		self.Primer3Error=stderr
		outer=open(os.path.join(self.RunDir,self.P3filename+".P3OUT"),'w')
		outer.write(stdout)
		outer.close()

	def FilterPrimers(self):
		# # sequence                       start ln  N   GC%     Tm any_th end_th   pin   sim   lity

		
		primerID=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_[A-Z]+\=")
		primerPlain=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}\=[0-9]{1,6}")
		primerseq=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_SEQUENCE\=")
		primerTM=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_TM\=")
		primerANYTH=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_SELF_ANY_TH\=")
		primerENDTH=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_SELF_END_TH\=")
		primerHAIRPINTH=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_HAIRPIN_TH\=")
		primerGC=re.compile("^PRIMER_(LEFT|RIGHT)_[0-9]{1,3}_GC_PERCENT\=")
		primerPAIR_ANYTH=re.compile("^PRIMER_PAIR_[0-9]{1,3}_COMPL_ANY_TH\=")
		primerPAIR_ENDTH=re.compile("^PRIMER_PAIR_[0-9]{1,3}_COMPL_END_TH\=")


	
		self.PrimerPairs={}
		for line in self.Primer3Output.split("\n"):
			if primerID.match(line):
				try:
					ID,lclClass
					if ID not in self.PrimerPairs:
						lclClass.CheckComplete()
						if lclClass.Complete:
							self.PrimerPairs[ID]=lclClass
					else:
						lclClass=PrimerPairing()
						ID=int(line.split("_")[2].split("=")[0])
				except NameError:
					lclClass=PrimerPairing()
					ID=int(line.split("_")[2].split("=")[0])
			if primerPlain.match(line):
				position,length=line.split("=")[1].split(",")
				if re.search("LEFT",line):
					lclClass.For_START=int(position)
					lclClass.For_LENGTH=int(length)
				else:
					lclClass.Rev_START=int(position)
					lclClass.Rev_LENGTH=int(length)
			elif primerseq.match(line):
				if re.search("LEFT",line):
					lclClass.For_SEQUENCE=line.split("=")[1].strip().upper()
				else:
					lclClass.Rev_SEQUENCE=line.split("=")[1].strip().upper()
			elif primerGC.match(line):
				if re.search("LEFT",line):
					lclClass.For_GC=float(line.split("=")[1].strip())
				else:
					lclClass.Rev_GC=float(line.split("=")[1].strip())
			elif primerANYTH.match(line):
				if re.search("LEFT",line):
					lclClass.For_ANYTH=float(line.split("=")[1].strip())
				else:
					lclClass.Rev_ANYTH=float(line.split("=")[1].strip())
			elif primerENDTH.match(line):
				if re.search("LEFT",line):
					lclClass.For_ENDTH=float(line.split("=")[1].strip())
				else:
					lclClass.Rev_ENDTH=float(line.split("=")[1].strip())
			elif primerTM.match(line):
				if re.search("LEFT",line):
					lclClass.For_TM=float(line.split("=")[1].strip())
				else:
					lclClass.Rev_TM=float(line.split("=")[1].strip())
			elif primerHAIRPINTH.match(line):
				if re.search("LEFT",line):
					lclClass.For_HAIRPINTH=float(line.split("=")[1].strip())
				else:
					lclClass.Rev_HAIRPINTH=float(line.split("=")[1].strip())
			elif primerPAIR_ANYTH.match(line):
				lclClass.PairANYTH=float(line.split("=")[1].strip())
			elif primerPAIR_ENDTH.match(line):
				lclClass.PairENDTH=float(line.split("=")[1].strip())
		##Collect the last
		try:
			lclClass
		except UnboundLocalError:
			print "No primers found"
			self.ValidPairs=False
			self.PrimerPairs[self.PrimerType]={}
			return

		lclClass.CheckComplete()
		if lclClass.Complete:
			self.PrimerPairs[ID]=lclClass
		print "Found total of ",len(self.PrimerPairs)

		WTdistanceBuffer=150
		EMdistanceBuffer=150
		self.ValidPairs=False
		localpairs=[]
		if len(self.PrimerPairs.keys()) == 1:
			ranges=self.PrimerPairs.keys()
		else:
			ranges=range(0,max(map(lambda i: int(i),self.PrimerPairs.keys())))
		for pair in ranges:
			l=self.PrimerPairs[pair]
			ForwardPosition=l.For_START
			ReversePosition=l.Rev_START
#			print self.PrimerType
#			print "FOR",ForwardPosition,self.upstream.cutsite-self.upstreamCUT, ForwardPosition < self.upstream.cutsite-self.upstreamCUT
#			print "REV",ReversePosition,self.downstream.cutsite -self.upstreamCUT, ReversePosition < self.downstream.cutsite -self.upstreamCUT
#			print l.For_SEQUENCE
#			print l.Rev_SEQUENCE
			if self.PrimerType == "WT":
				if self.multiUP and self.multiDN:
					##most common case
					#if (ForwardPosition < self.features['gRNA_U5'].cutsite - self.upstreamCUT-WTdistanceBuffer) and \
					if (ForwardPosition < self.features['gRNA_U5'].cutsite - self.upstreamCUT) and \
					( (ReversePosition > self.features['gRNA_U5'].end - self.upstreamCUT) and \
					(ReversePosition < self.features['gRNA_D5'].cutsite - self.upstreamCUT)):
						localpairs.append(pair)#,self.Forward[pair],self.Reverse[pair]])
		
				elif not self.multiUP and not self.multiDN:
					##second most common
					#if ForwardPosition < self.upstream.cutsite-self.upstreamCUT-WTdistanceBuffer and \
					if ForwardPosition < self.upstream.cutsite-self.upstreamCUT and \
						ReversePosition < self.downstream.cutsite -self.upstreamCUT:
						localpairs.append(pair)#,self.Forward[pair],self.Reverse[pair]])
				elif self.multiUP and not self.multiDN:
					#if ForwardPosition < self.features['gRNA_U5'].cutsite-self.upstreamCUT-WTdistanceBuffer and \
					if ForwardPosition < self.features['gRNA_U5'].cutsite-self.upstreamCUT and \
					( (ReversePosition > self.features['gRNA_U5'].end-self.upstreamCUT) and \
					  (ReversePosition < self.downstream.cutsite-self.upstreamCUT) ):
						localpairs.append(pair)#,self.Forward[pair],self.Reverse[pair]])
			elif self.PrimerType == "EM":
				if self.multiUP and self.multiDN:
					##most common case
					if (ForwardPosition < (self.features['gRNA_U5'].cutsite - self.upstreamCUT)) and \
					(ReversePosition > self.features['gRNA_D3'].cutsite - self.upstreamCUT): 
						localpairs.append(pair)#,self.Forward[pair],self.Reverse[pair]])
				elif not self.multiUP and not self.multiDN:
					##second most common
					if ForwardPosition < self.upstream.cutsite-self.upstreamCUT and \
					(ReversePosition > self.downstream.cutsite -self.upstreamCUT):
						localpairs.append(pair)#,self.Forward[pair],self.Reverse[pair]])
				elif self.multiUP and not self.multiDN:
					if ForwardPosition < self.features['gRNA_U5'].cutsite-self.upstreamCUT and \
					(ReversePosition > self.features['gRNA_D'].cutsite-self.upstreamCUT):
						localpairs.append(pair)#,self.Forward[pair],self.Reverse[pair]])
						
		if not len(localpairs):
			self.ValidPairs=False
		#	print "NO %s PAIRS FOUND" % (self.PrimerType)
			self.PrimerPairs[self.PrimerType]={}
		else:
			self.ValidPairs=True
			print "Found %s %s valid paired primers" % (len(localpairs),self.PrimerType)
			self.PrimerTypedPairs[self.PrimerType]={pair:self.PrimerPairs[pair] for pair in localpairs}
			for rank in localpairs:
				##Go from sequences to Rank, Rank to Sequences
				self.VerifiedSequencePairs[self.PrimerType+"_"+self.PrimerPairs[rank].For_SEQUENCE+"_"+self.PrimerPairs[rank].Rev_SEQUENCE]=self.PrimerPairs[rank]
				self.RankVerifiedPairs[self.PrimerType+"_"+self.PrimerPairs[rank].For_SEQUENCE+"_"+self.PrimerPairs[rank].Rev_SEQUENCE]=rank

	def Primerline(self):
		"""wrap for including a miscellaneous feature  ("AKA primer")"""
		WT_F1="wt_f1"
		if self.PrimerType == "WT":
			idFor="wt_F1"
			idRev="wt_R1"
		else:
			idFor="em_F1"
			idRev="em_R1"
		#print self.feature
		#print self.features
		if idFor not in self.features:

			#####FORWARD
			if self.FOR[self.Forward][-2] == False:
				self.fullfile.insert(self.originPoint,"""     misc_feature    %s..%s\n""" % (self.FOR[self.Forward][-1]+1,self.FOR[self.Forward][-1]+len(self.Forward)))
			else:
				self.fullfile.insert(self.originPoint,"""     misc_feature    complement(%s..%s)\n""" % (self.FOR[self.Forward]+1,self.FOR[self.Forward][-1]+len(self.Forward)))
			self.originPoint+=1
			self.fullfile.insert(self.originPoint,"""                     /locus_tag=\"%s\"\n""" % (idFor,))
			self.originPoint+=1	
			self.fullfile.insert(self.originPoint,"""                     /ApEinfo_label=\"%s\"\n""" % (idFor,))
			self.originPoint+=1	
			self.fullfile.insert(self.originPoint,"""                     /ApEinfo_fwdcolor=\"%s\"\n""" % ("#aa81ff",))
			self.originPoint+=1	
			self.fullfile.insert(self.originPoint,"""                     /ApEinfo_revcolor=\"%s\"\n""" % ("#aa81ff",))
			self.originPoint+=1	

		if idRev not in self.features:
		####REVERSE
			if self.REV[self.Reverse][-2] == False:
				self.fullfile.insert(self.originPoint,"""     misc_feature    %s..%s\n""" % (self.REV[self.Reverse][-1]+1,self.REV[self.Reverse][-1]+len(self.Reverse)))
			else:
				self.fullfile.insert(self.originPoint,"""     misc_feature    complement(%s..%s)\n""" % (self.REV[self.Reverse][-1]+1,self.REV[self.Reverse][-1]+len(self.Reverse)))
			self.originPoint+=1
			self.fullfile.insert(self.originPoint,"""                     /locus_tag=\"%s\"\n""" % (idRev,))
			self.originPoint+=1	
			self.fullfile.insert(self.originPoint,"""                     /ApEinfo_label=\"%s\"\n""" % (idRev,))
			self.originPoint+=1	
			self.fullfile.insert(self.originPoint,"""                     /ApEinfo_fwdcolor=\"%s\"\n""" % ("#aa81ff",))
			self.originPoint+=1	
			self.fullfile.insert(self.originPoint,"""                     /ApEinfo_revcolor=\"%s\"\n""" % ("#aa81ff",))
			self.originPoint+=1	

	def InsertLine(self,newline):
		"""We are wrapping insert with an automated incrememt"""
		self.fullfile.insert(self.insertionPoint,newline)
		self.insertionPoint+=1
		self.originPoint+=1

	def LocatePrimer(self,primerseq):
		cleansedSeq=self.Sequence.strip().upper()
		complement=False
		primerseq=primerseq.upper()
		try:
			position=cleansedSeq.index(primerseq)
		except ValueError:
			revseq=self.RevSeq(primerseq)
			try:
				position=cleansedSeq.index(revseq)
				complement=True
			except ValueError:
				print "NO PRIMER FOUND"
				print primerseq
				print revseq
				sys.exit()
		return [complement,position]


	def PickPrimers(self):

		WT_ranks=self.PrimerTypedPairs['WT'].keys()
		WT_ranks.sort()
		for wtrank in WT_ranks:
			try:
				self.VerifiedSequencePairs["WT"+"_"+self.PrimerTypedPairs['WT'][wtrank].For_SEQUENCE+"_"+self.PrimerTypedPairs['WT'][wtrank].Rev_SEQUENCE]
				WT_pair=self.PrimerTypedPairs['WT'][wtrank]
				break
			except KeyError:
				print "HELP ME"	
				print wtrank,"WT"+"_"+self.PrimerTypedPairs['WT'][wtrank].For_SEQUENCE+"_"+self.PrimerTypedPairs['WT'][wtrank].Rev_SEQUENCE
		
		EM_ranks=self.PrimerTypedPairs['EM'].keys()
		EM_ranks.sort()
		for emrank in EM_ranks:
			try:
				self.VerifiedSequencePairs["EM"+"_"+self.PrimerTypedPairs['EM'][emrank].For_SEQUENCE+"_"+self.PrimerTypedPairs['EM'][emrank].Rev_SEQUENCE]
				EM_pair=self.PrimerTypedPairs['EM'][emrank]
				break
			except KeyError:
				print "HELP ME"	
				print "EM"+"_"+self.PrimerTypedPairs['EM'][emrank].For_SEQUENCE+"_"+self.PrimerTypedPairs['EM'][emrank].Rev_SEQUENCE

		self.WT_pair=WT_pair

		strand,location=self.LocatePrimer(WT_pair.For_SEQUENCE)
		self.WT_pair.For_STRAND=strand
		self.WT_pair.For_LOCATION=location

		strand,location=self.LocatePrimer(WT_pair.Rev_SEQUENCE)
		self.WT_pair.Rev_STRAND=strand
		self.WT_pair.Rev_LOCATION=location

		self.EM_pair=EM_pair
		strand,location=self.LocatePrimer(EM_pair.For_SEQUENCE)
		self.EM_pair.For_STRAND=strand
		self.EM_pair.For_LOCATION=location

		strand,location=self.LocatePrimer(EM_pair.Rev_SEQUENCE)
		self.EM_pair.Rev_STRAND=strand
		self.EM_pair.Rev_LOCATION=location

		self.CommonForward=False
		if (WT_pair.For_SEQUENCE == EM_pair.For_SEQUENCE) and (WT_pair.For_START == EM_pair.For_START):
			self.CommonForward=True
			print "COMMON FORWARD PRIMER"
		#	print WT_pair.For_SEQUENCE, WT_pair.For_START,WT_pair.For_TM
		#	print EM_pair.For_SEQUENCE, EM_pair.For_START,EM_pair.For_TM



	def GrabOutput(self):
		
		lclREGEX=re.compile(self.symbol+"_"+self.PrimerType+"_chr[0-9XY]{1,2}.output")
		##Acad10_EM_chr10.output
		allPr=filter(lambda p: lclREGEX.match(p.split("/")[-1]),glob.glob(os.path.join(self.pcrdir,"*.output")))
		pcr_products={}
		pcr_info={}
		deletelist=[]


		for f in allPr:
			if os.path.getsize(f) == 0:
				os.remove(f)
			else:
				io=open(f,'r').readlines()
				for line in io:
					if line.startswith(">"):
						infoline=line.strip().split()
						rank=float(infoline[1].split(":")[1].split("_")[0])+.00
						isPCRsize=int(infoline[2].strip().rstrip("bp"))
						myPCRsize=int(infoline[1].split(":")[2])
						##This is awkward. But if we have more than 100 off-target products..
						##  lets face it, we won't care if a few get overwritten 
						if rank in pcr_info:
							deletelist.append(rank)		
						pcr_info[rank]=line.strip().split()+[isPCRsize,myPCRsize,isPCRsize==myPCRsize]
						pcr_products[rank]=""
					else:
						pcr_products[rank]+=line.strip()
		#print "Total pcr_products",len(pcr_products)
		deletelist=list(set(deletelist))
		for item in deletelist:
			del pcr_info[item]
			del pcr_products[item]

		#print pcr_info
		all_products=pcr_info.keys()
		all_products.sort()
		#self.unique_pairs=[]
		local_pairs=[]
		print "isPCR finds %s PCR products" % (len(all_products))
		for p in all_products:
			##the second part pcr_info[p][-2] -delSize < 1000 is checking for EM, but will always be true for WT anyways
			if pcr_info[p][-1] and (pcr_info[p][-2]-self.delSize) < 1200:
				##genmonic info,APE_INFO,genomic_DEL+bp,FORWARD,REVERSE,genomic_DEL,APE_DEL,genomic_DEL==APE_DEL?
				forwardSeq=pcr_info[p][3].upper()
				reverseSeq=pcr_info[p][4].upper()
				try:
					self.VerifiedSequencePairs[self.PrimerType+"_"+forwardSeq+"_"+reverseSeq]
					local_pairs.append(self.PrimerType+"_"+forwardSeq+"_"+reverseSeq)
				except KeyError:
					print "FAILED %s pair FOR: %s\t REV:%s " % (self.PrimerType,forwardSeq,reverseSeq)
				#n,sequence,start,length,n,gc,tm,any_th,three_th,hairpin,sim,quality=self.PrimerData[forwardSeq]
				#['27', 'gtcctggaactgaatttagcccttc', '2474', '25', '0', '48.00', '61.607', '0.00', '2.08', '0.00', '12.00', '1.607']
				#local_pairs.append([forwardSeq,reverseSeq])
		if not len(local_pairs):
			self.ValidPairs=False
			return

	def UniqueGenomic(self):
		"""here we run run blat and then call out to the function GrabOutput
		to ensure that we have a unique PCR product. If not we discard and restart process"""
		#print "Attempting to run isPCR - only works on Centos"
		if self.ValidPairs:
			pcrdir="/home/clarkg/PrimerDesign/PCRdir"
			PcrExec='/home/clarkg/PrimerDesign/isPcr'
			prfilename=self.symbol+"_"+self.PrimerType+".primers"
			prfile=open(os.path.join(pcrdir,prfilename),'w')
			self.localPairs=self.PrimerTypedPairs[self.PrimerType]
			for n,pairs in enumerate(self.localPairs):
				#print n,pairs
				cPair=self.PrimerPairs[pairs]
		
				rank=pairs
				ampSize=int(cPair.Rev_START)-int(cPair.For_START)+1
				ampSeq=self.cullSeq[int(cPair.For_START):int(cPair.Rev_START)]
				
				prfile.write("Rank:"+str(rank+1)+"_"+
					self.PrimerType+"_Del:"+
					str(ampSize)+
					":"+
					str(cPair.For_START)+
					"-"+str(cPair.Rev_START)
					+"\t"+cPair.For_SEQUENCE+"\t"+cPair.Rev_SEQUENCE+"\n")
			prfile.close()
			for chrF in self.ChromFiles:
				###OUTPUT FILE = <gene symbol>_[WT|EM]_chr[1-19XY].primers 
				lclout=self.symbol+"_"+self.PrimerType+"_"+chrF.split("/")[-1].split(".")[0]+".output"
				p = Popen([PcrExec,chrF,os.path.join(pcrdir,prfilename),os.path.join(pcrdir,lclout)],stdout=PIPE, stderr=PIPE)
				stdout, stderr = p.communicate()

			self.GrabOutput()


	def CommitPairFile(self):
		self.fullfile=open(self.APEgDNA,'r').readlines()
		lastcomment=re.compile("COMMENT[\W]+[A-z0-9]")
		seqstart=re.compile("^ORIGIN[\W]+$")
#		print self.RankVerifiedPairs
#		return
		for x in range(len(self.fullfile)):
			if lastcomment.match(self.fullfile[x]):
				self.insertionPoint=x
			elif seqstart.match(self.fullfile[x]):
				self.originPoint=x

		### Get proper sense for sequences
		##WT
		if self.WT_pair.For_STRAND:
			self.WT_pair.For_CALC_LOCATION=self.WT_pair.For_LOCATION-self.WT_pair.For_LENGTH
			self.WT_pair.For_SENSE_SEQUENCE=self.RevSeq(self.WT_pair.For_SEQUENCE)
		else:
			self.WT_pair.For_SENSE_SEQUENCE=self.WT_pair.For_SEQUENCE
			self.WT_pair.For_CALC_LOCATION=self.WT_pair.For_LOCATION
			
		
		if self.WT_pair.Rev_STRAND:
			self.WT_pair.Rev_CALC_LOCATION=self.WT_pair.Rev_LOCATION-self.WT_pair.Rev_LENGTH
			self.WT_pair.Rev_SENSE_SEQUENCE=self.RevSeq(self.WT_pair.Rev_SEQUENCE)
		else:
			#self.WT_pair.Rev_SENSE_SEQUENCE=self.WT_pair.Rev_SEQUENCE
			self.WT_pair.Rev_CALC_LOCATION=self.WT_pair.Rev_LOCATION
		###EM
		if self.EM_pair.For_STRAND:
			self.EM_pair.For_CALC_LOCATION=self.EM_pair.For_LOCATION-self.EM_pair.For_LENGTH
			self.EM_pair.For_SENSE_SEQUENCE=self.RevSeq(self.EM_pair.For_SEQUENCE)
		else:
			self.EM_pair.For_SENSE_SEQUENCE=self.EM_pair.For_SEQUENCE
			self.EM_pair.For_CALC_LOCATION=self.EM_pair.For_LOCATION

		if self.EM_pair.Rev_STRAND:
			self.EM_pair.Rev_CALC_LOCATION=self.EM_pair.Rev_LOCATION-self.EM_pair.Rev_LENGTH
			self.EM_pair.Rev_SENSE_SEQUENCE=self.RevSeq(self.EM_pair.Rev_SEQUENCE)
		else:
			self.EM_pair.Rev_CALC_LOCATION=self.EM_pair.Rev_LOCATION
			#self.EM_pair.Rev_SENSE_SEQUENCE=self.EM_pair.Rev_SEQUENCE

		### Finished with this Sense stuff

		print "**WT Primers**"
		print "We had to make %s changes to WT primer specification" % (self.AttemptType['WT'])
		print "Oligo\tSeq\t5'-3'\tGC\ttm"
		print "FOR",self.WT_pair.For_SEQUENCE,self.WT_pair.For_SENSE_SEQUENCE,self.WT_pair.For_TM,self.WT_pair.For_GC
		print "REV",self.WT_pair.Rev_SEQUENCE,self.WT_pair.Rev_SENSE_SEQUENCE,self.WT_pair.Rev_TM,self.WT_pair.Rev_GC


		print "**EM Primers**"
		print "We had to make %s changes to EM primer specification" % (self.AttemptType['EM'])
		print "Oligo\tSeq\t5'-3'\tGC\ttm"
		print "FOR",self.EM_pair.For_SEQUENCE,self.EM_pair.For_SENSE_SEQUENCE,self.EM_pair.For_TM,self.EM_pair.For_GC
		print "REV",self.EM_pair.Rev_SEQUENCE,self.EM_pair.Rev_SENSE_SEQUENCE,self.EM_pair.Rev_TM,self.EM_pair.Rev_GC
		#self.WT_pair_PRODUCT=(self.WT_pair.Rev_CALC_LOCATION+self.WT_pair.Rev_LENGTH)-self.WT_pair.For_CALC_LOCATION
		self.WT_pair_PRODUCT=(self.WT_pair.Rev_LOCATION+self.WT_pair.Rev_LENGTH)-self.WT_pair.For_LOCATION

		#self.WT_EM_pair_PRODUCT=(self.EM_pair.Rev_CALC_LOCATION+self.EM_pair.Rev_LENGTH)-self.EM_pair.For_CALC_LOCATION
		self.WT_EM_pair_PRODUCT=(self.EM_pair.Rev_LOCATION+self.EM_pair.Rev_LENGTH)-self.EM_pair.For_LOCATION

		#self.EM_pair_PRODUCT=(self.EM_pair.Rev_CALC_LOCATION+self.EM_pair.Rev_LENGTH)-self.EM_pair.For_CALC_LOCATION-self.delSize
		self.EM_pair_PRODUCT=(self.EM_pair.Rev_LOCATION+self.EM_pair.Rev_LENGTH)-self.EM_pair.For_LOCATION-self.delSize
		print "*****PRODUCTS*****"
		print "WT Product", self.WT_pair_PRODUCT
		print "EM-WT Procduct", self.WT_EM_pair_PRODUCT
		print "EM Product",self.EM_pair_PRODUCT
		print "Deletion Size",self.delSize
		print "\n\n" 
		#sys.exit()
		#for PT in ["WT","EM"]:
		for PT in ["WT","EM"]:
			self.InsertLine("""COMMENT     \n""")			
			if PT == "WT":
				if self.AttemptType["WT"] > 1:
					self.InsertLine("""COMMENT     ***Deviated from %s Primer Conditions, %s changes made\n\n""" % ("WT",self.AttemptType["WT"]-1,))
				self.InsertLine("""COMMENT     WT Product Size (wt_F1 - wtR1) = %s \n\n""" % (self.WT_pair_PRODUCT,))
				self.InsertLine("""COMMENT     Wild-Type PCR:\n""")
				self.InsertLine("""COMMENT     OLIGO     start     len     tm     gc%     any_th     3'th     hairpin     seq(5' to 3')\n""")

				#self.InsertLine("""COMMENT     LEFT       %s      %s     %s     %s     %s      %s     %s    %s\n""" 
				self.InsertLine("""COMMENT     LEFT       %s      %s     %s     %s     %s      %s     %s    %s\n""" 
					% (str(self.WT_pair.For_LOCATION),self.WT_pair.For_LENGTH,self.WT_pair.For_TM,self.WT_pair.For_GC,
					self.WT_pair.For_ANYTH,self.WT_pair.For_ENDTH,self.WT_pair.For_HAIRPINTH,self.WT_pair.For_SEQUENCE) )

				self.InsertLine("""COMMENT     RIGHT      %s      %s     %s     %s     %s      %s     %s    %s\n""" 
					% (str(self.WT_pair.Rev_LOCATION),self.WT_pair.Rev_LENGTH,self.WT_pair.Rev_TM,self.WT_pair.Rev_GC,
					self.WT_pair.Rev_ANYTH,self.WT_pair.Rev_ENDTH,self.WT_pair.Rev_HAIRPINTH,self.WT_pair.Rev_SEQUENCE))

				self.InsertLine("""COMMENT     \n""")	
				##
				#####FORWARD
				if self.WT_pair.For_STRAND == False:
					self.fullfile.insert(self.originPoint,"""     misc_feature    %s..%s\n""" % 
						(self.WT_pair.For_LOCATION+1,self.WT_pair.For_LOCATION+self.WT_pair.For_LENGTH))
				else:
					self.fullfile.insert(self.originPoint,"""     misc_feature    complement(%s..%s)\n""" % 
						(self.WT_pair.For_LOCATION+1,self.WT_pair.For_LOCATION+self.WT_pair.For_LENGTH))
				self.originPoint+=1
				self.fullfile.insert(self.originPoint,"""                     /locus_tag=\"%s\"\n""" % ("wt_F1",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_label=\"%s\"\n""" % ("wt_F1",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_fwdcolor=\"%s\"\n""" % ("#aa81ff",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_revcolor=\"%s\"\n""" % ("#aa81ff",))
				self.originPoint+=1	

				#print "::::",self.wtR1[1],wtR1_sense,self.wtR1[-2]
				####REVERSE
				if self.WT_pair.Rev_STRAND == False:
					self.fullfile.insert(self.originPoint,"""     misc_feature    %s..%s\n""" %
						(self.WT_pair.Rev_LOCATION+1,self.WT_pair.Rev_LOCATION+self.WT_pair.Rev_LENGTH ))
				else:
					self.fullfile.insert(self.originPoint,"""     misc_feature    complement(%s..%s)\n""" %
						(self.WT_pair.Rev_LOCATION+1,self.WT_pair.Rev_LOCATION+self.WT_pair.Rev_LENGTH))
				self.originPoint+=1
				self.fullfile.insert(self.originPoint,"""                     /locus_tag=\"%s\"\n""" % ("wt_R1",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_label=\"%s\"\n""" % ("wt_R1",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_fwdcolor=\"%s\"\n""" % ("#aa81ff",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_revcolor=\"%s\"\n""" % ("#aa81ff",))
				self.originPoint+=1	
				self.InsertLine("""COMMENT     \n""")	

			else:
				if self.AttemptType["EM"] > 1:
					self.InsertLine("""COMMENT     ***Deviated from %s Primer Conditions, %s changes made\n\n""" % ("EM",self.AttemptType["EM"]-1,))
				if self.CommonForward:
					self.InsertLine("""COMMENT     WT-EM Product Size (wt_F1 - em_R1) = %s \n""" % (self.WT_EM_pair_PRODUCT))
					self.InsertLine("""COMMENT     EM Product Size (wt_F1 - em_R1) = %s \n\n""" % (self.EM_pair_PRODUCT))
					self.InsertLine("""COMMENT     Deletion PCR:\n""")
					self.InsertLine("""COMMENT     \t\t----WT Forward primer is used----\n""")
				self.InsertLine("""COMMENT     OLIGO     start     len     tm     gc%     any_th     3'th     hairpin     seq(5' to 3')\n""")
				if self.CommonForward == False:
					self.InsertLine("""COMMENT     WT-EM Product Size (em_F1 - em_R1) = %s \n""" % (self.WT_EM_pair_PRODUCT))
					self.InsertLine("""COMMENT     EM Product Size (em_F1 - em_R1) = %s \n\n""" % (self.EM_pair_PRODUCT))
					self.InsertLine("""COMMENT     Deletion PCR:\n""")

					self.InsertLine("""COMMENT     LEFT       %s      %s     %s     %s     %s      %s     %s    %s\n"""  %
						 (str(self.EM_pair.For_LOCATION),self.EM_pair.For_LENGTH,self.EM_pair.For_TM,self.EM_pair.For_GC,
						self.EM_pair.For_ANYTH,self.EM_pair.For_ENDTH,self.EM_pair.For_HAIRPINTH,self.EM_pair.For_SEQUENCE) )

				self.InsertLine("""COMMENT     RIGHT      %s      %s     %s     %s     %s      %s     %s    %s\n""" 
					% (str(self.EM_pair.Rev_LOCATION),self.EM_pair.Rev_LENGTH,self.EM_pair.Rev_TM,self.EM_pair.Rev_GC,
					self.EM_pair.Rev_ANYTH,self.EM_pair.Rev_ENDTH,self.EM_pair.Rev_HAIRPINTH,self.EM_pair.Rev_SEQUENCE))
				self.InsertLine("""COMMENT     \n""")	

				if self.CommonForward == False:
					if self.EM_pair.For_STRAND == False:
						self.fullfile.insert(self.originPoint,"""     misc_feature    %s..%s\n""" % 
							(self.EM_pair.For_LOCATION+1,self.EM_pair.For_LOCATION+self.EM_pair.For_LENGTH ))
					else:
						self.fullfile.insert(self.originPoint,"""     misc_feature    complement(%s..%s)\n""" % 
							(self.EM_pair.For_LOCATION+1,self.EM_pair.For_LOCATION+self.EM_pair.For_LENGTH ))

					self.originPoint+=1
					self.fullfile.insert(self.originPoint,"""                     /locus_tag=\"%s\"\n""" % ("em_F1",))
					self.originPoint+=1	
					self.fullfile.insert(self.originPoint,"""                     /ApEinfo_label=\"%s\"\n""" % ("em_F1",))
					self.originPoint+=1	
					self.fullfile.insert(self.originPoint,"""                     /ApEinfo_fwdcolor=\"%s\"\n""" % ("#aa81ff",))
					self.originPoint+=1	
					self.fullfile.insert(self.originPoint,"""                     /ApEinfo_revcolor=\"%s\"\n""" % ("#aa81ff",))
					self.originPoint+=1	

				####REVERSE
				if self.EM_pair.Rev_STRAND == False:
					self.fullfile.insert(self.originPoint,"""     misc_feature    %s..%s\n""" % 
						(self.EM_pair.Rev_LOCATION+1,self.EM_pair.Rev_LOCATION+self.EM_pair.Rev_LENGTH ))
				else:
					self.fullfile.insert(self.originPoint,"""     misc_feature    complement(%s..%s)\n""" % 
						(self.EM_pair.Rev_LOCATION+1,self.EM_pair.Rev_LOCATION+self.EM_pair.Rev_LENGTH ))
				self.originPoint+=1
				self.fullfile.insert(self.originPoint,"""                     /locus_tag=\"%s\"\n""" % ("em_R1",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_label=\"%s\"\n""" % ("em_R1",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_fwdcolor=\"%s\"\n""" % ("#aa81ff",))
				self.originPoint+=1	
				self.fullfile.insert(self.originPoint,"""                     /ApEinfo_revcolor=\"%s\"\n""" % ("#aa81ff",))
				self.originPoint+=1	
				self.InsertLine("""COMMENT     \n""")	
	

		gRNA_only=filter(lambda p: re.search("gRNA",p),self.features)
		self.InsertLine("""COMMENT     \n""")
		
		self.InsertLine("""COMMENT     All Potential Deletion Sizes:\n""")
		self.CutSizes={}
		for x1 in range(len(gRNA_only)):
			g1=gRNA_only[x1]
			for x2 in range(x1+1,len(gRNA_only)):
				g2=gRNA_only[x2]
				lclDel=abs(self.features[g1].cutsite-self.features[g2].cutsite)
				#self.InsertLine("""COMMENT     \tDeletion : %s <-> %s = %s\n""" % (g1, g2,abs(self.features[g1].cutsite-self.features[g2].cutsite)))
				if self.CommonForward:
					self.InsertLine("""COMMENT     \tProduct Size  %s <-> %s = %s\n\n""" % (g1,g2,abs(self.EM_pair.For_START-self.EM_pair.Rev_START)+1-lclDel),)
				else:
					self.InsertLine("""COMMENT     \tProduct Size  %s <-> %s = %s\n\n""" % (g1,g2,abs(self.EM_pair.For_START-self.EM_pair.Rev_START)+1-lclDel),)
				
				#self.InsertLine("""COMMENT     \t\tDeletion : %s <-> %s = %s\n""" % (g1, g2,abs(self.features[g1].cutsite-self.features[g2].cutsite)))
				#self.InsertLine()	
				#print g,filein.features[g].cutsite
					
		self.InsertLine("""COMMENT     \n""")
		username=getpass.getuser()
		current=str(datetime.date.today())
		if username == 'root':	
			if not os.path.isdir(os.path.join(self.fileDIR,current+"_"+"AutomatedDesign")):
				os.makedirs(os.path.join(self.fileDIR,current+"_"+"AutomatedDesign"))
			#else:
			
			local=open(os.path.join(self.fileDIR,current+"_"+"AutomatedDesign",self.symbol+"_gDNA_masked.ape" ),"w")
		else:
			local=open(os.path.join("/home/clarkg/PrimerDesign",self.symbol+"_gDNA_masked.ape" ),"w")
		local.write("".join(self.fullfile))
		local.close()

class APEfeatures(APEfile):
	def __init__(self,name):
		self.name=name #this will be the locus first...then we change to gRNA_U5m,etc

	def grabFeatures(self,remaininglines):
		start,end=self.name.split("..")
		direction="+"
		if self.name.endswith("_rev"):
			end=end.rstrip("_rev")
			direction="-"
		#print start,end,direction
		self.start=int(start)
		self.end=int(end)
		#
		self.direction=direction
		for x,line in enumerate(remaininglines):
			if self.misc.match(line) or self.exontag.match(line) or line.startswith("ORIGIN"):
				#just save time/space
				break
			if self.locus.match(line):
				self.feature_desc=line.strip().split("=")[1].strip("\"")
			elif self.label.match(line):
				self.feature_desc=line.strip().split("=")[1].strip("\"")
			if self.misc.match(line) or self.exonlabel.match(line) or line.startswith("ORIGIN"):
				#just save time/space
				break
			elif self.exonlabel.match(line):
				#print "LINE HERE",line
				self.feature_desc=line.strip().split("=")[1].strip("\"")
	
		self.FeatureSequence=filein.Sequence[(self.start-1):self.end]
		if direction == "+" and re.search("gRNA",self.feature_desc):
			#print self.feature_desc,direction,filein.Sequence[(self.start-1):self.end+3]
			self.cutsite=self.end-3
		elif direction =="-" and re.search("gRNA",self.feature_desc):
			#print self.feature_desc,direction,filein.Sequence[(self.start-4):self.end]
			self.cutsite=self.start+2


if __name__ == "__main__":
	import time
	import argparse
	#parser = argparse.ArgumentParser(description="Add symbols as arguments to file")
	#parser.add_argument('--genes', type=list,nargs="?",help='Add gene symbols to design as a comma-separated list')
	#for symbol in ["Sept11","Zcwpw1"]:
	#for symbol in ["C2cd3","Dnah2"]:
	#for symbol in ["Hist1h2bc","Peli2","Tnk2","Kdm3b","Klhl12",

	#for symbol in ["A2ml1","Slco6d1","Polr3d","Serpine1","Slc35a1","Pcdhb13"]:
	for symbol in ["Irx2","F8","Pdxdc1"]:#"Gm15155"]:

		print symbol
		filein=APEfile(symbol)
		filein.findFile()

		filein.readAPESequence()
		filein.readAPEFeatures()
		filein.ValidPairs=False
		filein.PrimerSelect()

		if filein.fileFound:
			print "SYMBOL:",symbol
			filein.symbol=symbol
			for primertype in ["WT","EM"]:
				filein.PrimerType=primertype
				filein.ValidPairs=False
				print "Parsing %s primers" % (primertype,)
				## Just instantiating
				pr=Primer(symbol)
				filein.Attempt=0		
				while filein.Attempt < 9 and not filein.ValidPairs:
					if filein.Attempt >= 1:
						print "Changed %s Primer Parameters, stage %s" % (filein.PrimerType,filein.Attempt)
					#### Create input file for primer3_core
					##   In CreatePrimer3File, we will sometimes shorten the Sequence.
					##   This must be reflected in the co-ordinates of the APE file features.
					##   They should all be shorted by the distance to the upstream gRNA start point.
					##   i.e. gRNA_U5, gRNA_U, or gRNA_E2U(?)
					filein.CreatePrimer3File(pr)

					#### Primer3exe will make the stdout/stderr wrapped subprocess call out to primer 3
					##   The results will be in the "RunData" directory 
					filein.Primer3exe()

					#### FilterPrimers will grab the output from primer3 and ensure that the primer
					##   pairs are suitable to our co-ordinates
					filein.FilterPrimers()

					#### UniqueGenomic makes a stdout/stderr subprocess call to isPCR
					##   If a primer pair is found not to be unique it is removed from the 
					##   list of available primer pairs and we continue (if that list is populated)
					filein.UniqueGenomic()

					filein.Attempt+=1
					filein.AttemptType[primertype]+=1
					#### We have either found the primers and can proceed
					##   If we do not have any valid pairs, we increase our Attempt number
					##   Increasing the attempt number will change the primer3 input file parameters
	
		#### At this point we have a list of EM and WT primers that provide unique products
		if not filein.ValidPairs:	
			print "FAILED TO FIND PRIMERS"
		else:
			filein.PickPrimers()
			filein.CommitPairFile()
