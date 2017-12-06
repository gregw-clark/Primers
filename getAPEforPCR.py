#!/usr/bin/env python

import string, re, sys
import glob,os
import subprocess
from primerDesign import Primer
from subprocess import Popen, PIPE
import operator



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
		self.ValidPairs=False

	def WriteFile(self,pr):
		##pr is our primer instance. it should be changing with repeated attempts
		io=open(self.filename,'w')
		io.write(pr.SEQUENCE_ID+"\n")
		io.write(pr.SEQUENCE_TEMPLATE)
		io.write("SEQUENCE_TARGET="+str(self.initCUT)+","+str(self.target_length)+"\n")
		io.write(pr.GENERIC_PARAMS)
		io.write(pr.PRIMER_MIN_GC)
		io.write(pr.PRIMER_MAX_GC)
		io.write(pr.PRIMER_MIN_SIZE)
		io.write(pr.HAIRPIN)
		io.write(pr.SELF_TH)
		io.write(pr.END_TH)
		io.write(pr.PRIMER_MIN_TM)
		io.write(pr.PRIMER_MAX_TM)
		io.write(pr.PRIMER_PRODUCT_SIZE_RANGE)
		io.write(pr.MASKING)
		io.write(pr.ENDFILE)
		io.close()

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
		pairs={'A':'T' , 'C':'G', 'G':'C' , 'T':'A', 'N':'N'}
		lclS=sequence[::-1]
		rev=[]
		for s in lclS:
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
					print "COORDINATES MESSED"
					print line
					sys.exit()
			# elif line.startswith("ORIGIN"):
			# 	self.OriginStart=True
			# if self.OriginStart:
			# 	seqline=line.split()[1:]
			# 	if len(seqline):
			# 		self.Sequence+="".join(seqline)


	def PrimerSelect(self):

		##We want to capture the largest potential cutsite
		##	APE files have a variety of ways of defining gRNA sites (via the user)
		##	So we are necessarily handcuffed here at defining these sites
		try:
			self.upstream=self.features['gRNA_U5']
			self.multiUP=True
			self.upstream2=self.features['gRNA_U3']
		except KeyError:
			try:
				self.upstream=self.features['gRNA_U']
				self.multiUP=False
			except KeyError:
				pass
		try:
			self.downstream=self.features['gRNA_D5']
			self.multiDN=True
			self.downstream2=self.features['gRNA_D3']
		except KeyError:
			try:
				self.downstream=self.features['gRNA_D']
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

		if self.multiUP:
			self.wtStart=self.upstream.start
			self.wtEnd=self.upstream2.end
			#self.exclude=downstream.end
		elif not self.multiUP and self.multiDN:
			self.wtStart=self.upstream.end
			self.wtEnd=self.downstream.start
			#self.exclude=downstream2.start
		elif not self.multiUP and not self.multiDN:
			self.wtStart=self.upstream.start
			self.wtEnd=self.upstream.end


	def CreatePrimer3File(self):
	
		if self.PrimerType == "WT":
			self.target_length=filein.wtEnd-filein.wtStart
			pr=Primer(symbol)
			pr.updateParams(self.Attempt)
			if filein.wtStart > 2000:
				self.initCUT=2000
				self.upstreamCUT=filein.wtStart-2000

			else:
				self.upstreamCUT=0
				self.initCUT=filein.wtStart
				self.initCUT=filein.wtStart
			if len(filein.Sequence) - filein.wtEnd > 2000:
				downstreamCUT=filein.wtEnd+2000
			else:
				downstreamCUT=len(filein.Sequence)

			newSeq=filein.Sequence[(self.upstreamCUT):downstreamCUT]
			pr.SEQUENCE_TEMPLATE+=newSeq+"\n"
			self.filename=os.path.join(os.getcwd(),symbol+".txt")
			io=open(self.filename,'w')
			io.write(pr.SEQUENCE_ID+"\n")
			io.write(pr.SEQUENCE_TEMPLATE)
			io.write("SEQUENCE_TARGET="+str(self.initCUT)+","+str(self.target_length)+"\n")
			io.write(pr.GENERIC_PARAMS)
			io.write(pr.HAIRPIN+"\n")
			io.write(pr.SELF_TH+"\n")
			io.write(pr.END_TH+"\n")
			io.write(pr.MASKING+"\n")
			io.write(pr.ENDFILE)
			io.close()
	
	def Primer3exe(self):
		cmd='/usr/bin/primer3_core < '+self.filename
		p = Popen(cmd, shell=True,stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()

	def FormatPrimers(self):
		# # sequence                       start ln  N   GC%     Tm any_th end_th   pin   sim   lity
		FOR={}
		REV={}
		pairs=[]

		if self.multiUP and self.multiDN:
			##most common case
		 	self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.features['gRNA_U5'].start - self.upstreamCUT),self.ForwardFile)
		 	self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.features['gRNA_U3'].end - self.upstreamCUT) and int(e.strip().split()[2]) < (self.features['gRNA_D5'].end - self.upstreamCUT),self.ReverseFile) 
		elif not self.multiUP and not self.multiDN:
		 	##second most common
		 	self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.upstream.start-self.upstreamCUT),self.ForwardFile)
		 	self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) < (self.downstream.start -self.upstreamCUT),self.ReverseFile)
		elif self.multiUP and not self.multiDN:
			self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.features['gRNA_U5'].start-self.upstreamCUT),self.ForwardFile)
			self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.features['gRNA_U3'].start-self.upstreamCUT) and int(e.strip().split()[2]) < (self.downstream.start-self.upstreamCUT),self.ReverseFile)
		elif self.multiDN and not multiUP:
			print "DO SOMETHING ABOUT THIS"
		
			newSeq=filein.Sequence[(self.upstreamCUT):downstreamCUT]
			pr.SEQUENCE_TEMPLATE+=newSeq+"\n"
			self.filename=os.path.join(os.getcwd(),symbol+"_em.txt")
			self.delSize = 0


		elif self.PrimerType == "EM":
			try:
				firstcut=filein.upstream.cutsite
				secondcut=filein.downstream2.cutsite
			except AttributeError:
				secondcut=filein.downstream.cutsite

			delSize=secondcut-firstcut
			#print "Deletion size is = %s" % (delSize,)
			self.target_length=delSize	##to ensure we are 200 bp away from either cut site
			pr=Primer(symbol+"_em")

			pr.updateParams(self.Attempt)

			#delSize=wtEnd-wtStart
			if filein.wtStart > 2000:
				self.initCUT=1800
				self.upstreamCUT=filein.wtStart-2000

			else:
				self.upstreamCUT=0
				self.initCUT=filein.wtStart-200

			if len(filein.Sequence) - filein.wtEnd > 3000:
				self.downstreamCUT=filein.wtEnd+3000
				#print filein.wtEnd,"WE RIGHT HERE"
			else:
				self.downstreamCUT=len(filein.Sequence)	

		#	print self.initCUT
			self.delSize=delSize
			maxSize=delSize+2000
		#	sys.exit()
			newrange=str(100)+"-"+str(maxSize)+"\n"
			#print maxSize
			pr.PRIMER_PRODUCT_SIZE_RANGE="PRIMER_PRODUCT_SIZE_RANGE= "+newrange
			newSeq=filein.Sequence[(self.upstreamCUT):self.downstreamCUT]
			#print newSeq
			#sys.exit()
			pr.SEQUENCE_TEMPLATE+=newSeq+"\n"
			self.filename=os.path.join(os.getcwd(),symbol+"_em.txt")
		filein.WriteFile(pr)

	def Primer3exe(self):
		#print self.filename
		cmd='/usr/bin/primer3_core < '+self.filename
		p = Popen(cmd, shell=True,stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
		#print stderr

	def FormatPrimers(self):
		# # sequence                       start ln  N   GC%     Tm any_th end_th   pin   sim   lity
		self.FOR={}
		self.REV={}
		pairs=[]

		if self.PrimerType == "WT":

			if self.multiUP and self.multiDN:
				##most common case
			 	self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.features['gRNA_U5'].start - self.upstreamCUT),self.ForwardFile)
			 	self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.features['gRNA_U3'].end - self.upstreamCUT) and int(e.strip().split()[2]) < (self.features['gRNA_D5'].end - self.upstreamCUT),self.ReverseFile) 
			elif not self.multiUP and not self.multiDN:
			 	##second most common
			 	self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.upstream.start-self.upstreamCUT),self.ForwardFile)
			 	self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) < (self.downstream.start -self.upstreamCUT),self.ReverseFile)
			elif self.multiUP and not self.multiDN:
				self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.features['gRNA_U5'].start-self.upstreamCUT),self.ForwardFile)
				self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.features['gRNA_U3'].start-self.upstreamCUT) and int(e.strip().split()[2]) < (self.downstream.start-self.upstreamCUT),self.ReverseFile)
			elif self.multiDN and not multiUP:
				print "DO SOMETHING ABOUT THIS"

		elif self.PrimerType == "EM":
			if self.multiUP and self.multiDN:
				##most common case
			 	self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.features['gRNA_U5'].start - self.upstreamCUT-200),self.ForwardFile)
			 	self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.features['gRNA_D5'].end - self.upstreamCUT+200),self.ReverseFile) 
			elif not self.multiUP and not self.multiDN:
			 	##second most common
			 	self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.upstream.start-self.upstreamCUT-200),self.ForwardFile)
			 	self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.downstream.start -self.upstreamCUT+200),self.ReverseFile)
			elif self.multiUP and not self.multiDN:
				self.ForwardFile=filter(lambda s: int(s.strip().split()[2]) < (self.features['gRNA_U5'].start-self.upstreamCUT-200),self.ForwardFile)
				self.ReverseFile=filter(lambda e: int(e.strip().split()[2]) > (self.features['gRNA_U3'].start-self.upstreamCUT+200) and int(e.strip().split()[2]) < (self.downstream.start-self.upstreamCUT),self.ReverseFile)
			elif self.multiDN and not multiUP:
				pass

		for k in range(len(self.ForwardFile)):
			fordata=self.ForwardFile[k].strip().split()
			frank,fsequence,fstart,flength,fN,fGC,fTm,fany_th,fend_th,fhairpin,flibsim,fquality=fordata
			fstart=int(fstart)
			fGC=abs(50-float(fGC))
			fTm=abs(60-float(fTm))
			flength=abs(25-int(flength))
			fany_th=float(fany_th)
			fend_th=float(fend_th)
			fhairpin=float(fhairpin)
			fquality=float(fquality)
			fsequence=fsequence.strip().upper()
			self.FOR[fsequence]=fordata

			for l in range(len(self.ReverseFile)):
				revdata=self.ReverseFile[l].strip().split()

				rrank,rsequence,rstart,rlength,rN,rGC,rTm,rany_th,rend_th,rhairpin,rlibsim,rquality=revdata
				rstart=int(rstart)
				#if rstart < self.downstream:
				rGC=abs(50-float(rGC))
				rTm=abs(60-float(rTm))
				rlength=abs(25-int(rlength))
				rany_th=float(rany_th)
				rend_th=float(rend_th)
				rhairpin=float(rhairpin)
				rquality=float(rquality)
				rsequence=rsequence.strip().upper()
				self.REV[rsequence]=revdata
				thermoScore=round(sum([rGC+fGC+rTm+fTm]),1)
				structScore=round(sum([rend_th+rany_th+fend_th+rend_th+rhairpin+fhairpin]),1)

				combinedquality=round(rquality+fquality,1)
				distance=rstart-fstart
				if distance > 100 and distance < 1000  and combinedquality < 7 and thermoScore < 10 and structScore < 5:
					#print combinedquality,distance, thermoScore,structScore,"\t",fstart,fsequence,"\t",rstart,rsequence
					pairs.append([combinedquality,distance,thermoScore,structScore,fstart,fsequence,rstart,rsequence])
				if self.PrimerType =="WT":
					if distance > 200 and distance < 1000  and combinedquality < 10 and thermoScore < 10 and structScore < 10:
						#print combinedquality,distance, thermoScore,structScore,"\t",fstart,fsequence,"\t",rstart,rsequence
						pairs.append([combinedquality,distance,thermoScore,structScore,fstart,fsequence,rstart,rsequence])
				elif self.PrimerType=="EM":
					if distance > 200 and distance < (1000+ self.delSize)  and combinedquality < 10 and thermoScore < 10 and structScore < 10:
					#print combinedquality,distance, thermoScore,structScore,"\t",fstart,fsequence,"\t",rstart,rsequence
						pairs.append([combinedquality,distance,thermoScore,structScore,fstart,fsequence,rstart,rsequence])

		pairs=sorted(pairs,key=operator.itemgetter(0,2,3,1))
		if not len(pairs):
			self.ValidPairs=False
		else:
			self.PrimerPair=pairs[0]
			if pairs[0][4] > pairs[0][6]:
				self.Forward=pairs[0][7]
				self.ForwardStart=pairs[0][6]+ self.upstreamCUT
				self.Reverse=pairs[0][5]
				self.ReverseStart=pairs[0][4]+ self.upstreamCUT
			else:
				self.Forward=pairs[0][5]
				self.ForwardStart=pairs[0][4]+ self.upstreamCUT
				self.Reverse=pairs[0][7]
				self.ReverseStart=pairs[0][6] + self.upstreamCUT
			self.ValidPairs=True

	def GetPrimerLists(self):
		if os.path.exists(self.symbol+".for"):
			self.ForwardFile=open(self.symbol+".for",'r').readlines()
			self.ForwardFile=self.ForwardFile[3:]
		else:
			self.ForwardFile=False
		if os.path.exists(self.symbol+".rev"):
			self.ReverseFile=open(self.symbol+".rev",'r').readlines()



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

	def GetPrimerLists(self):
		if self.PrimerType == "WT":
			forfile=self.symbol+".for"
			revfile=self.symbol+".rev"
		elif self.PrimerType == "EM":
			forfile=self.symbol+"_em.for"
			revfile=self.symbol+"_em.rev"

		if os.path.exists(forfile):
			self.ForwardFile=open(forfile,'r').readlines()
			self.ForwardFile=self.ForwardFile[3:]
		else:
			self.ForwardFile=False
		if os.path.exists(revfile):
			self.ReverseFile=open(revfile,'r').readlines()
			self.ReverseFile=self.ReverseFile[3:]
		else:
			self.ReverseFile=False
	
		if self.ForwardFile and self.ReverseFile:
			# # sequence                       start ln  N   GC%     Tm any_th end_th   pin   sim   lity
			self.FormatPrimers()
		else:
			self.ValidPairs=False

	def InsertLine(self,newline):
		"""We are wrapping insert with an automated incrememt"""
		self.fullfile.insert(self.insertionPoint,newline)
		self.insertionPoint+=1
		self.originPoint+=1

	def LocatePrimer(self,primerseq):
		cleansedSeq=self.Sequence.strip().upper()
		complement=False
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


	def CommitFile(self):
		#self.fullfile=open(self.APEgDNA,'r').readlines()
		lastcomment=re.compile("COMMENT[\W]+[A-z0-9]")
		seqstart=re.compile("^ORIGIN[\W]+$")
		for x in range(len(self.fullfile)):
			if lastcomment.match(self.fullfile[x]):
				self.insertionPoint=x
			elif seqstart.match(self.fullfile[x]):
				self.originPoint=x
		try:
			forwardstats=self.FOR[self.Forward]
		except KeyError:
			try:
				forwardstats=self.REV[self.Forward]
			except KeyError:
				forwardstats=None
				pass
		try:
			reversestats=self.REV[self.Reverse]
		except KeyError:
			try:
				reversestats=self.FOR[self.Reverse]
			except KeyError:
				reversestats=None

		#print self.Sequence
		#upseq=
		##Lets add the location that WE find for the features

		self.FOR[self.Forward]+=self.LocatePrimer(self.Forward)
		self.REV[self.Reverse]+=self.LocatePrimer(self.Reverse)


		self.InsertLine("""COMMENT     \n""")			
		#self.InsertLine("""COMMENT     %s PRODUCT SIZE: %s\n""" % (self.PrimerType, self.product))
		if self.PrimerType == "WT":
			self.InsertLine("""COMMENT     Wild-Type PCR:\n""")
		else:
			self.InsertLine("""COMMENT     Deletion PCR:\n""")

		self.InsertLine("""COMMENT     OLIGO     start     len     tm     gc%     any_th     3'th     hairpin     seq\n""")
		self.InsertLine("""COMMENT     LEFT       %s      %s     %s     %s     %s      %s     %s    %s\n""" 
			% (str(self.ForwardStart),forwardstats[3],forwardstats[6],forwardstats[5],forwardstats[7],forwardstats[8],forwardstats[9],forwardstats[1]) )

		self.InsertLine("""COMMENT     RIGHT      %s      %s     %s     %s     %s      %s     %s    %s\n""" 
			% (str(self.ReverseStart),reversestats[3],reversestats[6],reversestats[5],reversestats[7],reversestats[8],reversestats[9],reversestats[1])       )
		self.InsertLine("""COMMENT     \n""")			

		gRNA_only=filter(lambda p: re.search("gRNA",p),self.features)
		#print gRNA_only
		if self.PrimerType != "WT":
			self.InsertLine("""COMMENT     %s PRODUCT SIZE: %s\n""" % ("WT_EM", self.WTEMproduct))
			self.InsertLine("""COMMENT     %s PRODUCT SIZE: %s\n""" % (self.PrimerType, self.EMproduct))
			self.InsertLine("""COMMENT     \n""")

			self.InsertLine("""COMMENT     All Potential Deletion Sizes:\n""")
			for x1 in range(len(gRNA_only)):
				g1=gRNA_only[x1]
				for x2 in range(x1+1,len(gRNA_only)):
					g2=gRNA_only[x2]
					self.InsertLine("""COMMENT     %s <-> %s : %s\n""" % (g1, g2,abs(self.features[g1].cutsite-self.features[g2].cutsite)))
					#self.InsertLine()	
					#print g,filein.features[g].cutsite
		else:
			self.InsertLine("""COMMENT     %s PRODUCT SIZE: %s\n""" % ("WT", self.WTproduct))
						
		self.InsertLine("""COMMENT     \n""")

		self.Primerline()
		#fullfile.insert(insertionPoint,"""COMMENT     %s Reverse \t %s\n""" %(primertype,filein.Reverse))
	#	print reversestats
	#	if self.PrimerType != "WT":
	#		local=open("/home/clarkg/IntentShare/TargetedTEST.ape",'w')
	#	else:
		local=open(os.path.join("/home/clarkg/IntentShare/",self.symbol+"_AltNull.ape" ),"w")
		local.write("".join(self.fullfile))
		local.close()
		#print len(fullfile)
	#	sys.exit()


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
				print "LINE HERE",line
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
	for symbol in ["Tcf12"]:
		filein=APEfile(symbol)
		filein.findFile()

		if filein.fileFound:
			print "SYMBOL:",symbol+":"
			filein.symbol=symbol
			for primertype in ["WT","EM"]:
				filein.readAPESequence()
				filein.readAPEFeatures()
				#print filein.APEgDNA,primertype
				filein.ValidPairs=False
				filein.PrimerType=primertype
				filein.PrimerSelect()
				#print filein.upstream.start
				#print filein.downstream.start
				####In CreatePrimer3File, we will sometimes shorten the Sequence.
				##		This must be reflected in the co-ordinates of the APE file features.
				##		They should all be shorted by the distance to the upstream gRNA start point.
				## 		i.e. gRNA_U5, gRNA_U, or gRNA_E2U(?)
				while filein.Attempt < 5 and not filein.ValidPairs:
					if filein.Attempt > 1:
						print "Changed Primer Parameters, stage %s" % (filein.Attempt)
					filein.CreatePrimer3File()
					filein.Primer3exe()
					filein.GetPrimerLists()
					filein.Attempt+=1


				if filein.ValidPairs:
					if filein.Attempt > 1:
						print "Primer parameters for %s changed, %s attempts passed.\n" % (primertype,filein.Attempt)
					#print filein.PrimerPair
					print primertype+" Forward",filein.Forward,filein.ForwardStart,"\t",filein.upstream.cutsite
					print primertype+" Reverse",filein.Reverse,filein.ReverseStart,"\t",filein.downstream.cutsite
					try:
						firstcut=filein.upstream.cutsite
						secondcut=filein.downstream2.cutsite
					except AttributeError:
						secondcut=filein.downstream.cutsite

					filein.delSize=secondcut-firstcut
					if primertype == "WT":
						filein.WTproduct=abs(filein.ReverseStart- filein.ForwardStart) +1
					if primertype == "EM":
						filein.WTEMproduct=abs(filein.ReverseStart - filein.ForwardStart) +1
						filein.EMproduct=abs(filein.ReverseStart - filein.ForwardStart) - filein.delSize +1
					#print "%s PRODUCT SIZE = %s" % (primertype,product)
					filein.delSeq=filein.Sequence[:firstcut]+filein.Sequence[secondcut:]
					filein.CommitFile()
					#filein.Valid
					#print filein.PrimerType,delSize
	
					# gRNA_only=filter(lambda p: re.search("gRNA",p),filein.features)
					# print gRNA_only
					# for g1 in range(len(gRNA_only)):
					# 	for g2 in range(g1,len(gRNA_only)):

					# 	print g,filein.features[g].cutsite
					#for k,j in filein.features.iteritems():
					#	pass
					#	if filein.features[k].start >= firstcut:
					 	#print "Redefine feature",delSize,k,filein.features[k].start-delSize,filein.features[k].end-delSize
					# 	if re.search("gRNA",k) and (filein.features[k].start > filein.ForwardStart and filein.features[k].end < filein.ReverseStart):
					# 		print "Amplifying",k,filein.features[k].start,filein.features[k].end

				else:
					print "FAILED TO FIND PRIMERS After much effort"
				time.sleep(4)
				print "\n"
		break
