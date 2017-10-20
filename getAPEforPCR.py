#!/usr/bin/env python

import string, re, sys
import glob,os
import subprocess
from primerDesign import Primer
from subprocess import Popen, PIPE
import operator

#import APEread

class APEfile:
	if not glob.glob("/mnt/cmmr/"):
		print "Must mount /mnt/cmmr/ using 'rootmount.sh' in home directory. Ensure that passwords for user have no changed"
		sys.exit()
	misc=re.compile("[\W]+misc_feature[\W]+")
	locus=re.compile("^[\W]+/locus_tag=")
	label=re.compile("^[\W]+/label=")
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
				self.fileFound=True

	def readAPESequence(self):
		self.OriginStart=False
		fullfile=open(self.APEgDNA,'r').readlines()
		self.Sequence=""
		for x,line in enumerate(fullfile):
			if line.startswith("ORIGIN"):
				self.OriginStart=True
			if self.OriginStart:
				seqline=line.split()[1:]
				if len(seqline):
					self.Sequence+="".join(seqline)

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
			if len(filein.Sequence) - filein.wtEnd > 2000:
				downstreamCUT=filein.wtEnd+2000
			else:
				downstreamCUT=len(filein.Sequence)
			newSeq=filein.Sequence[(self.upstreamCUT):downstreamCUT]
			pr.SEQUENCE_TEMPLATE+=newSeq+"\n"
			self.filename=os.path.join(os.getcwd(),symbol+"_em.txt")

		elif self.PrimerType == "EM":
			try:
				firstcut=filein.upstream.cutsite
				secondcut=filein.downstream2.cutsite
			except AttributeError:
				secondcut=filein.downstream.cutsite

			delSize=secondcut-firstcut
			print delSize
			self.target_length=500 +delSize	##to ensure we are 200 bp away from either cut site
			pr=Primer(symbol+"_em")
			pr.
			pr.updateParams(self.Attempt)
			#delSize=wtEnd-wtStart
			if filein.wtStart > 2000:
				self.initCUT=1800
				self.upstreamCUT=filein.wtStart-2000

			else:
				self.upstreamCUT=0
				self.initCUT=filein.wtStart-200
			if len(filein.Sequence) - filein.wtEnd > 2000:
				downstreamCUT=filein.wtEnd+2000
			else:
				downstreamCUT=len(filein.Sequence)	
			if len(str(delSize+500)) == 4:
				minSize=int(round(delSize+500,-3))
			elif len(str(delSize+500)) == 3:
				minSize=int(round(delSize+500),-2)
			print minSize
			maxSize=minSize+1000
			newranges=
			
			pr.self.PRIMER_PRODUCT_SIZE_RANGE="PRIMER_PRODUCT_SIZE_RANGE= "
			newSeq=filein.Sequence[(self.upstreamCUT):downstreamCUT]
			pr.SEQUENCE_TEMPLATE+=newSeq+"\n"
			self.filename=os.path.join(os.getcwd(),symbol+"_em.txt")
		filein.WriteFile(pr)

	def Primer3exe(self):
		print self.filename
		cmd='/usr/bin/primer3_core < '+self.filename
		p = Popen(cmd, shell=True,stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
		print stderr

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
			FOR[frank+"_"+fsequence]=fordata

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
				REV[rrank+"_"+rsequence]=revdata
				thermoScore=round(sum([rGC+fGC+rTm+fTm]),1)
				structScore=round(sum([rend_th+rany_th+fend_th+rend_th+rhairpin+fhairpin]),1)

				combinedquality=round(rquality+fquality,1)
				distance=rstart-fstart
				if distance > 100 and distance < 1000  and combinedquality < 7 and thermoScore < 10 and structScore < 5:
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
		#for not self.ValidPairs:


	def FormatEMPrimers(self):
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
			FOR[frank+"_"+fsequence]=fordata

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
				REV[rrank+"_"+rsequence]=revdata
				thermoScore=round(sum([rGC+fGC+rTm+fTm]),1)
				structScore=round(sum([rend_th+rany_th+fend_th+rend_th+rhairpin+fhairpin]),1)

				combinedquality=round(rquality+fquality,1)
				distance=rstart-fstart
				if distance > 100 and distance < 1000  and combinedquality < 7 and thermoScore < 10 and structScore < 5:
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
	
		self.direction=direction
		for x,line in enumerate(remaininglines):
			if self.misc.match(line) or self.exontag.match(line) or line.startswith("ORIGIN"):
				#just save time/space
				break
			if self.locus.match(line):
				self.feature_desc=line.strip().split("=")[1].strip("\"")
			elif self.label.match(line):
				self.feature_desc=line.strip().split("=")[1].strip("\"")
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
	#for symbol in ["Atp2c2","Plin4","Hist1h2bc","Peli2","Sspo","Tnk2","Kdm3b","Klhl12","Mmadhc"]:
	#for symbol in ["Peli2","Sspo","Tnk2","Kdm3b","Klhl12","Mmadhc"]:
	for symbol in ["Acad10","Aass","Cct5","Cox20","Ifit3","Mcm6"]:
		filein=APEfile(symbol)
		filein.findFile()
		if filein.fileFound:
			print "SYMBOL:",symbol+":"
			print filein.APEgDNA
			filein.readAPESequence()
			filein.readAPEFeatures()

			filein.PrimerType="EM"
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
				#print filein.PrimerPair
				print "FORWARD",filein.Forward,filein.ForwardStart,"\t",filein.upstream.cutsite
				print "REVERSE",filein.Reverse,filein.ReverseStart,"\t",filein.downstream.cutsite
				try:
					firstcut=filein.upstream.cutsite
					secondcut=filein.downstream2.cutsite
				except AttributeError:
					secondcut=filein.downstream.cutsite

				delSize=secondcut-firstcut
				delSeq=filein.Sequence[:firstcut]+filein.Sequence[secondcut:]
			
				# for k,j in filein.features.iteritems():
				# 	if filein.features[k].start >= firstcut:
				# 		print "Redefine feature",delSize,k,filein.features[k].start-delSize,filein.features[k].end-delSize
				# 	if re.search("gRNA",k) and (filein.features[k].start > filein.ForwardStart and filein.features[k].end < filein.ReverseStart):
				# 		print "Amplifying",k,filein.features[k].start,filein.features[k].end

			else:
				print "FAILED TO FIND PRIMERS After much effort"
			#break
			time.sleep(4)
			print "\n"
			#break