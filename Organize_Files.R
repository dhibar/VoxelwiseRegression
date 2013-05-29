###########################################################################
#
# This script performs association tests of a single SNP, a set of SNPs individually, or a set of SNPs as a group at each point within a user-provided 
#     full brain mask. This script is appropriate for testing SNP effects on Jacobian maps outputted from tensor-based morphometry (TBM),
#     grey matter density maps from voxel-based morphometry (VBM), fractional anisotropy (FA) maps (or related measures) from diffusion imaging (DTI),
#     contrast maps from fMRI tasks or any phenotype that you can iterate through given a mask.
#
# These scripts require that you have the following software installed and available on your system: plink, R, R package dcemriS4
#
# These scripts were tested and working on Mac OS X and Ubuntu Linux, modifications may be required to run on a Windows machine
#
# # # #  In order to run this script, you will need to edit the options below and generate the appropriate files as detailed in the Manual.txt
#
#  1) table file #req columns: subjectID
Table="/Users/dhibar/Desktop/Regression/table.txt"
#  2) full path to mask file: mask file has to be in the same file format as the rest of your images 
Maskfile="/Users/dhibar/Desktop/Regression/brain_mask.nii.gz" #supported formats: .nii , .nii.gz , .img/.hdr, and .img.gz/.hdr.gz
#  3) path to plink file: can be in binary format (.bed) or text format (.ped)
PlinkFile="/Users/dhibar/Desktop/Regression/my_subjects_genotypes.ped"
#  4) snps of interest: semicolon seperated list of snps 
SNP="rs6265;rs11030101;rs76327806"
#  5) if you have listed multiple SNPs and would like to perform a combined association analysis with Principal Components Regression choose "no"
Individually="no"  #if you only have one SNP or would like to test each SNP effect individually set this to "yes"
#  6) what model should be used to code the SNPs (usually Additive, but can also be Dominant or Recessive)
Modeltype="Additive"
#  7) full path for the output directory for all of the data no / at the end of the folder name
Outputdirectory="/Users/dhibar/Desktop/Regression/BDNFgene"
#
#
#  Filtering criteria
#  8) column name from the table file with the path to the images
Imagepaths="imagepaths"
#  9) list of column names from the table file to be used as covariates
Covariates="Age;Sex" #semicolon separated list of table headers
#  10) filters for doing more complex subsetting of the table file
Filters=""
#
# Written by Derrek Hibar, Neda Jahanshad, and Jason Stein (Version 0.1; May 21, 2013)
#
##############################################################################

#Create output directory and subfolder to store voxelwise output
dir.create(Outputdirectory);
OutputdirectoryEX = paste(Outputdirectory,'/vLM/',sep='');
dir.create(OutputdirectoryEX);

#Open the table with all the info
Table0<-read.table(Table,header=T,blank.lines.skip = TRUE)
Columnnames = colnames(Table0);
DesignMatrix<-data.frame(subjectID=Table0$subjectID, imgpaths=Table0[,which(Columnnames==Imagepaths)])

# Find number of subjects, number of covariates, number of filters
Nsub=length(as.vector(DesignMatrix$subjectID))
Ncov=length(parse(text=Covariates))
Nfilters=length(parse(text=Filters))

## read and parse covariates if they exist
## add to design Matrix
if (Ncov > 0) {
	Cov<-as.matrix(0,nrow=Nsub,ncol=Ncov)
	ParsedCov=parse(text=Covariates)

	for (nc in 1:Ncov) {
			CovName<-as.character(ParsedCov[nc])
			print(CovName)
			Origcolnames = colnames(DesignMatrix);
			DesignMatrix[,length(DesignMatrix)+1] = Table0[,which(Columnnames==CovName)]
			colnames(DesignMatrix) = c(Origcolnames,CovName);
	}
	print(paste("Done adding covariates",ParsedCov))
}

## read and parse filters if they exist
## add to design Matrix
if (Nfilters > 0) {
	ParsedFilt=parse(text=Filters)

	for (nf in 1:Nfilters) {
		FiltName<-as.character(ParsedFilt[nf])
		print(FiltName)
		Origcolnames = colnames(DesignMatrix);
		DesignMatrix[,length(DesignMatrix)+1] = Table0[,which(Columnnames==FiltName)]
		colnames(DesignMatrix) = c(Origcolnames,FiltName);
	}
	print(paste("Done adding filter",ParsedFilt))
}

## Drop subjects with missing values or NA's in a filter
i=which(apply(DesignMatrix,1,function(x)any(is.na(x))));
print(which(apply(DesignMatrix,1,function(x)any(is.na(x)))))
##get rid of all rows with NAs in them
if (length(i) > 0){
DesignMatrix<-DesignMatrix[-which(apply(DesignMatrix,1,function(x)any(is.na(x)))),]
}

#Store column names from design matrix, will be using this file from now on
Columnnames = colnames(DesignMatrix);

# Store the new total of subjects after removing NA's
Nsub=length(as.vector(DesignMatrix$subjectID)) 

#######
#######

# determine plink file type
PlinkFiletype=substr(PlinkFile,nchar(PlinkFile)-2,nchar(PlinkFile))
PlinkFile=substr(PlinkFile,1,nchar(PlinkFile)-4)

# parse the SNP data
AllSNP=parse(text=SNP)
Snplist=as.matrix(as.character(AllSNP))
write.table(Snplist,file=paste(Outputdirectory,'/',"SNP.list",sep=''),quote=F,row.names=F,col.names=F)

# call plink and generate additive raw files, this step will error if none of your SNPs are found in your plink files
if (PlinkFiletype=="bed") {
	system(paste('plink --noweb  --bfile ',PlinkFile,' --recodeA --extract ',paste(Outputdirectory,'/',"SNP.list",sep=''), ' --out ',paste(Outputdirectory,'/',"SNPs",sep='')));
}
if (PlinkFiletype=="ped") {
	system(paste('plink --noweb --file ',PlinkFile,' --recodeA --extract ',paste(Outputdirectory,'/',"SNP.list",sep=''),' --out ',paste(Outputdirectory,'/',"SNPs",sep='')));
}

cat('Plink output file found here:\n')
print(paste(Outputdirectory,'/',"SNPs",'.raw',sep=''))

## add SNP info to design matrix
SNPinfo<-read.table(paste(Outputdirectory,'/',"SNPs",'.raw',sep=''),header=T,blank.lines.skip = TRUE)
PlinkColumnnames = colnames(SNPinfo);
Origcolnames = colnames(DesignMatrix);

## Get row locations for the subjects in the DesignMatrix that also have genetic data
matchind = match(DesignMatrix$subjectID,SNPinfo$IID); #must have plink file coded with unique IIDs

# loop through the SNPs in the PLINK file and add them to the DesignMatrix
# reformat the SNPs based on the 
for(tt in 7:ncol(SNPinfo)){
	if (Modeltype=='Additive') {
	        print('    Using an additive model\n');
	        DesignMatrix[,length(DesignMatrix)+1] = SNPinfo[,tt][matchind];
	} else if (Modeltype == 'Dominant') {
	        print('    Using a dominant model\n');
			DesignMatrix[,length(DesignMatrix)+1] = SNPinfo[,tt][matchind];
	        zerosnp = length(which(DesignMatrix[,ncol(DesignMatrix)]==0));
	        twosnp = length(which(DesignMatrix[,ncol(DesignMatrix)]==2));
	        if (zerosnp > twosnp) {
				DesignMatrix[which(DesignMatrix[,ncol(DesignMatrix)]==2),] = 1
	        } else {
				DesignMatrix[which(DesignMatrix[,ncol(DesignMatrix)]==0),] = 1
				DesignMatrix[which(DesignMatrix[,ncol(DesignMatrix)]==2),] = 0
			}
	} else if (Modeltype == 'Recessive') {
	        print('    Using a recessive model for minor allele\n');
			DesignMatrix[,length(DesignMatrix)+1] = SNPinfo[,tt][matchind];
	        zerosnp = length(which(DesignMatrix[,ncol(DesignMatrix)]==0));
	        twosnp = length(which(DesignMatrix[,ncol(DesignMatrix)]==2));
	        if (zerosnp > twosnp) {
				DesignMatrix[which(DesignMatrix[,ncol(DesignMatrix)]==0),] = 1
				DesignMatrix[which(DesignMatrix[,ncol(DesignMatrix)]==2),] = 0
	        } else {
				DesignMatrix[which(DesignMatrix[,ncol(DesignMatrix)]==2),] = 1
			}
	}
}


snpname=PlinkColumnnames[7:length(PlinkColumnnames)];
colnames(DesignMatrix) = c(Origcolnames,PlinkColumnnames[7:length(PlinkColumnnames)]);
numsubjects = length(DesignMatrix$subjectID);
cat('    There are',numsubjects,'subjects in the Design Matrix\n')

#######################################################################
#write out some useful info into txt file:
zz <- file(paste(Outputdirectory,'/',"SNP_information_file.txt",sep=""),"w")
writeLines(paste("   Plink file is:",PlinkFile),con=zz,sep="\n")
writeLines(paste("   Modeltype is:",Modeltype),con=zz,sep="\n")
writeLines(paste("   SNP is:",SNP),con=zz,sep="\n")
writeLines(paste('   There are',numsubjects,'subjects with image data'),con=zz,sep="\n")
writeLines(paste("   the image paths are: "),con=zz,sep="\n")
writeLines(paste("              ",DesignMatrix$imgpaths),con=zz,sep="\n")

if (Ncov > 0){
	for (n in 1: Ncov){
		writeLines(paste("   Covariate:",ParsedCov[n]),con=zz,sep="\n")
	}
}
if (Nfilters > 0){
	for (n in 1: Nfilters){
		writeLines(paste("   Filter:",ParsedFilt[n]),con=zz,sep="\n")
	}
}
close(zz)
#########################################################################

# write out a text version of the design matrix for debugging
write.table(DesignMatrix,paste(Outputdirectory,'/DesignMatrix.txt',sep=''),sep="\t",row.names=F,col.names=T);

# generate files needed to run the association tests
#######################################################################
#Create files to run one per slice for each of the SNPs

library(dcemriS4)

#Are there multiple SNPs?
if(length(AllSNP) > 1){
	#Run SNPs as a combined test?
	if(Individually == "no"){
		nSNPs=1
		CombTest=TRUE
	} else {
		nSNPs=length(AllSNP)
		CombTest=FALSE
	}
} else {
	nSNPs=1
	CombTest=FALSE
}

#Find out the file type
ngz=substr(Maskfile,nchar(Maskfile)-2,nchar(Maskfile))
wgz=substr(Maskfile,nchar(Maskfile)-5,nchar(Maskfile)-3)
if(ngz == "nii" || wgz == "nii"){
	readimage = readNIfTI
} else {
	readimage = readANALYZE
}

#Read in the maskfile
Mask = readimage(Maskfile)

#Get the number of slices from the z direction in the mask
numslices = dim(Mask)[3]

#Store the image values
imagefnames = DesignMatrix$imgpaths

for(q in 1:nSNPs){
	if(CombTest == FALSE){
		
		covind=match(ParsedCov,colnames(DesignMatrix))
		snp=DesignMatrix[,(max(covind)+nSNPs)]
		xs=as.matrix(cbind(DesignMatrix[,covind],snp))
		
		#Loop over each slice
		for (slice in 1:numslices) {
			
			#Output to the user
			cat('   Working on slice ',slice,' of ',numslices,' for SNP ', AllSNP[q],'\n',sep="");

			#create an output slice
			outslice = nifti(matrix(data=1,nrow=dim(Mask)[1],ncol=dim(Mask)[2]),datatype=16);
			outsliceT = nifti(matrix(data=0,nrow=dim(Mask)[1],ncol=dim(Mask)[2]),datatype=16);
			outsliceB = nifti(matrix(data=0,nrow=dim(Mask)[1],ncol=dim(Mask)[2]),datatype=16);
	
			#Get the chosen slice from the mask
			Ymask = Mask[,,slice];

			#find all the nonzero voxels
			Ymaskind = which(Ymask>0);

			#If no nonzero voxels, loop again
			if (length(Ymaskind)==0) {
				cat("skip\n");
				next;
			}
	
			#save the current workspace
			save.image(file=paste(OutputdirectoryEX,slice,"_",AllSNP[q],".Rdata",sep=""));

			#Write out the process for each slice to a text file
			write("library(dcemriS4);",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""));
			write(paste("load('",OutputdirectoryEX,slice,"_",AllSNP[q],".Rdata')",sep=""),file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write(paste("slice <- ",slice,sep=""),file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("#Loop over subjects",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("for (subj in 1:length(imagefnames)) {",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    #create a matrix for the Jac files in each slice",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    if (subj==1) {",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("        Yjac = array(0,dim=c(length(imagefnames),length(Ymaskind)));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    }",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    #read in each subject",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write(paste("    deblankfile = sub('^[[:space:]]*(.*?)[[:space:]]*$', '","\\","\\","1', as.character(imagefnames[subj]), perl=TRUE);",sep=''),file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    cat('Working on file:',as.character(imagefnames[subj]),'\n');",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    Ynewjac = readimage(imagefnames[subj]);",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    Ynewjac = Ynewjac[,,slice];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    #Get all the values within the mask",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    Yjac[subj,1:length(Ymaskind)] = Ynewjac[Ymaskind];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("}",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("#Loop over voxels within the mask in the slice to run the regression",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);

			write("for (vox in 1:length(Ymaskind)) {",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    full <- summary(lm(Yjac[,vox] ~ xs))", file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    outslice[Ymaskind[vox]] = full$coefficients[ncol(xs)+1,4];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    outsliceT[Ymaskind[vox]] = full$coefficients[ncol(xs)+1,3];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    outsliceB[Ymaskind[vox]] = full$coefficients[ncol(xs)+1,1];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("}",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write(" ",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);

			write("#Write the output image",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("writeNIfTI(outslice,paste(OutputdirectoryEX,'AssocSlice','_',AllSNP[q],slice,sep=''));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("writeNIfTI(outsliceT,paste(OutputdirectoryEX,'AssocSliceT','_',AllSNP[q],slice,sep=''));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("writeNIfTI(outsliceB,paste(OutputdirectoryEX,'AssocSliceB','_',AllSNP[q],slice,sep=''));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
      
			#Write an output file to run the R script
			write("#!/bin/bash",file=paste(OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""));
			write(" ",file=paste(OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""),append=T);
			write(paste("R --no-save --quiet < ",OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),file=paste(OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""),append=T);
			system(paste("chmod a+x ",OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""));
		}
	} else {
		
		AllSNP="CombSNPs"
		
		covind=match(ParsedCov,colnames(DesignMatrix))
		DesignMatrix=DesignMatrix[complete.cases(DesignMatrix),]
		ped=DesignMatrix[,(max(covind)+1):ncol(DesignMatrix)]
		
		#Store the image values
		imagefnames = DesignMatrix$imgpaths
		
		#######################################################################
		#Calculate all of the variances for each snp in the gene
		z.inter = 0;
		for (ttx in c(1:(ncol(ped)))){
			var.all = (ped[,ttx] - mean(ped[,ttx]))/sd(ped[,ttx]);
			z.inter = cbind(z.inter, var.all);
		}
		z.vars = z.inter[,2:(ncol(z.inter))];

		#Run singular value decomposition on the snp matrix and extract the columns
		#representing 80% of the variance
		z.svd = svd(z.vars);
		denom = sum(z.svd$d);
		numer = 0;
		count = 1;
		var.svd = 0;

		while (var.svd <= 0.8){
			numer = numer + z.svd$d[count];
			var.svd = numer/denom;
			count = count + 1;
		}
		to.include=count-1;
		
		xs = cbind(DesignMatrix[,covind],z.svd$u[,1:to.include])
		xs.red = cbind(DesignMatrix[,covind])
		
		#Loop over each slice
		for (slice in 1:numslices) {

			#Output to the user
			cat('   Working on slice ',slice,' of ',numslices,' for SNP ', AllSNP[q],'\n',sep="");

			#create an output slice
			outslice = nifti(matrix(data=1,nrow=dim(Mask)[1],ncol=dim(Mask)[2]),datatype=16);
			outsliceF = nifti(matrix(data=0,nrow=dim(Mask)[1],ncol=dim(Mask)[2]),datatype=16);
	
			#Get the chosen slice from the mask
			Ymask = Mask[,,slice];

			#find all the nonzero voxels
			Ymaskind = which(Ymask>0);

			#If no nonzero voxels, loop again
			if (length(Ymaskind)==0) {
				cat("skip\n");
				next;
			}
	
			#save the current workspace
			save.image(file=paste(OutputdirectoryEX,slice,"_",AllSNP[q],".Rdata",sep=""));

			#Write out the process for each slice to a text file
			write("library(dcemriS4);",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""));
			write(paste("load('",OutputdirectoryEX,slice,"_",AllSNP[q],".Rdata')",sep=""),file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write(paste("slice <- ",slice,sep=""),file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("#Loop over subjects",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("for (subj in 1:length(imagefnames)) {",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    #create a matrix for the Jac files in each slice",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    if (subj==1) {",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("        Yjac = array(0,dim=c(length(imagefnames),length(Ymaskind)));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    }",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    #read in each subject",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write(paste("    deblankfile = sub('^[[:space:]]*(.*?)[[:space:]]*$', '","\\","\\","1', as.character(imagefnames[subj]), perl=TRUE);",sep=''),file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    cat('Working on file:',as.character(imagefnames[subj]),'\n');",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    Ynewjac = readimage(imagefnames[subj]);",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    Ynewjac = Ynewjac[,,slice];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    #Get all the values within the mask",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    Yjac[subj,1:length(Ymaskind)] = Ynewjac[Ymaskind];",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("}",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("#Loop over voxels within the mask in the slice to run the regression",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("for (vox in 1:length(Ymaskind)) {",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    full <- lm(Yjac[,vox] ~ ., data=as.data.frame(xs))", file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    reduced <- lm(Yjac[,vox] ~ ., data=as.data.frame(xs.red))", file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    results = anova(full, reduced);", file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    outslice[Ymaskind[vox]] = results$\"Pr(>F)\"[2];", file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("    outsliceF[Ymaskind[vox]] = results$F[2];", file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("}",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write(" ",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("#Write the output image",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("writeNIfTI(outslice,paste(OutputdirectoryEX,'AssocSlice','_',AllSNP[q],slice,sep=''));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
			write("writeNIfTI(outsliceF,paste(OutputdirectoryEX,'AssocSliceF','_',AllSNP[q],slice,sep=''));",file=paste(OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),append=T);
      
			#Write an output file to run the R script
			write("#!/bin/bash",file=paste(OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""));
			write(" ",file=paste(OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""),append=T);
			write(paste("R --no-save --quiet < ",OutputdirectoryEX,"AssocSlice",slice,"_",AllSNP[q],".R",sep=""),file=paste(OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""),append=T);
			system(paste("chmod a+x ",OutputdirectoryEX,"RunAssocSlice",slice,"_",AllSNP[q],".sh",sep=""));
		}
	}
}

