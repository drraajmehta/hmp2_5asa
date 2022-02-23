# hmp2_5asa
R/python/unix scripts for 5-ASA analysis 

##### READ ME FOR FILES #####

#Data processing:
#---------------# 
#metadata_create.R		--> makes metadata file from input, NB: LOCF and NOCB used for steroid data 
#hmp2omics_preprocessingbb3.R	--> reads files from ibdmdb website and trims/etc 
###uniref_proptable_dnabb3.R
###uniref_proptable_rnabb3.R
#/PRISM/prismprocessing.R 	--> precursor to Fig 2c

#Overview (results, section 1):
#------------------------------#
#table1_v3.R			--> makes table 1 + supplementary Figure 1
#5ASAfigsv5.pptx		--> Fig 1 

#Metabolomics analysis (results, section 2): 
#------------------------------------------#
#metabPCOA_clean.R 		--> Fig2A [X]
#DA_test_species.R 		--> Supp Fig 2A
###graphlan_prep 		   --> precursor to supp fig 2
#graphlan_1217fileforplot.txt 	--> Supp fig 2B (done on galaxy server)
#ROCmetabs_newclean.R 		--> Fig2B [X], Supp Fig 3a, 3c, 3d, 4, Supp table 2, Supp Fig 5 (
	#reads in my_adj_list85.txt, mergemethod2.txt for network) 
###hallaprep_clean.R		   --> precursor to supp fig 3b
#hallabb3.py 			--> supp fig 3b (done on FAS bc software is installed there) 
#biotransf_clean.R		--> Fig2C, Supp Fig 6 [X] --> reads in metab_biotransf3.txt + PRISMfile.txt 
#/PRISM/prismmetab_reviewed.R 	--> Fig 2C, Supp Fig 7 [X]

#metabs_newuserextra.R --> this is info for text about new steroid and new biologic users 
#chemdraw was used to make the structures 

#Gene discovery (results, section 3): 
#-----------------------------------#
#DA_test_uniref_clean.R		--> Fig 3a, 3c
#sen_spec_clean		--> Fig 3b
#jalview.txt			--> Fig 3d,e
#5ASAfigsv5.ppt		--> Fig 3f

#diamond3.sh			--> negative homology search (done on FAS) 
###diamondoverlap.R		--> provides overlap 	
	#--> output out4.tsv, overlap -> diamond_uniref90.txt 


###### LIST OF UNIREFS
####### uniprot-yourlist_M20210323A94466D2655679D1FD8953E075198DA8185E589.txt

#Biochemistry (results, section 4): 
#-----------------------------------#
#genecontext.R 		--> Supp Fig 10
###proteindomains.R
#chimerascript2.txt		--> Fig 4d

#chemdraw was used to make structures
#jared made panels 4a,b,c supp fig 11 

#Clinical (results, section 5): 
#-----------------------------#
#clinical_clean.R		--> Fig 5a, 5b, supp tables 5,6
### 
#####predicated upon 
########1) /home/rsm34/5asabb3/archive/metadata_create.R #to define steroid users 
########2) preprocessing + metadata_load (unadulterated) 
########3) /home/rsm34/5asabb3/genomics/hmp2genomes.R #for host genetics 
########4) /home/rsm34/5asabb3/clinical_nonuser_clean.R #for prevalence among non users/ORs 
