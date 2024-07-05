#### AMPLISAS OUTPUT ANALYSIS
#### Tomek Gaczorek
#### tomasz.gaczorek@doctoral.uj.edu.pl
#### 17.11.20

#packages - additionally ape, msa and stringr (not loaded generally as they cover some other functions)
library(ggtree)
library(ggplot2)
library(dplyr)
library(seqinr)
library(phangorn)

source("/home/tomek/Dropbox/MHC_projekt/Shared_MHC_project/Amplisas/R_code/Functions_to_be_loaded.R")
# ^ ADJUST!!! - path to the script with created functions de novo
#setwd("~/Dropbox/MHC_projekt/Shared_MHC_project/Podarcis/Amplisas_results/MHC_I/Podarcis_MHCI_171120_first_nm")
# ^ ADJUST!!! - folder with the general result from the amplisas
analysed_file <- "2_second_Natrix_MHCII_080423.xlsx"
# ^ ADJUST!!! - the general result from the amplisas

## BODY #### - GO STEP BY STEP ####

#create log file
log_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_log.txt",sep = "")
writeLines(c(analysed_file,""),log_name)

#reading
data_list <- read_and_adjust(analysed_file) #old sequencing method
#data_list <- read_and_adjust_II(analysed_file,split_var = "-",pattern_var = "20|21|22")
# ^ additional arguments to set colnames to the individual IDs only
# ^ split_var - value that separate the elements of previous name
# ^ pattern_var - the pattern that will unambigously distinguish the proper IDs
# ^ 17 - Anguis; 16 - Salmo; 16 - Solea
data <- data_list[[1]] # actual data
#colnames(data) <- sapply(sapply(colnames(data),strsplit,"_"),"[",1)
#View(data) #check if everything is ok
appendix <- data_list[[2]] # table with information about coverage
#colnames(appendix) <- sapply(sapply(colnames(appendix),strsplit,"_"),"[",1)

#descriptive
write_to_log("DESCRIPTION",body = paste("Number of individuals: ",ncol(data)-6,"\nNumber of sequences: ",
                                        nrow(data),"\nMean number of alleles per individual: ",
                                        round(mean(as.numeric(apply(data[,7:ncol(data)],2,function(x){sum(!is.na(x))}))),2),
                                        "\nMean coverage per individual: ",round(mean(as.numeric(appendix[1,2:ncol(appendix)])),2)),log_name)

#sort by length
data <- arrange(data,LENGTH)
most_frequent(data$LENGTH) #printing most frequent value
#View(data) #check the proper length, usually the most frequent one
desired_length <- most_frequent(data$LENGTH)
#desired_length <- 224

#filtering by length
data <- filter_by_length(data,desired_length,2,rem_shifts = F,log_name) 
# ^ digit refers to the acceptable deviation from the proper length (number od codons)

#filtering by amplicon depth (total coverage of an individual)
quantile(as.numeric(appendix[1,2:ncol(appendix)]),probs = seq(0,1,0.05)) #distribution of total coverage
#View(appendix)
data <- filter_by_amplicon_depth(data,appendix,log_name,minimum_depth = 1000)
# ^ usually more than 300; discarded should not constitute more than 5 %

#limit data to species
#desired_species <- c("cristatus","dobrogicus","macedonicus",'ivanbureschi',"anatolicus")
#data <- limit_to_species(data,"../../Triturus_ids.csv",desired_species,log_name)

#check if there are some reverse complements
#you must have a reference
#data <- check_reversal(data,analysed_file) #without reference
data <- check_reversal_alt(data,analysed_file)

#### here CONTAMINATION (delation can change the PAFs)
#contamination_checking(data,"../../Triturus_ids.csv",dup_symbols = "d")

#additionally removed individuals (OPTIONAL, you must have good reason to delete individual e.g. contamination)
#ids <- seq(16117,16251,1)
#imiona <- colnames(data)[7:ncol(data)]
#keep <- c(ids,paste(ids,"d",sep = ""))
#rem <- imiona[!(imiona %in% keep)]
#rem <- out
# data <- remove_inds(data,c(14791,14627,13157,14624,14663,14784,"14784d"),log_name,
#                     mess = "INDIVIDUALS REMOVED BASED ON TOO HIGH ALLELES NUMBER OR SEQUENCING PROBLEMS")
# data <- remove_inds(data,c(18132),log_name,
#                     mess = "INDIVIDUALS REMOVED BASED ON TOO HIGH ALLELES NUMBER OR SEQUENCING PROBLEMS")

# adding alleles frequencies
list_of_things <- adding_freqs_info(data,appendix)
#list_of_things <- adding_freqs_info_allele(data,appendix) #allele based PAF
data <- list_of_things[[1]]
freqs_data <- list_of_things[[2]]
write.csv(freqs_data,"freqs_data.csv",row.names = F,fileEncoding = "UTF-8")

# plotting max frequencies
#plot_name <-  paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_freqPlot.png",sep = "")
#plot_max_freq(data,plot_name,w = 40,t = .1) # w refers to width of picture, the highest the better resolution,
# ^ t - threshold of plotted points (lower than)
data <- filter_by_tag_switch(data,freqs_data,threshold = .25,log_name) #single threshold
#data <- filter_by_tag_switch_dual(data,freqs_data,threshold = 4.7,dual = 0,log_name)
# ^ filter out the sequences which are poorly cover (probable artefacts)
# ^ dual - additional threshold for alleles within individuals (alleles above PAF)

#modifying names
#data <- modify_names(data,"MHCIex2")
# ^ 3 first letters from original file + provided pattern + number
data <- modify_names_new(data,abr = "Nat",gene = "MHCIIex2",order = 2)

#are there any Ns in sequences?
Are_Ns(data)

#optional - remove some sequences
data <- remove_some_seqs(data,alleles = c("Nat_MHCIIex2_499","Nat_MHCIIex2_240","Nat_MHCIIex2_358","Nat_MHCIIex2_536"),log = log_name,pat_or_vec = "VEC")

#writing sequences to fasta
fasta_name <- save_to_fasta(data,paste("withFrameShifts",analysed_file,sep = "_"))

# OPTIONAL - writing modified file (with_frame shits)
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_modified_FS.csv",sep = "")
write.csv(data,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)
#writing appendix
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_appendix_FS.csv",sep = "")
write.csv(appendix,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)

#optional filtering based on frame shifts
frame_shifts <- detect_frame_shifts(data,desired_length) #checking frame shifts based on length
data <- remove_some_seqs(data,alleles = frame_shifts,log = log_name,pat_or_vec = "VEC")


#repeatability (OPTIONAL, you must have individuals sequenced twice;
#by convention we add 'd' at the end of ID to distinguish duplicated one)
repeatability(data,log_name)

#descriptive
write_to_log("DESCRIPTION",body = paste("Number of individuals: ",ncol(data)-9,"\nNumber of sequences: ",
                                        nrow(data),"\nMean number of alleles per individual: ",
                                        round(mean(as.numeric(apply(data[,10:ncol(data)],2,function(x){sum(!is.na(x))}))),2)),log_name)

#writing modified file (after whole filtering)
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_modified.csv",sep = "")
write.csv(data,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)
#writing appendix
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_appendix.csv",sep = "")
write.csv(appendix,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)

#writing sequences to fasta
fasta_name <- save_to_fasta(data,analysed_file)

#fasta work
#fasta_bin <- seqinr::read.fasta("4_Emys_MHCII_160521_forth_nm.fasta") ### backdoor

allignment_name_file <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned.fasta",sep = "")
#alligned <- allign(fasta_bin,allignment_name_file) #alligned sequences with ClustalW
#alligned <- align_mafft(fasta_name,allignment_name_file,rev_comp = F)
#### MUST BE CHECKED MANUALLY AND PROBLEMS WITH TRANSLATIONS MUST BE RESOLVED

fasta_bin <- seqinr::read.fasta(allignment_name_file,forceDNAtolower = F)
frame_translate <- 2
with_stop_codons <- detect_stop_codons(fasta_bin,log_name,frame_val = frame_translate,stop_verify = F)
#with_stop_codons <- c()
#with_stop_codons <- c(with_stop_codons,"Emy_MHCIex2_068") #add some stops
# ^ names with stop codons - CHANGE FRAME IF NEEDED!!! - you can expect only few such cases
frame_shifts <- frame_shifts[!(frame_shifts %in% with_stop_codons)]
# ^ from this point sequences with stop codons are not treated as frame shifts

#building tree
#alligned <- ape::read.FASTA("1_Lissotriton_MHCI_first_080321_om_alligned.fasta") ### backdoor
tree_nm <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree.new",sep = "")
# not sure how it would behave under Windows as the argument mc_cores within bootstrap.pml function 
# (look into build_tree()) should be used only under UNIX based systems. 
# Without this limit my computer was stacked. In case of it maybe turn off 
# parralelization and just wait longer (multi argument)
alligned <- ape::read.FASTA(allignment_name_file)
tree <- build_tree(alligned,tree_nm,multi = T,boot = F) # the heaviest part of the code

#showing proper tree
#tree <- read.tree("3z_third_Myodes_080321_nm_DQB_tree.new")
#tree_image <- "3z_third_Myodes_080321_nm_DQB_tree.png"
tree_image <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree.png",sep = "")
create_gg_tree(tree,with_stop_codons,frame_shi = frame_shifts,
               image_file = tree_image,x_lim = 1,height_im = 25)
# ^ x_lim - size of the x axis - often needs to be adjusted depending on the branches length
# ^ height_im - height of the imag  e - needs to be increased when there is many alleles

### remove non-classical
stops <- with_stop_codons # seqs with stop codons

#to_remove <- get_alleles_names("Ang_MHCIIex2",c(250,199,240,234,57,193,246)) # this if additional removed based on tree
# ^ it can happen that e.g. one proper is clustered with stop codons, you can consider removing it
#to_remove <- c("Tri_MHCIIex2_276","Tri_MHCIIex2_271") # this if nothing else removed
to_remove <- c()

#one more data modification
data2 <- remove_some_seqs(data,alleles = c(stops,to_remove,frame_shifts),log = log_name,pat_or_vec = "VEC")
out_name <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_modified_withoutSTOPS.csv",sep = "")
write.csv(data2,out_name,quote = F,fileEncoding = "UTF-8",row.names = F)

repeatability(data2,log_name)

#descriptive
write_to_log("DESCRIPTION",body = paste("Number of individuals: ",ncol(data)-9,"\nNumber of sequences: ",
                                        nrow(data),"\nMean number of alleles per individual: ",
                                        round(mean(as.numeric(apply(data[,10:ncol(data)],2,function(x){sum(!is.na(x))}))),2)),log_name)

#removed <- subset(alligned,!(labels(alligned) %in% c(stops,to_remove))) ### old - do not use!
removed <- alligned[!(labels(alligned) %in% c(stops,to_remove,frame_shifts))]

#modified tree
# tree_nm <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree_wihoutSTOPS.new",sep = "")
# tree <- build_tree(removed,tree_nm,multi = T,boot = F) 
# tree_image <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_tree_withoutSTOPS.png",sep = "")
# create_gg_tree(tree,with_stop_codons,frame_shi = frame_shifts,
#                image_file = tree_image,x_lim = 0.3,height_im = 14)

#removed <- alligned # if nothing is removed (above should work for both cases)
allignment_func <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned_functional.fasta",sep = "")
write.FASTA(x = removed,file = allignment_func) #final allignment

#translate to protein
allignment_aa <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned_aa.fasta",sep = "")
protein_alignment(unalligned_DNA = fasta_bin,not_needed = c(),
                  framing = frame_translate,out_file = allignment_aa)
allignment_func_aa <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_alligned_functional_aa.fasta",sep = "")
protein_alignment(unalligned_DNA = fasta_bin,not_needed = c(stops,to_remove,frame_shifts),
                  framing = frame_translate,out_file = allignment_func_aa) #protein allignment

#change to matrix with genotypes
#data <- read.csv("Emys_MHCI_ex2_combined_modified.csv",header = T,check.names = F) ### backdoor
gen_file <- paste(strsplit(analysed_file,split = "\\.")[[1]][1],"_genotypes.csv",sep = "")
genotypes <- change_to_genotypes(data,removed,gen_file) #table with genotypes

