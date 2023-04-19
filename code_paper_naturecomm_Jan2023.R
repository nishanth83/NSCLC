## code to reproduce some of the key analysis used in the Nair et al. paper titled "A landscape of response to drug combinations in non-small cell lung cancer". 
# author: Nishanth Ulhas Nair

library(data.table)
library(ggplot2)
library(rio)
library(openxlsx)


load("for_keeping_only_zone2_and_3_correct_June2019.RData")
load("store_data_code.RData")

###### figure 1 is an overview figure and does not have any code or results


######  analysis related to figure 2

originalscreen <- readRDS("original_screen_data.RDS")


originalscreen_AZD7762 <- originalscreen[which(originalscreen$Anchor_Drug_Name == "AZD7762"),] # data for combinations involving AZD7762
originalscreen_AZD7762_HSAvalues <- originalscreen_AZD7762[,c(1:3,31:36)] # HSA values for combinations involving AZD7762
## export(originalscreen_AZD7762_HSAvalues, file="originalscreen_AZD7762_HSAvalues_store.txt")

originalscreen_AZD7762_HSAvalues_Adavosertib <- originalscreen_AZD7762_HSAvalues[which(originalscreen_AZD7762_HSAvalues$Library_Drug_Name == "Adavosertib"),] # this is an example for getting the HSA data for a single drug combination involving AZD7762 and Adavosertib across all cell lines. 
## export(originalscreen_AZD7762_HSAvalues_Adavosertib, file="originalscreen_AZD7762_HSAvalues_Adavosertib_store.txt")

originalscreen_AZD7762_HSAvalues_Adavosertib_NCIH460 <- originalscreen_AZD7762_HSAvalues_Adavosertib[which(originalscreen_AZD7762_HSAvalues_Adavosertib$Cell_Line == "NCI-H460"),] # this is an example for getting the HSA data for a single drug combination involving AZD7762 and Adavosertib for one cell line NCI-H460. 


originalscreen_AZD7762_synergy_values <- originalscreen_AZD7762[,c(1:3,19:24)] # Synergy values for combinations involving AZD7762
## export(originalscreen_AZD7762_synergy_values, file="originalscreen_AZD7762_synergy_values_store.txt")

originalscreen_AZD7762_synergy_values_Adavosertib <- originalscreen_AZD7762_synergy_values[which(originalscreen_AZD7762_synergy_values$Library_Drug_Name == "Adavosertib"),] # this is an example for getting the synergy data for a single drug combination involving AZD7762 and Adavosertib across all cell lines. 
## export(originalscreen_AZD7762_synergy_values_Adavosertib, file="originalscreen_AZD7762_synergy_values_Adavosertib_store.txt")

originalscreen_AZD7762_synergy_values_Adavosertib_NCIH460 <- originalscreen_AZD7762_synergy_values_Adavosertib[which(originalscreen_AZD7762_synergy_values_Adavosertib$Cell_Line == "NCI-H460"),] # this is an example for getting the synergy data for a single drug combination involving AZD7762 and Adavosertib for one cell line NCI-H460. 

######  analysis related to figure 3

# compute for each drug combination the percentage of highly synergistic cell lines (top 5% chosen as the threshold of high synergy)
threshold_synergy = quantile(originalscreen$synergy_second_best_ratio, na.rm=TRUE, prob=0.05)
store_synergycombos_percent = c()
store_anchor = c()
store_lib = c()
for (aa in 1:length(uaid_keep)) {
	for (ll in 1:length(ulid_keep)) {
		synviab_latest_zone23_keep_percombo = originalscreen[which((uaid_keep[aa]==originalscreen$AnchorID) & (ulid_keep[ll]==originalscreen$LibraryID)),]
		selcid = synviab_latest_zone23_keep_percombo$CellID
		syntemp = synviab_latest_zone23_keep_percombo$synergy_second_best_ratio[match(ucid_keep, selcid)]
		store_anchor = c(store_anchor, uaid_keep[aa])
		store_lib = c(store_lib, ulid_keep[ll])
		store_synergycombos_percent = c(store_synergycombos_percent, ((sum(syntemp < threshold_synergy, na.rm=TRUE)/sum(!is.na(syntemp)))*100))
	}
}

store_anchorname = unq_anchor_name[match(store_anchor, uaid_keep)]
store_libraryname = unq_library_name[match(store_lib, ulid_keep)]

dat = data.frame(cbind(store_anchor, store_lib, store_anchorname, store_libraryname, store_synergycombos_percent))
dat$store_synergycombos_percent = as.numeric(as.character(dat$store_synergycombos_percent))
dat1 = dat[order(dat$store_synergycombos_percent, decreasing=TRUE), ]
## saveRDS(dat1, file= "drugcombinations_percentage_highsynergy.RDS")

synorder <- readRDS("drugcombinations_percentage_highsynergy.RDS") ## file containing the percentage of highly synergistic cell lines in each drug combination. 

synorder_AZD7762 <- synorder[which(synorder$store_anchorname == "AZD7762"),]
synorder_Olaparib <- synorder[which(synorder$store_anchorname == "Olaparib"),]
synorder_Pemetrexed <- synorder[which(synorder$store_anchorname == "Pemetrexed"),]
synorder_Trametinib <- synorder[which(synorder$store_anchorname == "Trametinib"),]
synorder_Alpelisib <- synorder[which(synorder$store_anchorname == "Alpelisib"),]
synorder_Navitoclax <- synorder[which(synorder$store_anchorname == "Navitoclax"),]

## export(synorder_AZD7762, file="percentage_of_high_synergistic_cell_lines_for_drug_combination_anchor_AZD7762.txt")
## export(synorder_Olaparib, file="percentage_of_high_synergistic_cell_lines_for_drug_combination_anchor_Olaparib.txt")
## export(synorder_Pemetrexed, file="percentage_of_high_synergistic_cell_lines_for_drug_combination_anchor_Pemetrexed.txt")
## export(synorder_Trametinib, file="percentage_of_high_synergistic_cell_lines_for_drug_combination_anchor_Trametinib.txt")
## export(synorder_Alpelisib, file="percentage_of_high_synergistic_cell_lines_for_drug_combination_anchor_Alpelisib.txt")
## export(synorder_Navitoclax, file="percentage_of_high_synergistic_cell_lines_for_drug_combination_anchor_Navitoclax.txt")



###### analysis related to figure 4


originalscreen_fig4 <- readRDS("original_screen_data.RDS")
originalscreen_fig4$log_synergy_ratio_D1 = log(originalscreen_fig4$synergy_ratio_D1)
originalscreen_fig4$log_synergy_ratio_D2 = log(originalscreen_fig4$synergy_ratio_D2)
originalscreen_fig4$log_synergy_ratio_D3 = log(originalscreen_fig4$synergy_ratio_D3)
originalscreen_fig4$log_synergy_ratio_D4 = log(originalscreen_fig4$synergy_ratio_D4)
originalscreen_fig4$log_synergy_ratio_D5 = log(originalscreen_fig4$synergy_ratio_D5)
originalscreen_fig4$log_synergy_second_best_ratio = log(originalscreen_fig4$synergy_second_best_ratio)


originalscreen_Trametinib <- originalscreen_fig4[which(originalscreen_fig4$Anchor_Drug_Name == "Trametinib"),] # data for combinations involving Trametinib
originalscreen_Trametinib_synergy_values <- originalscreen_Trametinib[,c(1:3,19:24)] # Synergy values for combinations involving Trametinib
originalscreen_Trametinib_synergy_values_Imatinib = originalscreen_Trametinib_synergy_values[which(originalscreen_Trametinib_synergy_values$Library_Drug_Name == "Imatinib"),]
originalscreen_Trametinib_synergy_values_Nilotinib = originalscreen_Trametinib_synergy_values[which(originalscreen_Trametinib_synergy_values$Library_Drug_Name == "Nilotinib"),]
# similar code as above can be used for any other drug combination used in Fig. 4 analysis

## export(originalscreen_Trametinib_synergy_values_Imatinib, file="originalscreen_Trametinib_synergy_values_Imatinib_log_store.txt") ## stores synergy values for anchor drug Trametinib and library drug Imatinib
## export(originalscreen_Trametinib_synergy_values_Nilotinib, file="originalscreen_Trametinib_synergy_values_Nilotinib_log_store.txt") ## stores synergy values for anchor drug Trametinib and library drug Nilotinib


synorder_Tozasertib_RO3306 <- synorder[which((synorder$store_anchorname == "Tozasertib") & (synorder$store_libraryname == "RO-3306")),] ## file containing the percentage of highly synergistic cell lines in each drug combination Tozasertib & RO-3306. Similarly you can choose other drug combinations. 



###### analysis related to figure 5


load("best_matrix_validationdata_cancers_Jan2023.RData") ## validation screen on cancer cell lines
# best_matrix_DBSumNeg can be used for synergy scores (best score out of 4 replicates) for each drug combination and cell line (rows are cell lines and columns are drug combinations). This can be used to create main figure 5. 
# best_matrix_ExcessHSA can be used for HSA scores (best score out of 4 replicates) for each drug combination and cell line (rows are cell lines and columns are drug combinations). This can be used to create main figure 5.



load("originaldata_subset_for_validation_Jan2023.RData") ## subset of original screen for the purpose of validation. We use only second best synergy score (out of 5 library doses) for each drug combination for each cell line. 


thr1 = 0.2
org_threshold_synergy = quantile(original_secondbestsynergy, na.rm=TRUE, prob=thr1)
org_threshold_hsa = quantile(original_secondbestHSA, na.rm=TRUE, prob=thr1)
valthr1 = -1000

rank_syergy_org_val = matrix(NA,dim(original_secondbestsynergy)[2],2)
rank_hsa_org_val = matrix(NA,dim(original_secondbestsynergy)[2],2)
for (dd in 1:dim(original_secondbestsynergy)[2]) {
	ord_rank = sum(original_secondbestsynergy[,dd] < org_threshold_synergy, na.rm=T)/sum(!is.na(original_secondbestsynergy[,dd]))
	val_rank = sum(best_matrix_ExcessHSA[,dd] < valthr1, na.rm=T)/sum(!is.na(best_matrix_ExcessHSA[,dd]))
	rank_syergy_org_val[dd,] = c(ord_rank, val_rank)

	ord_rank = sum(original_secondbestHSA[,dd] < org_threshold_hsa, na.rm=T)/sum(!is.na(original_secondbestHSA[,dd]))
	rank_hsa_org_val[dd,] = c(ord_rank, val_rank)

}
cor.test(rank_syergy_org_val[,1], rank_syergy_org_val[,2], method="spearman") ## correlation of percentage of highly synergistic cell lines between the drug combinations in the original screen and validation screen


cors= cor.test(rank_syergy_org_val[,1]*100, rank_syergy_org_val[,2]*100, method="spearman")


new_unq_combid_few = matrix(NA, length(colnames(best_matrix_ExcessHSA)), 1)
new_unq_combid_few[which((rank_syergy_org_val[,1] >= 0.10) & (rank_syergy_org_val[,2] >= 0.10))] = colnames(best_matrix_ExcessHSA)[which((rank_syergy_org_val[,1] >= 0.10) & (rank_syergy_org_val[,2] >= 0.10))]


df1 = data.frame(original=rank_syergy_org_val[,1]*100, validation=rank_syergy_org_val[,2]*100, combo = new_unq_combid_few)
 
pdf("scatter_synergy_compare_original_validation_screen.pdf", width=10,height=10)  ## this figure generates Supp. figure 10. 
p1 <- ggplot(df1, aes(x=original, y=validation)) + geom_point(shape=2,size=3, color="red") + geom_smooth(method=lm) + theme_bw() + theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=25), axis.title=element_text(size=25), plot.title = element_text(size = 25, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none')  + ylab("Drug combinations (validation screen) \n (% of highly synergistic cell lines)") + xlab("Drug combinations (original screen) \n (% of highly synergistic cell lines)") + ggtitle("Correlation of high synergies between \n original screen and validation screen") + annotate("text", x = c(30,30), y = c(35,30), label = c(paste("rho =",format(cors$estimate,digits=2)), paste("P =",format(cors$p.value,digits=3))), size = 9, angle=0) 
p1 <- p1 + ggrepel::geom_text_repel(aes(angle =0,  hjust=-0.6,  label=combo), size=5, point.padding = 1, force=1, color="blue") 
p1
dev.off()



###### analysis related to figure 6


load("store_another_code.RData")

anchordrugid = unique(anchordrugfile$"Compound#")
relevant_anchorgenes = anchordrugfile$"targets_latest"[match(anchordrugid,anchordrugfile$"Compound#")]
anchor_drugs_targets_map = data.frame(cbind(relevant_anchorgenes, anchordrugid))

librarydrugid = unique(librarydrugfile$"CMT_CPD#")
relevant_librarygenes = librarydrugfile$"TARGETGENESLatest"[match(librarydrugid,librarydrugfile$"CMT_CPD#")]


library_drugs_targets_map = data.frame(cbind(relevant_librarygenes, librarydrugid))

all_libgenes = unique(unlist(strsplit(relevant_librarygenes,",")))
all_anchgenes = unique(unlist(strsplit(relevant_anchorgenes,",")))


### for figure 6A

target_pairs_per_como = matrix(NA,1,dim(syn_per_combo_order)[1])
targets_per_como = matrix(NA,1,dim(syn_per_combo_order)[1])
for (i in 1:dim(syn_per_combo_order)[1]){
	anchdrugclassval = drugclass[match(syn_per_combo_order$AnchorID[i], drugclass$"CMT_CPD#")[1],6:14]
	anchdrugclassval_names = names(anchdrugclassval)[which(!is.na(anchdrugclassval))]
	libdrugclassval = drugclass[match(syn_per_combo_order$LibraryID[i], drugclass$"CMT_CPD#")[1],6:14]
	libdrugclassval_names = names(libdrugclassval)[which(!is.na(libdrugclassval))]
	
	if (!((sum(anchdrugclassval_names == "Chemotherapeutics") > 0) | (sum(libdrugclassval_names == "Chemotherapeutics") > 0))) {

		anchtarget = anchor_drugs_targets_map$relevant_anchorgenes[match(syn_per_combo_order$AnchorID[i], anchor_drugs_targets_map$anchordrugid)]
		libtarget = library_drugs_targets_map$relevant_librarygenes[match(syn_per_combo_order$LibraryID[i], library_drugs_targets_map$librarydrugid)]
		anchtargetunlist = unlist(strsplit(as.character(anchtarget), ","))
		anchtargetunlist = anchtargetunlist[anchtargetunlist %in% humangenename$Genename]
		libtargetunlist = unlist(strsplit(as.character(libtarget), ","))
		libtargetunlist = libtargetunlist[libtargetunlist %in% humangenename$Genename]
		if (length(libtargetunlist) > 0) {
			targetpairs1 = cbind(rep(anchtargetunlist,each=length(libtargetunlist)),rep(libtargetunlist,length(anchtargetunlist)))
			if (length(targetpairs1) > 2){
				targetpairs = targetpairs1[targetpairs1[,1]!=targetpairs1[,2],]
				targets_per_como[i] = length(unique(c(targetpairs)))
				if (length(targetpairs)>2) {
					target_pairs_per_como[i] = dim(targetpairs)[1]
					} else {
						target_pairs_per_como[i] = 1
					}
				} else {
					if (sum(targetpairs1[1]!=targetpairs1[2]) > 0) {
						targets_per_como[i] = length(unique(c(targetpairs1)))
						target_pairs_per_como[i] = 1
					} else {
						target_pairs_per_como[i] = NA
					}
				}
			}
		}
}


cors = cor.test(t(targets_per_como),  syn_per_combo_order$percentage_synergistic_celllines, method="spearman")

df = data.frame(targets=t(targets_per_como), persyn=syn_per_combo_order$percentage_synergistic_cellline)
 
pdf("newfigure_scatter_synergycombos_totaltargets_correctdata_Aug2019_UPDATED_Feb2020.pdf", width=8,height=8) 
p1 <- ggplot(df, aes(x=targets, y=persyn)) + geom_point(shape=1, color="blue") + geom_smooth(method=lm) +  theme_bw() + theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30), axis.title=element_text(size=30), plot.title = element_text(size = 35, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + ylab("Percentage of synergistic cell lines") + xlab("No. of total targets") + annotate("text", x = c(15,15), y = c(55,50), label = c(paste("rho =",format(cors$estimate,digits=2)),paste("P =",format(cors$p.value,digits=3))), size = 10) 
p1
dev.off()

## for figure 6D


hsasyn_385 = hsasyn[which(hsasyn$AnchorID == 385),]

thresholdsvarious = c(0.1,0.2,0.5,0.6,0.7,0.732,0.8,0.9,1)

thresholdsvarious_cellcount = array(NA,c(length(thresholdsvarious),length(ulid_keep),5))
for (i in 1:length(thresholdsvarious)) {
	cn = 0 
	# for (aa in 1:length(uaid_keep)) {
		for (ll in 1:length(ulid_keep)) {
			synviab_latest_zone23_keep_ll = hsasyn_385[which((hsasyn_385$LibraryID == ulid_keep[ll])),]
			if ((dim(synviab_latest_zone23_keep_ll)[1]) > 0) {
				cn = cn + dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,1] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D1 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,2] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D2 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,3] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D3 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,4] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D4 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,5] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D5 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
			}
		}
	# }
}


dose_vec = c("D1", "D2", "D3", "D4", "D5")
datdf_store = c()
for (i in 1:dim(thresholdsvarious_cellcount)[1]) {
	thresholdsvarious_cellcount_tt = thresholdsvarious_cellcount[i,,]
	datdf = data.frame(cbind(count=c(t(thresholdsvarious_cellcount_tt)), thresh=rep(thresholdsvarious[i], each=length(thresholdsvarious_cellcount_tt)), dose=rep(dose_vec, dim(thresholdsvarious_cellcount_tt)[1])))
	datdf$count = as.numeric(datdf$count)
	datdf$thresh = as.numeric(datdf$thresh)
	datdf_store = rbind(datdf_store, datdf)
}


library(RColorBrewer)
myColors <- brewer.pal(9,"Paired")
names(myColors) <- levels(thresholdsvarious)
colScale <- scale_colour_manual(name = "grp",values = myColors)


datdf_store$thresh = factor(datdf_store$thresh, levels=sort(unique(datdf_store$thresh), decreasing=T))

pdf("figure_original_different_synergy_anchor_385_thresholds_countcellline_density_Jan2021.pdf", width=12,height=10) ## for Navitoclax
p1 <- ggplot(datdf_store, aes(x=count*100, after_stat(density), colour=thresh)) +  geom_freqpoly(size= 2, bins = 25) + theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), axis.title.x=element_text(size=40), axis.title.y=element_text(size=40), plot.title = element_text(size = 40), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=30), legend.text=element_text(size=30), legend.position = c(0.7, 0.6)) + ylab("Density") + xlab("Percentage of total screened")  + coord_cartesian(xlim = c(0, 100), expand=FALSE) + guides(colour=guide_legend(title="Threshold")) +  scale_colour_manual(values = c("purple", "darkblue", "skyblue", "orange", "grey", "black", "red", "green", "blue"))
p1
dev.off()


### for figure 6E


thresholdsvarious = c(0.1,0.2,0.5,0.6,0.7,0.732,0.8,0.9,1)

thresholdsvarious_cellcount = array(NA,c(length(thresholdsvarious),length(ulid_keep),5))
for (i in 1:length(thresholdsvarious)) {
	cn = 0 
	# for (aa in 1:length(uaid_keep)) {
		for (ll in 1:length(ulid_keep)) {
			synviab_latest_zone23_keep_ll = hsasyn[which((hsasyn$LibraryID == ulid_keep[ll])),]
			if ((dim(synviab_latest_zone23_keep_ll)[1]) > 0) {
				cn = cn + dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,1] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D1 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,2] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D2 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,3] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D3 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,4] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D4 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,5] = sum(synviab_latest_zone23_keep_ll$synergy_ratio_D5 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
			}
		}
	# }
}


dose_vec = c("D1", "D2", "D3", "D4", "D5")
datdf_store = c()
for (i in 1:dim(thresholdsvarious_cellcount)[1]) {
	thresholdsvarious_cellcount_tt = thresholdsvarious_cellcount[i,,]
	datdf = data.frame(cbind(count=c(t(thresholdsvarious_cellcount_tt)), thresh=rep(thresholdsvarious[i], each=length(thresholdsvarious_cellcount_tt)), dose=rep(dose_vec, dim(thresholdsvarious_cellcount_tt)[1])))
	datdf$count = as.numeric(datdf$count)
	datdf$thresh = as.numeric(datdf$thresh)
	datdf_store = rbind(datdf_store, datdf)
}


library(RColorBrewer)
myColors <- brewer.pal(9,"Paired")
names(myColors) <- levels(thresholdsvarious)
colScale <- scale_colour_manual(name = "grp",values = myColors)


datdf_store$thresh = factor(datdf_store$thresh, levels=sort(unique(datdf_store$thresh), decreasing=T))

pdf("figure_original_different_synergy_thresholds_countcellline_density_Jan2021.pdf", width=12,height=10)
p1 <- ggplot(datdf_store, aes(x=count*100, after_stat(density), colour=thresh)) +  geom_freqpoly(size= 2, bins = 25) + theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), axis.title.x=element_text(size=40), axis.title.y=element_text(size=40), plot.title = element_text(size = 40), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=30), legend.text=element_text(size=30), legend.position = c(0.7, 0.6)) + ylab("Density") + xlab("Percentage of total screened")  + coord_cartesian(xlim = c(0, 100), expand=FALSE) + guides(colour=guide_legend(title="Threshold")) +  scale_colour_manual(values = c("purple", "darkblue", "skyblue", "orange", "grey", "black", "red", "green", "blue"))
p1
dev.off()


## figure 6F
thresholdsvarious = c(-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0)

thresholdsvarious_cellcount = array(NA,c(length(thresholdsvarious),length(ulid_keep),5))
for (i in 1:length(thresholdsvarious)) {
	cn = 0 
	# for (aa in 1:length(uaid_keep)) {
		for (ll in 1:length(ulid_keep)) {
			synviab_latest_zone23_keep_ll = hsasyn[which((hsasyn$LibraryID == ulid_keep[ll])),]
			if ((dim(synviab_latest_zone23_keep_ll)[1]) > 0) {
				cn = cn + dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,1] = sum(synviab_latest_zone23_keep_ll$HSA_D1 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,2] = sum(synviab_latest_zone23_keep_ll$HSA_D2 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,3] = sum(synviab_latest_zone23_keep_ll$HSA_D3 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,4] = sum(synviab_latest_zone23_keep_ll$HSA_D4 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
				thresholdsvarious_cellcount[i,ll,5] = sum(synviab_latest_zone23_keep_ll$HSA_D5 < thresholdsvarious[i], na.rm=T) / dim(synviab_latest_zone23_keep_ll)[1]
			}
		}
	# }
}


dose_vec = c("D1", "D2", "D3", "D4", "D5")
datdf_store = c()
for (i in 1:dim(thresholdsvarious_cellcount)[1]) {
	thresholdsvarious_cellcount_tt = thresholdsvarious_cellcount[i,,]
	datdf = data.frame(cbind(count=c(t(thresholdsvarious_cellcount_tt)), thresh=rep(thresholdsvarious[i], each=length(thresholdsvarious_cellcount_tt)), dose=rep(dose_vec, dim(thresholdsvarious_cellcount_tt)[1])))
	datdf$count = as.numeric(datdf$count)
	datdf$thresh = as.numeric(datdf$thresh)
	datdf_store = rbind(datdf_store, datdf)
}


datdf_store$thresh = factor(datdf_store$thresh, levels=sort(unique(datdf_store$thresh), decreasing=T))

pdf("figure_original_different_hsa_thresholds_countcellline_density_Jan2021.pdf", width=12,height=10)
p1 <- ggplot(datdf_store, aes(x=count*100, after_stat(density), colour=thresh)) +  geom_freqpoly(size= 2, bins = 25) + theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), axis.title.x=element_text(size=40), axis.title.y=element_text(size=40), plot.title = element_text(size = 40), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=30), legend.text=element_text(size=30), legend.position = c(0.7, 0.6)) + ylab("Density") + xlab("Percentage of total screened")  + coord_cartesian(xlim = c(0, 100), expand=FALSE) + guides(colour=guide_legend(title="Threshold")) +  scale_colour_manual(values = c("purple", "darkblue", "skyblue", "orange", "grey", "black", "red", "green", "blue"))
p1
dev.off()


## figure 6G



greedyfunction_minsetcover <- function(mat1, minthreshpercent, cidsel) {
	store_maxcellline_ids = c()
	flag = 0
	cn_drug = 0
	alldrugscount = dim(mat1)[1]
	allcellcount = dim(mat1)[2]
	while (flag == 0) {
		cn_drug = cn_drug + 1
		celllinedistr = apply(mat1,1,sum,na.rm=TRUE)
		maxcellline_num = max(celllinedistr,na.rm=TRUE)
		index_maxdrug = which.max(celllinedistr)
		relevant_cid = cidsel[which(mat1[index_maxdrug,] == 1)]
		relevant_cid_index = which(mat1[index_maxdrug,] == 1)
		store_maxcellline_ids = union(store_maxcellline_ids, relevant_cid)
		percent_ratio_celllines = (length(store_maxcellline_ids)/allcellcount)*100
		if ((percent_ratio_celllines >=  minthreshpercent) | (cn_drug >= alldrugscount)) {
			flag = 1
		}
		cidsel = cidsel[setdiff(1:dim(mat1)[2],relevant_cid_index)]
		mat1 = mat1[,setdiff(1:dim(mat1)[2],relevant_cid_index)]
		if (!is.matrix(mat1)){
			mat1 = as.matrix(t(mat1), nrow=1)
		}
	}
	return(c(cn_drug, percent_ratio_celllines))
}

threshold_synergy = quantile(synviab_latest_zone23_keep$synergy_second_best_ratio, na.rm=TRUE, prob=0.05)
store_setcover_greedy_allanchors = array(NA,c(length(uaid_keep),11,2))
for (aa in 1:length(uaid_keep)) {

	print(aa)
	synviab_anch = synviab_latest_zone23_keep[which(synviab_latest_zone23_keep$AnchorID == uaid_keep[aa]),]

	synviab_anch$highsynergy = as.numeric(synviab_anch$synergy_second_best_ratio < threshold_synergy)

	cid_highsyn = synviab_anch$CellID[which(synviab_anch$highsynergy == 1)]
	unique_cid_highsyn = unique(cid_highsyn)
	unique_allcid = unique(synviab_anch$CellID)
	# length(unique_cid_highsyn) / length(unique_allcid)

	# paste(synviab_anch$LibraryID, synviab_anch$CellID, sep="_")

	synviab_anch_lidcid = matrix(NA,length(ulid_keep),length(ucid_keep))
	for (ll in 1:length(ulid_keep)) {
		for (cc in 1:length(ucid_keep)) {
			if (length(which((synviab_anch$LibraryID == ulid_keep[ll]) & (synviab_anch$CellID == ucid_keep[cc]))) > 0) {
				synviab_anch_lidcid[ll,cc] = as.numeric(synviab_anch$synergy_second_best_ratio[which((synviab_anch$LibraryID == ulid_keep[ll]) & (synviab_anch$CellID == ucid_keep[cc]))] < threshold_synergy)
			}
		}
	}
	ucid_keep_selected = ucid_keep[colSums(is.na(synviab_anch_lidcid))<nrow(synviab_anch_lidcid)]
	synviab_anch_lidcid <- synviab_anch_lidcid[,colSums(is.na(synviab_anch_lidcid))<nrow(synviab_anch_lidcid)]
	store_setcover_greedy = matrix(NA,11,2)
	tt1 = 1
	for (thresh_minpercent in c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100)){
		setcover_greedy <- greedyfunction_minsetcover(synviab_anch_lidcid, thresh_minpercent, ucid_keep_selected)
		if ((setcover_greedy[1]==dim(synviab_anch_lidcid)[1]) & (setcover_greedy[2]<thresh_minpercent)) {
			setcover_greedy = c(NA, setcover_greedy[2])
		}
		store_setcover_greedy[tt1,] = setcover_greedy
		tt1 = tt1 + 1
	}
	store_setcover_greedy_allanchors[aa,,] = store_setcover_greedy
}


numbers_store_setcover_greedy_allanchors = store_setcover_greedy_allanchors[,,1]
datf = data.frame(cbind(uaid_keep, unq_anchor_name, numbers_store_setcover_greedy_allanchors))
colnames(datf) = c("anchor_ID", "anchor_name", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", "100")

# write(colnames(datf), file="minimum_no_of_librarydrugs_for_each_anchor_for_percentage_of_highsynergycelllines_Feb2020.txt", ncol=length(colnames(datf)), append=FALSE, sep="\t")
# write(t(datf), file="minimum_no_of_librarydrugs_for_each_anchor_for_percentage_of_highsynergycelllines_Feb2020.txt", ncol=length(colnames(datf)), append=TRUE, sep="\t")

dt1 = fread("minimum_no_of_librarydrugs_for_each_anchor_for_percentage_of_highsynergycelllines_Feb2020.txt", header=TRUE)
dt1 = dt1[order(dt1$"80",decreasing=FALSE),]
dt1$anchor_name = factor(dt1$anchor_name, levels = unique(dt1$anchor_name))

dt1$consider = dt1$"80"

pdf("newfigure_minimum_no_of_librarydrugs_names_for_each_anchor_for_percentage_of_highsynergycelllines_UPDATED_Feb2020.pdf", width=12,height=10)
p1 <- ggplot(dt1, aes(x=dt1$anchor_name, y=consider, fill=anchor_name)) + geom_bar(stat="identity", position=position_dodge(), width = 0.5) +  theme_bw() + theme(axis.text.x=element_text(size=30, angle=90, hjust = 0.95, vjust = 0.5), axis.text.y=element_text(size=30), axis.title=element_text(size=30), plot.title = element_text(size = 35, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none')  + xlab("Anchor drugs") + ylab("Minimum no. of library drugs to \n cover at least 80% of cell lines") 
p1
dev.off()



dt1 = fread("minimum_no_of_librarydrugs_for_each_anchor_for_percentage_of_highsynergycelllines_Feb2020.txt", header=TRUE)
dt1 = dt1[order(dt1$"80",decreasing=FALSE),]
dt1$anchor_ID = factor(dt1$anchor_ID, levels = unique(dt1$anchor_ID))

dt1$consider = dt1$"80"

pdf("newfigure_minimum_no_of_librarydrugs_ID_for_each_anchor_for_percentage_of_highsynergycelllines_UPDATED_Feb2020.pdf", width=12,height=10)
p1 <- ggplot(dt1, aes(x=dt1$anchor_ID, y=consider, fill=anchor_name)) + geom_bar(stat="identity", position=position_dodge(), width = 0.5) +  theme_bw() + theme(axis.text.x=element_text(size=30, angle=90, hjust = 0.95, vjust = 0.5), axis.text.y=element_text(size=30), axis.title=element_text(size=30), plot.title = element_text(size = 35, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none')  + xlab("Anchor drugs") + ylab("Minimum no. of library drugs to \n cover at least 80% of cell lines") 
p1
dev.off()


############ analysis related to figure 7


## super sensitive cell line analysis for figure 7A

load("viabzone23_leep_anchlibonly.RData")

best_resistant_celllines_minus_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
best_sensitive_celllines_minus_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
best_all_celllines_minus_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
mean_all_celllines_minus_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
resistant_celllines_compared_to_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
sensitive_celllines_compared_to_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
all_celllines_compared_to_most_sensitive = matrix(NA,length(uaid_keep),length(ulid_keep))
best_across_viab = matrix(NA,length(uaid_keep),length(ulid_keep))
best_across_anch = matrix(NA,length(uaid_keep),1)
for (aa in 1:length(uaid_keep)) {
	for (ll in 1:length(ulid_keep)) {
		synviab_latest_zone23_keep_percombo = synviab_latest_zone23_keep[which((uaid_keep[aa]==synviab_latest_zone23_keep$AnchorID) & (ulid_keep[ll]==synviab_latest_zone23_keep$LibraryID)),]
		selcid = synviab_latest_zone23_keep_percombo$CellID
		viab_AL_D4 = synviab_latest_zone23_keep_percombo$AL_by_OO_D4[match(ucid_keep, selcid)]
	
		anchviabtemp = viabanchonly_latest_zone23_keep[which(uaid_keep[aa]==viabanchonly_latest_zone23_keep$AnchorID),]
		anchviab = anchviabtemp$AO_by_OO_median[match(ucid_keep, anchviabtemp$CellID)]

		libviabtemp = viablibonly_latest_zone23_keep[which(ulid_keep[ll]==viablibonly_latest_zone23_keep$LibraryID),]
		libviab = libviabtemp$OL_by_OO_D4_median[match(ucid_keep, libviabtemp$CellID)]

		min_of_two_viab = pmin(anchviab,libviab)
		bestviab = min(anchviab,libviab, na.rm=T)
		bestviab_cellid = ucid_keep[which(min_of_two_viab==bestviab)]

		best_across_viab[aa,ll] = bestviab

		resistant_celllines = ucid_keep[order(min_of_two_viab, decreasing=TRUE)[1:10]]
		sensitive_celllines = ucid_keep[order(min_of_two_viab, decreasing=FALSE)[1:11]]
		sensitive_celllines = setdiff(sensitive_celllines, bestviab_cellid) ## remove most sensitive cell line
		allcelllines_expect_best = setdiff(ucid_keep, bestviab_cellid)
		# allcelllines_expect_best = ucid_keep

		resistant_celllines_compared_to_most_sensitive[aa,ll] = (sum(viab_AL_D4[match(resistant_celllines, ucid_keep)] < bestviab, na.rm=TRUE)/sum(!is.na(viab_AL_D4[match(resistant_celllines, ucid_keep)])))*100
		sensitive_celllines_compared_to_most_sensitive[aa,ll] = (sum(viab_AL_D4[match(sensitive_celllines, ucid_keep)] < bestviab, na.rm=TRUE)/sum(!is.na(viab_AL_D4[match(sensitive_celllines, ucid_keep)])))*100
		all_celllines_compared_to_most_sensitive[aa,ll] = (sum(viab_AL_D4[match(allcelllines_expect_best, ucid_keep)] < bestviab, na.rm=TRUE)/sum(!is.na(viab_AL_D4[match(allcelllines_expect_best, ucid_keep)])))*100


		best_resistant_celllines_minus_most_sensitive[aa,ll] = (min(viab_AL_D4[match(resistant_celllines, ucid_keep)], na.rm=TRUE) - bestviab)
		best_sensitive_celllines_minus_most_sensitive[aa,ll] = (min(viab_AL_D4[match(sensitive_celllines, ucid_keep)], na.rm=TRUE) - bestviab)
		best_all_celllines_minus_most_sensitive[aa,ll] = (min(viab_AL_D4[match(allcelllines_expect_best, ucid_keep)], na.rm=TRUE) - bestviab)
		viabdifff = (viab_AL_D4[match(allcelllines_expect_best, ucid_keep)] - bestviab)
		mean_all_celllines_minus_most_sensitive[aa,ll] = mean(sort(viabdifff)[1:5], na.rm=TRUE)


	}
	best_across_anch[aa] = min(anchviab, na.rm=T)

}

rownames(best_across_anch) = uaid_keep

resistant_celllines_compared_to_most_sensitive.ecdf = ecdf(resistant_celllines_compared_to_most_sensitive)


df = data.frame(percent=(c(sum(resistant_celllines_compared_to_most_sensitive==0),sum(sensitive_celllines_compared_to_most_sensitive==0),sum(all_celllines_compared_to_most_sensitive==0), sum((resistant_celllines_compared_to_most_sensitive>0) & (resistant_celllines_compared_to_most_sensitive <= 10)),sum((sensitive_celllines_compared_to_most_sensitive>0) & (sensitive_celllines_compared_to_most_sensitive <= 10)),sum((all_celllines_compared_to_most_sensitive>0) & (all_celllines_compared_to_most_sensitive <= 10)), sum(resistant_celllines_compared_to_most_sensitive>10),sum(sensitive_celllines_compared_to_most_sensitive>10),sum(all_celllines_compared_to_most_sensitive>10))/length(all_celllines_compared_to_most_sensitive))*100, cellline=c("resistant cell lines", "sensitive cell lines", "all cell lines", "resistant cell lines", "sensitive cell lines", "all cell lines", "resistant cell lines", "sensitive cell lines", "all cell lines"), count=c("0", "0", "0", "0 to 10", "0 to 10", "0 to 10", "> 10", "> 10", "> 10"))
df$count = factor(df$count, levels = unique(df$count))

pdf("newfigure_histogram_cyril_hypothesis_resistant_sensitive_all_compared_to_most_sensitive_doseD4_correcteddata_Aug2019_UPDATED_Feb2020.pdf", width=12,height=10)

p1 <- ggplot(df, aes(x=count, y=percent, fill=cellline)) + geom_bar(stat="identity", position=position_dodge(), width = 0.5) +  theme_bw() + theme(axis.text.x=element_text(size=30, angle=0, hjust = 0.95, vjust = 0.5), axis.text.y=element_text(size=30), axis.title=element_text(size=30), plot.title = element_text(size = 35, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), legend.title=element_blank(), legend.text=element_text(size=30))  + ylab("Percentage of tested combinations") + xlab("Super sensitive cell lines (%)") + ylim(c(0,100))
p1
dev.off()

## for figure 7B

store_aid_all = c()
store_lid_all = c()
store_all_celllines_compared_to_most_sensitive_all = c()
store_bestindividual_viab_all = c()
store_aid1 = c()
store_lid1 = c()
store_aid1_low = c()
store_lid1_low = c()
store_aid2 = c()
store_lid2 = c()
store_aid3 = c()
store_lid3 = c()
store_aid4 = c()
store_lid4 = c()
store_all_celllines_compared_to_most_sensitive = c()
store_all_celllines_compared_to_most_sensitive_low = c()
store_best_all_celllines_minus_most_sensitive = c()
store3_all_celllines_compared_to_most_sensitive = c()
store3_best_all_celllines_minus_most_sensitive = c()
store_mean_all_celllines_minus_most_sensitive = c()
store_bestindividual_viab1 = c()
store_bestindividual_viab1_low = c()
store_bestindividual_viab2 = c()
store_bestindividual_viab3 = c()
store_bestindividual_viab4 = c()
for (aa in 1:length(uaid_keep)) {
	for (ll in 1:length(ulid_keep)) {
		if (all_celllines_compared_to_most_sensitive[aa,ll] >= 10) {
			store_aid1 = c(store_aid1, uaid_keep[aa])
			store_lid1 = c(store_lid1, ulid_keep[ll])
			store_all_celllines_compared_to_most_sensitive = c(store_all_celllines_compared_to_most_sensitive, all_celllines_compared_to_most_sensitive[aa,ll])
			store_bestindividual_viab1 = c(store_bestindividual_viab1, best_across_viab[aa,ll])
		}
		if (best_all_celllines_minus_most_sensitive[aa,ll] < (-0.1)) {
			store_aid2 = c(store_aid2, uaid_keep[aa])
			store_lid2 = c(store_lid2, ulid_keep[ll])
			store_best_all_celllines_minus_most_sensitive = c(store_best_all_celllines_minus_most_sensitive, best_all_celllines_minus_most_sensitive[aa,ll])
			store_bestindividual_viab2 = c(store_bestindividual_viab2, best_across_viab[aa,ll])
		}
		if (mean_all_celllines_minus_most_sensitive[aa,ll] < (-0.1)) {
			store_aid4 = c(store_aid4, uaid_keep[aa])
			store_lid4 = c(store_lid4, ulid_keep[ll])
			store_mean_all_celllines_minus_most_sensitive = c(store_mean_all_celllines_minus_most_sensitive, mean_all_celllines_minus_most_sensitive[aa,ll])
			store_bestindividual_viab4 = c(store_bestindividual_viab4, best_across_viab[aa,ll])
		}		
		if ((all_celllines_compared_to_most_sensitive[aa,ll] >= 5) & (best_all_celllines_minus_most_sensitive[aa,ll] < (-0.1))) {
			store_aid3 = c(store_aid3, uaid_keep[aa])
			store_lid3 = c(store_lid3, ulid_keep[ll])
			store3_all_celllines_compared_to_most_sensitive = c(store3_all_celllines_compared_to_most_sensitive, all_celllines_compared_to_most_sensitive[aa,ll])
			store3_best_all_celllines_minus_most_sensitive = c(store3_best_all_celllines_minus_most_sensitive, best_all_celllines_minus_most_sensitive[aa,ll])
			store_bestindividual_viab3 = c(store_bestindividual_viab3, best_across_viab[aa,ll])
		}
		if (all_celllines_compared_to_most_sensitive[aa,ll] == 0) {
			store_aid1_low = c(store_aid1_low, uaid_keep[aa])
			store_lid1_low = c(store_lid1_low, ulid_keep[ll])
			store_all_celllines_compared_to_most_sensitive_low = c(store_all_celllines_compared_to_most_sensitive_low, all_celllines_compared_to_most_sensitive[aa,ll])
			store_bestindividual_viab1_low = c(store_bestindividual_viab1_low, best_across_viab[aa,ll])
		}
		store_aid_all = c(store_aid_all, uaid_keep[aa])
		store_lid_all = c(store_lid_all, ulid_keep[ll])
		store_all_celllines_compared_to_most_sensitive_all = c(store_all_celllines_compared_to_most_sensitive_all, all_celllines_compared_to_most_sensitive[aa,ll])
		store_bestindividual_viab_all = c(store_bestindividual_viab_all, best_across_viab[aa,ll])
		
		
	}
}

drugnames_anchor_all = unq_anchor_name[match(store_aid_all, uaid_keep)]
drugnames_lid_all = unq_library_name[match(store_lid_all, ulid_keep)]
df_write_all = data.frame(cbind(anchorid=store_aid_all, libraryid=store_lid_all, anchorname=drugnames_anchor_all, libraryname=drugnames_lid_all, percent_cellines=store_all_celllines_compared_to_most_sensitive_all, bestindividual_viab = store_bestindividual_viab_all))
df_write_all$percent_cellines = as.numeric(as.character(df_write_all$percent_cellines))
df_write_all_sort = df_write_all[order(df_write_all$percent_cellines, decreasing=TRUE), ]

# write(colnames(df_write_all_sort), file="combinations_supersensitive_list_greaterthanbestindividual_all_combos_correcteddata_Aug2019.txt", append=FALSE, sep="\t", ncolumns=6)
# write(t(df_write_all_sort), file="combinations_supersensitive_list_greaterthanbestindividual_all_combos_correcteddata_Aug2019.txt", append=TRUE, sep="\t", ncolumns=6)

datf = fread("combinations_supersensitive_list_greaterthanbestindividual_all_combos_correcteddata_Aug2019.txt", header=TRUE)
datf$anchorname = new_anchor_name_June2020$AnchorName[match(datf$anchorid, new_anchor_name_June2020$AnchorID)]
datf$libraryname = new_library_name_June2020$Drug.Name[match(datf$libraryid, new_library_name_June2020$"Cmt.Cpd#")]
## write(names(datf), file="combinations_July2020_supersensitive_list_greaterthanbestindividual_all_combos_correcteddata_Aug2019.txt", ncolumns=dim(datf)[2], sep="\t", append=FALSE)
## write(t(datf), file="combinations_July2020_supersensitive_list_greaterthanbestindividual_all_combos_correcteddata_Aug2019.txt", ncolumns=dim(datf)[2], sep="\t", append=TRUE)



drugnames_anchor1 = unq_anchor_name[match(store_aid1, uaid_keep)]
drugnames_lid1 = unq_library_name[match(store_lid1, ulid_keep)]
df_write1 = data.frame(cbind(anchorid=store_aid1, libraryid=store_lid1, anchorname=drugnames_anchor1, libraryname=drugnames_lid1, percent_cellines=store_all_celllines_compared_to_most_sensitive, bestindividual_viab = store_bestindividual_viab1))
df_write1$percent_cellines = as.numeric(as.character(df_write1$percent_cellines))
df_write1_sort = df_write1[order(df_write1$percent_cellines, decreasing=TRUE), ]

drugnames_anchor2 = unq_anchor_name[match(store_aid2, uaid_keep)]
drugnames_lid2 = unq_library_name[match(store_lid2, ulid_keep)]
df_write2 = data.frame(cbind(anchorid=store_aid2, libraryid=store_lid2, anchorname=drugnames_anchor2, libraryname=drugnames_lid2, viab_diff=store_best_all_celllines_minus_most_sensitive, bestindividual_viab = store_bestindividual_viab2))
df_write2$viab_diff = as.numeric(as.character(df_write2$viab_diff))
df_write2_sort = df_write2[order(df_write2$viab_diff, decreasing=FALSE), ]

drugnames_anchor3 = unq_anchor_name[match(store_aid3, uaid_keep)]
drugnames_lid3 = unq_library_name[match(store_lid3, ulid_keep)]
df_write3 = data.frame(cbind(anchorid=store_aid3, libraryid=store_lid3, anchorname=drugnames_anchor3, libraryname=drugnames_lid3, percent_cellines=store3_all_celllines_compared_to_most_sensitive, viab_diff=store3_best_all_celllines_minus_most_sensitive, bestindividual_viab = store_bestindividual_viab3))
df_write3$viab_diff = as.numeric(as.character(df_write3$viab_diff))
df_write3_sort = df_write3[order(df_write3$viab_diff, decreasing=FALSE), ]


drugnames_anchor4 = unq_anchor_name[match(store_aid4, uaid_keep)]
drugnames_lid4 = unq_library_name[match(store_lid4, ulid_keep)]
df_write4 = data.frame(cbind(anchorid=store_aid4, libraryid=store_lid4, anchorname=drugnames_anchor4, libraryname=drugnames_lid4, viab_diff=store_mean_all_celllines_minus_most_sensitive, bestindividual_viab = store_bestindividual_viab4))
df_write4$viab_diff = as.numeric(as.character(df_write4$viab_diff))
df_write4_sort = df_write4[order(df_write4$viab_diff, decreasing=FALSE), ]


drugnames_anchor1_low = unq_anchor_name[match(store_aid1_low, uaid_keep)]
drugnames_lid1_low = unq_library_name[match(store_lid1_low, ulid_keep)]
df_write1_low = data.frame(cbind(anchorid=store_aid1_low, libraryid=store_lid1_low, anchorname=drugnames_anchor1_low, libraryname=drugnames_lid1_low, percent_cellines=store_all_celllines_compared_to_most_sensitive_low, bestindividual_viab = store_bestindividual_viab1_low))
df_write1_low$percent_cellines = as.numeric(as.character(df_write1_low$percent_cellines))
df_write1_low_sort = df_write1_low[order(df_write1_low$percent_cellines, decreasing=TRUE), ]


drugnames_anchor1 = unq_anchor_name[match(store_aid1, uaid_keep)]
drugnames_lid1 = unq_library_name[match(store_lid1, ulid_keep)]
df_write1 = data.frame(cbind(anchorid=store_aid1, libraryid=store_lid1, anchorname=drugnames_anchor1, libraryname=drugnames_lid1, percent_cellines=store_all_celllines_compared_to_most_sensitive, bestindividual_viab = store_bestindividual_viab1))
df_write1$percent_cellines = as.numeric(as.character(df_write1$percent_cellines))
df_write1_sort = df_write1[order(df_write1$percent_cellines, decreasing=TRUE), ]



# write(colnames(df_write1_sort), file="combinations_supersensitive_list_greaterthanbestindividual_atleast_10_percent_cellines_correcteddata_Aug2019.txt", append=FALSE, sep="\t", ncolumns=6)
# write(t(df_write1_sort), file="combinations_supersensitive_list_greaterthanbestindividual_atleast_10_percent_cellines_correcteddata_Aug2019.txt", append=TRUE, sep="\t", ncolumns=6)


supsen = fread("combinations_supersensitive_list_greaterthanbestindividual_atleast_10_percent_cellines_correcteddata_Aug2019.txt", header=TRUE) 
supsen_sel = supsen[,c(1:5)]
# export(supsen_sel, file="selected_combinations_supersensitive_list_greaterthanbestindividual_atleast_10_percent_cellines_correcteddata_Aug2019.txt") ## this file can be used to generate the graph in Fig. 7B


## for figure 7C

supsen = fread("combinations_supersensitive_list_greaterthanbestindividual_atleast_10_percent_cellines_correcteddata_Aug2019.txt", header=TRUE)



syn_per_combo_order_synergy_id = paste(syn_per_combo_order$AnchorID, syn_per_combo_order$LibraryID, sep="_")

supsen = supsen[which(supsen$anchorid != supsen$libraryid),]
supsen_id = paste(supsen$anchorid, supsen$libraryid, sep="_")
other_id = setdiff(syn_per_combo_order_synergy_id, supsen_id)

wlp = wilcox.test(syn_per_combo_order$percentage_synergistic_celllines[match(supsen_id, syn_per_combo_order_synergy_id)], syn_per_combo_order$percentage_synergistic_celllines[match(other_id, syn_per_combo_order_synergy_id)], alternative="greater")$p.value


dat = data.frame(cbind(synval=c(syn_per_combo_order$percentage_synergistic_celllines[match(supsen_id, syn_per_combo_order_synergy_id)], syn_per_combo_order$percentage_synergistic_celllines[match(other_id, syn_per_combo_order_synergy_id)]), type=c(rep("super-sensitizors",length(supsen_id)), rep("others",length(other_id)))))
dat$synval = as.numeric(as.character(dat$synval))

seg_df <- data.frame(x = c(1.1), y = c(21), xend = c(1.9), yend = c(21))

pdf("newfigure_box_synergy_supersensitive_list_correctdata_Aug2019.pdf", width=10,height=10)
p1 <- ggplot(dat, aes(x=type, y=synval, fill=type)) + 
  geom_boxplot(width=0.3,position=position_dodge(0.6)) + theme_bw() + theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30), axis.title.x=element_blank(), axis.title.y=element_text(size=35), plot.title = element_text(size = 30), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none') + ylab("Cell lines with synergy (%)") + coord_cartesian(xlim = NULL, ylim = c(0,27), expand = TRUE, default = FALSE, clip = "on") +  annotate("text", x = c(1.5), y = c(22), label = c(paste("P =",format(wlp,digits=3))), size = 11) +  geom_segment(data=seg_df, mapping=aes(x, y, xend=xend, yend=yend), size=2, color="red", inherit.aes = FALSE ) 
p1
dev.off()


## for figure 7D



topcombos_percent_per_dose_type3_dose4 = matrix(NA,length(uaid_keep), length(ulid_keep))
for (aa in 1:length(uaid_keep)) {
	for (ll in 1:length(ulid_keep)) {
		synviab_latest_zone23_keep_percombo = synviab_latest_zone23_keep[which((uaid_keep[aa]==synviab_latest_zone23_keep$AnchorID) & (ulid_keep[ll]==synviab_latest_zone23_keep$LibraryID)),]
		topcombos_percent_per_dose_type3_dose4[aa,ll] = (sum(((synviab_latest_zone23_keep_percombo$AO_by_OO > 0.75) & (synviab_latest_zone23_keep_percombo$OL_by_OO_D4 > 0.75) & (synviab_latest_zone23_keep_percombo$AL_by_OO_D4 < 0.4)),na.rm=TRUE))/dim(synviab_latest_zone23_keep_percombo)[1]*100

	}	
}

sel_aid = c()
sel_lid = c()
persel = c()
for (aa in 1:length(uaid_keep)) {
	for (ll in 1:length(ulid_keep)) {
		if (topcombos_percent_per_dose_type3_dose4[aa,ll] >= 4) {
			sel_aid = c(sel_aid, uaid_keep[aa])
			sel_lid = c(sel_lid, ulid_keep[ll])
			persel = c(persel, topcombos_percent_per_dose_type3_dose4[aa,ll])
		}
	}
}


df = data.frame(cbind(anchor=sel_aid, library=sel_lid, percentsel=persel))

drugnames_anchor = unq_library_name[match(df$anchor, ulid_keep)]
drugnames_lib = unq_library_name[match(df$library, ulid_keep)]
drugnames_anchor[43:46] = unq_anchor_name[match(308, uaid_keep)]
df_write = data.frame(cbind(anchorid=sel_aid, libraryid=sel_lid, anchorname=drugnames_anchor, libraryname=drugnames_lib, percentsel=persel))
## write(c("anchor_ID", "library_ID", "anchor_name", "library_name", "percent_cell_lines_satisfy_criteria_doseD4"), file="drug_combinations_selected_dose_D4_cyril_hypothesis_higheffectivecombos_viab_lowviabindividually_control_norm3_normalization_with_individualweakdrugs_corrected_Aug2019.txt", ncolumns=5, sep="\t", append=FALSE)
## write(t(df_write), file="drug_combinations_selected_dose_D4_cyril_hypothesis_higheffectivecombos_viab_lowviabindividually_control_norm3_normalization_with_individualweakdrugs_corrected_Aug2019.txt", ncolumns=5, sep="\t", append=TRUE)

dff_cyril = fread("drug_combinations_selected_dose_D4_cyril_hypothesis_higheffectivecombos_viab_lowviabindividually_control_norm3_normalization_with_individualweakdrugs_corrected_Aug2019.txt", header=TRUE)

dff_cyril$anchor_name = new_anchor_name_June2020$AnchorName[match(dff_cyril$anchor_ID, new_anchor_name_June2020$AnchorID)]
dff_cyril$library_name = new_library_name_June2020$Drug.Name[match(dff_cyril$library_ID, new_library_name_June2020$"Cmt.Cpd#")]
## write(names(dff_cyril), file="drug_July2020_combinations_selected_dose_D4_cyril_hypothesis_higheffectivecombos_viab_lowviabindividually_control_norm3_normalization_with_individualweakdrugs_corrected_Aug2019.txt", ncolumns=dim(dff_cyril)[2], sep="\t", append=FALSE)
## write(t(dff_cyril), file="drug_July2020_combinations_selected_dose_D4_cyril_hypothesis_higheffectivecombos_viab_lowviabindividually_control_norm3_normalization_with_individualweakdrugs_corrected_Aug2019.txt", ncolumns=dim(dff_cyril)[2], sep="\t", append=TRUE) ### file to generate figure 7D


## for figure 7E


## PPI network (STRING dataset downloaded on Aug. 8, 2019). 

ppi9606 = fread("9606.protein.info.v11.0.txt", header=TRUE)
ppi9606_links = fread("9606.protein.links.v11.0.txt", header=TRUE)
ppi9606_links_id1 = paste(ppi9606_links$protein1, ppi9606_links$protein2, sep="_")

ppi9606_action = fread("9606.protein.actions.v11.0.txt", header=TRUE)
ppi9606_action = ppi9606_action[which(ppi9606_action$mode == "binding"),]
ppi9606_action_id1 = paste(ppi9606_action$item_id_a, ppi9606_action$item_id_b, sep="_")


cn = 1
ppi_score_median_highsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines >= 10)))
ppi_score_max_highsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines >= 10)))
ppi_score_action_median_highsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines >= 10)))
ppi_score_action_max_highsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines >= 10)))
for (i in which(syn_per_combo_order$percentage_synergistic_celllines >= 10)){
	anchtarget = anchor_drugs_targets_map$relevant_anchorgenes[match(syn_per_combo_order$AnchorID[i], anchor_drugs_targets_map$anchordrugid)]
	libtarget = library_drugs_targets_map$relevant_librarygenes[match(syn_per_combo_order$LibraryID[i], library_drugs_targets_map$librarydrugid)]
	anchtargetunlist = unlist(strsplit(as.character(anchtarget), ","))
	libtargetunlist = unlist(strsplit(as.character(libtarget), ","))
	anchor_protein = ppi9606$protein_external_id[match(anchtargetunlist, ppi9606$preferred_name)]
	lib_protein = ppi9606$protein_external_id[match(libtargetunlist, ppi9606$preferred_name)]
	drugpairid1 = paste(anchor_protein, lib_protein, sep="_")
	drugpairid2 = paste(lib_protein, anchor_protein, sep="_")
	ppi_scorestore = ppi9606_links$combined_score[match(drugpairid1,ppi9606_links_id1)]
	ppi_score_median_highsyn[cn] = median(ppi_scorestore, na.rm=TRUE)
	ppi_score_max_highsyn[cn] = max(ppi_scorestore, na.rm=TRUE)

	drugpairid_comb = union(drugpairid1, drugpairid2)
	ppi_action_scorestore = ppi9606_action$score[match(drugpairid_comb,ppi9606_action_id1)]
	ppi_score_action_median_highsyn[cn] = median(ppi_action_scorestore, na.rm=TRUE)
	ppi_score_action_max_highsyn[cn] = max(ppi_action_scorestore, na.rm=TRUE)

	cn = cn + 1
}

cn = 1
ppi_score_median_lowsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines <= 0)))
ppi_score_max_lowsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines <= 0)))
ppi_score_action_median_lowsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines <= 0)))
ppi_score_action_max_lowsyn = matrix(NA,1,length(which(syn_per_combo_order$percentage_synergistic_celllines <= 0)))
for (i in which(syn_per_combo_order$percentage_synergistic_celllines <= 0)){
	anchtarget = anchor_drugs_targets_map$relevant_anchorgenes[match(syn_per_combo_order$AnchorID[i], anchor_drugs_targets_map$anchordrugid)]
	libtarget = library_drugs_targets_map$relevant_librarygenes[match(syn_per_combo_order$LibraryID[i], library_drugs_targets_map$librarydrugid)]
	anchtargetunlist = unlist(strsplit(as.character(anchtarget), ","))
	libtargetunlist = unlist(strsplit(as.character(libtarget), ","))
	anchor_protein = ppi9606$protein_external_id[match(anchtargetunlist, ppi9606$preferred_name)]
	lib_protein = ppi9606$protein_external_id[match(libtargetunlist, ppi9606$preferred_name)]
	drugpairid1 = paste(anchor_protein, lib_protein, sep="_")
	drugpairid2 = paste(lib_protein, anchor_protein, sep="_")
	ppi_scorestore = ppi9606_links$combined_score[match(drugpairid1,ppi9606_links_id1)]
	ppi_score_median_lowsyn[cn] = median(ppi_scorestore, na.rm=TRUE)
	ppi_score_max_lowsyn[cn] = max(ppi_scorestore, na.rm=TRUE)

	drugpairid_comb = union(drugpairid1, drugpairid2)
	ppi_action_scorestore = ppi9606_action$score[match(drugpairid_comb,ppi9606_action_id1)]
	ppi_score_action_median_lowsyn[cn] = median(ppi_action_scorestore, na.rm=TRUE)
	ppi_score_action_max_lowsyn[cn] = max(ppi_action_scorestore, na.rm=TRUE)

	cn = cn + 1
}

ppi_score_max_highsyn_removeinf = ppi_score_max_highsyn[which(ppi_score_max_highsyn != -Inf)]
ppi_score_max_lowsyn_removeinf = ppi_score_max_lowsyn[which(ppi_score_max_lowsyn != -Inf)]


ppi_score_action_max_highsyn_removeinf = ppi_score_action_max_highsyn[which(ppi_score_action_max_highsyn != -Inf)]
ppi_score_action_max_lowsyn_removeinf = ppi_score_action_max_lowsyn[which(ppi_score_action_max_lowsyn != -Inf)]


wlp = wilcox.test(ppi_score_max_highsyn_removeinf, ppi_score_max_lowsyn_removeinf, alternative="greater")$p.value

df = data.frame(ppiscore=c(ppi_score_max_highsyn_removeinf, ppi_score_max_lowsyn_removeinf), synstat=c(rep("top-synergy",length(ppi_score_max_highsyn_removeinf)), rep("no-synergy",length(ppi_score_max_lowsyn_removeinf))))
seg_df <- data.frame(x = c(1.1), y = c(990), xend = c(1.9), yend = c(990))

pdf("newfigure_boxplot_topbottom_ge10_le0_synergy_PPI_STRING_all_score_max_correctdata_Aug2019.pdf", width=10,height=8) ## figure 7E
ggplot(df, aes(x=synstat, y=(ppiscore), fill=synstat)) + 
   geom_boxplot(width=0.3) + theme_bw() + theme(axis.text.x=element_text(size=35), axis.text.y=element_text(size=35), axis.title.y=element_text(size=35), axis.title.x=element_blank(), plot.title = element_text(size = 40), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none') + ylab("PPI score (max)")  + annotate("text", x = c(1.5), y = c(1020), label = c(paste("P =",format(wlp,digits=3))), size = 11, angle=0) +  geom_segment(data=seg_df, mapping=aes(x, y, xend=xend, yend=yend), size=2, color="red", inherit.aes = FALSE ) 
dev.off() 


## for figure 7F


library("org.Hs.eg.db")

anchordrugid = unique(anchordrugfile$"Compound#")
relevant_anchorgenes = anchordrugfile$"targets_latest"[match(anchordrugid,anchordrugfile$"Compound#")]

anchor_drugs_targets_map = data.frame(cbind(relevant_anchorgenes, anchordrugid))

librarydrugid = unique(librarydrugfile$"CMT_CPD#")
relevant_librarygenes = librarydrugfile$"TARGETGENESLatest"[match(librarydrugid,librarydrugfile$"CMT_CPD#")]


library_drugs_targets_map = data.frame(cbind(relevant_librarygenes, librarydrugid))

all_libgenes = unique(unlist(strsplit(relevant_librarygenes,",")))
all_anchgenes = unique(unlist(strsplit(relevant_anchorgenes,",")))



no_syn = dim(syn_per_combo_order)[1]


avoidcertainindex = c(75)
avoidgenes = c("MT-ATP6", "MT-ATP8")

common_pathway_store = matrix(NA,no_syn,1)
uncommon_pathway_store = matrix(NA,no_syn,1)
commonKegg_store_unq_store = c()
uncommonKegg_store_unq_store = c()
for (i in 1:no_syn){
	anchtarget = anchor_drugs_targets_map$relevant_anchorgenes[match(syn_per_combo_order$AnchorID[i], anchor_drugs_targets_map$anchordrugid)]
	libtarget = library_drugs_targets_map$relevant_librarygenes[match(syn_per_combo_order$LibraryID[i], library_drugs_targets_map$librarydrugid)]
	anchtargetunlist = unlist(strsplit(as.character(anchtarget), ","))
	libtargetunlist = unlist(strsplit(as.character(libtarget), ","))
	if ((sum(anchtargetunlist %in% humangenename$Genename) > 0) & (sum(libtargetunlist %in% humangenename$Genename) > 0)) {
		if (syn_per_combo_order$AnchorID[i] != syn_per_combo_order$LibraryID[i]) {
			commonKegg_store = c()
			uncommonKegg_store = c()
			for (aa in 1:length(anchtargetunlist)) {
				for (ll in 1:length(libtargetunlist)) {
					twoGenes <- c(anchtargetunlist[aa], libtargetunlist[ll])
					if ((sum(twoGenes %in% humangenename$Genename) == 2) & (sum(twoGenes %in% avoidgenes) == 0)) {
						twoGenesEG <- unlist(mget(twoGenes, org.Hs.egSYMBOL2EG))
						twoGenesEGkegg <- mget(twoGenesEG, org.Hs.egPATH)
						commonKegg <- intersect(twoGenesEGkegg[[1]], twoGenesEGkegg[[2]])
						uncommonKegg <- setdiff(union(twoGenesEGkegg[[1]], twoGenesEGkegg[[2]]), intersect(twoGenesEGkegg[[1]], twoGenesEGkegg[[2]]))
						commonKegg_store = c(commonKegg_store, commonKegg)
						uncommonKegg_store = c(uncommonKegg_store, uncommonKegg)
					}
				}
			}
			commonKegg_store_unq = unique(commonKegg_store)
			commonKegg_store_unq = commonKegg_store_unq[!is.na(commonKegg_store_unq)]
			uncommonKegg_store_unq = unique(uncommonKegg_store)
			uncommonKegg_store_unq = uncommonKegg_store_unq[!is.na(uncommonKegg_store_unq)]

			common_pathway_store[i] = length(commonKegg_store_unq)
			uncommon_pathway_store[i] = length(uncommonKegg_store_unq)
			commonKegg_store_unq_store = c(commonKegg_store_unq_store, commonKegg_store_unq)
			uncommonKegg_store_unq_store = c(uncommonKegg_store_unq_store, uncommonKegg_store_unq)
		}
	}
}

combined_pathway_store = common_pathway_store + uncommon_pathway_store

no_of_anchor_targets=matrix(NA,no_syn,1)
no_of_library_targets=matrix(NA,no_syn,1)
for (i in 1:no_syn){
	anchtarget = anchor_drugs_targets_map$relevant_anchorgenes[match(syn_per_combo_order$AnchorID[i], anchor_drugs_targets_map$anchordrugid)]
	libtarget = library_drugs_targets_map$relevant_librarygenes[match(syn_per_combo_order$LibraryID[i], library_drugs_targets_map$librarydrugid)]
	anchtargetunlist = unlist(strsplit(as.character(anchtarget), ","))
	libtargetunlist = unlist(strsplit(as.character(libtarget), ","))
	no_of_anchor_targets[i] = length(anchtargetunlist)
	no_of_library_targets[i] = length(libtargetunlist)
}


## save(common_pathway_store, uncommon_pathway_store, commonKegg_store_unq_store, uncommonKegg_store_unq_store, combined_pathway_store, no_of_anchor_targets, no_of_library_targets, file="commonpathway_kegg_correcteddata_Aug2019.RData")



st_aid_high = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)]
st_lid_high = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)]

st_aid_low = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)]
st_lid_low = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)]


wlp1 = wilcox.test(common_pathway_store[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)],common_pathway_store[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)], alternative="greater")$p.value

wlp2 = wilcox.test(uncommon_pathway_store[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)],uncommon_pathway_store[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)], alternative="greater")$p.value




no_of_combined_targets = no_of_anchor_targets + no_of_library_targets
pos_set = which((syn_per_combo_order$percentage_synergistic_celllines >= 10) & (syn_per_combo_order$AnchorID != syn_per_combo_order$LibraryID))
neg_set = which((syn_per_combo_order$percentage_synergistic_celllines <= 0) & (syn_per_combo_order$AnchorID != syn_per_combo_order$LibraryID))
wlp1 = wilcox.test(((common_pathway_store[pos_set])/(no_of_combined_targets[pos_set])),((common_pathway_store[neg_set])/(no_of_combined_targets[neg_set])), alternative="greater")$p.value


df = data.frame(pathnum=c(((common_pathway_store[pos_set])/(no_of_combined_targets[pos_set])),((common_pathway_store[neg_set])/(no_of_combined_targets[neg_set]))), synstat=c(rep("top-synergy",length(pos_set)), rep("no-synergy",length(neg_set))))
seg_df <- data.frame(x = c(1.1), y = c(0.75), xend = c(1.9), yend = c(0.75))

pdf("newfigure_boxplot_topbottom_ge10_le0_synergy_commonpathwaynumbers_normalized_by_no_of_targets_corrected_Aug2019.pdf", width=10,height=8) ## figure 7F
ggplot(df, aes(x=synstat, y=asinh(pathnum), fill=synstat)) + 
   geom_boxplot(width=0.3) + theme_bw() + theme(axis.text.x=element_text(size=35), axis.text.y=element_text(size=35), axis.title.y=element_text(size=35), axis.title.x=element_blank(), plot.title = element_text(size = 40), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none') + ylab("No. of shared pathways \n normalized by no. of targets") + coord_cartesian(xlim = NULL, ylim = c(0,1), expand = TRUE,
  default = FALSE, clip = "on") + annotate("text", x = c(1.5), y = c(0.8), label = c(paste("P =",format(wlp1,digits=3))), size = 11, angle=0) +  geom_segment(data=seg_df, mapping=aes(x, y, xend=xend, yend=yend), size=2, color="red", inherit.aes = FALSE )
dev.off()



## for figure 7G


load("genenames_tcga_extended.RData")

load("sless1_store_cyrilproject_lungcancertcga.RData")


sl.ess.summary = c()
sl.ess1 = sl.ess1_store
unqgene = unique(sl.ess1_store[,1])
    sl.tot = sl.ess1
    FDR=0.2
    ix=which(p.adjust(sl.tot[,3],"BH")<FDR)
    sl.ess1=sl.tot[ix,]
    flag=0                          
    if (length(ix)>1){                                  # perform FDR of ISLE
        # iy=which(p.adjust(sl.ess1[,5],"BH")<FDR | p.adjust(sl.ess1[,6],"BH")<FDR)   # step I: ## can be & and also... underrepresentation
        iy=which(p.adjust(sl.ess1[,5],"BH")<FDR & p.adjust(sl.ess1[,6],"BH")<FDR)   # step I:  underrepresentation

        if (length(iy)>=1) { 
            sl.ess1=sl.ess1[iy,]                                                    # step II: better prognosis
            if (length(sl.ess1)==24) {
                sl.ess1 = matrix(sl.ess1,1,24)
            }
            coef=sl.ess1[,7];x=sl.ess1[,8];pv=ifelse(coef < 0, x/2, 1-x/2) 
            iz1=which(p.adjust(pv,"BH")<FDR)
            coef=sl.ess1[,15];x=sl.ess1[,16];pv=ifelse(coef < 0, x/2, 1-x/2) 
            iz2=which(p.adjust(pv,"BH")<FDR)
                coef=sl.ess1[iz1,13];x=sl.ess1[iz1,14];pv=ifelse(coef < 0, x/2, 1-x/2) 
            iz1=iz1[which(p.adjust(pv,"BH")<FDR)]
                coef=sl.ess1[iz2,21];x=sl.ess1[iz2,22];pv=ifelse(coef < 0, x/2, 1-x/2) 
            iz2=iz2[which(p.adjust(pv,"BH")<FDR)]
            iz=union(iz1,iz2)
            if (length(iz)>1) {
                sl.ess1=sl.ess1[iz,]        
                ix=which(sl.ess1[,24]< .5)                                          # step III: phylogeny
                if (length(ix)>1)
                    flag=1
            } 
            if (length(iz)==1) {
                sl.ess1=sl.ess1[iz,] 
                ix=which(sl.ess1[24]< .5)                                                 
            } 
            if (length(iz)==0) {
                sl.ess1=sl.ess1[iz,] 
                ix=which(sl.ess1[,24]< .5)                                                 
            } 
        } else {
            ix=c()
        }
    }

    if ((length(ix)>1)) {
        sl.ess1=sl.ess1[ix,]
    }
    if ((length(ix)==1) & (length(sl.ess1)==24)) {
        sl.ess1 = sl.ess1
    }
    if ((length(ix)==1) & (length(sl.ess1)>24)) {    
        sl.ess1=matrix(sl.ess1[ix,],1,24)
    }
    if (length(ix)==0){
        sl.ess1 = c()
    }
  

    if (length(ix)>1) {
        sl.ess2=cbind(genenames_tcga[sl.ess1[,1]],genenames_tcga[sl.ess1[,2]],
                    sl.ess1[,3], apply(cbind(sl.ess1[,5:6]),1,min,na.rm=T),
                    apply(cbind(sl.ess1[,c(8,14,16,22)]),1,min,na.rm=T),
                    sl.ess1[,24]
                    )
        sl.ess.summary = rbind(sl.ess.summary, sl.ess2)
    }
    if (length(ix)==1) {
        sl.ess2=cbind(genenames_tcga[sl.ess1[1]],genenames_tcga[sl.ess1[2]],
                    sl.ess1[3], apply(cbind(sl.ess1[5:6]),2,min,na.rm=T),
                    apply(cbind(sl.ess1[c(8,14,16,22)]),2,min,na.rm=T),
                    sl.ess1[24]
                    )
        sl.ess.summary = rbind(sl.ess.summary, sl.ess2)
    }

sl.ess.summary = data.frame(sl.ess.summary)
colnames(sl.ess.summary) = c("gene1", "gene2", "cell_line_pvalue", "molecular_screen_pvalue", "clinical_survival_pvalue", "phylogenetic_score")
# save(sl.ess.summary, file="sl.ess.summary_lungcancertcga_FDR_0p2_correct_Aug2019.RData")
#write(colnames(sl.ess.summary), file="text.sl.ess.summary_lungcancertcga_FDR_0p2_correct_Aug2019.txt", append=FALSE, sep="\t", ncolumns=6)
#write(t(sl.ess.summary), file="text.sl.ess.summary_lungcancertcga_FDR_0p2_correct_Aug2019.txt", append=TRUE, sep="\t", ncolumns=6)

sllsum = fread("text.sl.ess.summary_lungcancertcga_FDR_0p2_correct_Aug2019.txt", header=TRUE)


anchordrugid = unique(anchordrugfile$"Compound#")
relevant_anchorgenes = anchordrugfile$"targets_latest"[match(uaid_keep,anchordrugfile$"Compound#")]

anchor_drugs_targets_map = data.frame(cbind(relevant_anchorgenes, anchordrugid))


librarydrugid = unique(librarydrugfile$"CMT_CPD#")
relevant_librarygenes = librarydrugfile$"TARGETGENESLatest"[match(ulid_keep,librarydrugfile$"CMT_CPD#")]
combinegens = union(relevant_anchorgenes, relevant_librarygenes)
combinegens = unique(unlist(strsplit(combinegens,",")))

sllsum_select = sllsum[which((sllsum$gene1 %in% combinegens) & (sllsum$gene2 %in% combinegens)),]

#write(colnames(sllsum_select), file="text.sl.ess.summary_July2020_lungcancertcga_FDR_0p2_correct_Aug2019.txt", append=FALSE, sep="\t", ncolumns=6)
#write(t(sllsum_select), file="text.sl.ess.summary_July2020_lungcancertcga_FDR_0p2_correct_Aug2019.txt", append=TRUE, sep="\t", ncolumns=6)

anchordrugid = unique(anchordrugfile$"Compound#")
relevant_anchorgenes = anchordrugfile$"targets_latest"[match(anchordrugid,anchordrugfile$"Compound#")]

anchor_drugs_targets_map = data.frame(cbind(relevant_anchorgenes, anchordrugid))


librarydrugid = unique(librarydrugfile$"CMT_CPD#")
relevant_librarygenes = librarydrugfile$"TARGETGENESLatest"[match(librarydrugid,librarydrugfile$"CMT_CPD#")]

library_drugs_targets_map = data.frame(cbind(relevant_librarygenes, librarydrugid))


sl.ess.id = paste(sl.ess.summary[,1], sl.ess.summary[,2], sep="_")

st_aid_high = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)]
st_lid_high = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)]

st_aid_low = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)]
st_lid_low = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)]

st_id = c(paste(st_aid_high, st_lid_high, sep="_"), paste(st_aid_low, st_lid_low, sep="_")) 
class_id = as.factor(c(rep(1,length(st_aid_high)), rep(-1,length(st_aid_low))))

sum_target_pairs_sl = matrix(NA,length(anchordrugid),length(librarydrugid))
anchordrugid_SL = c()
librarydrugid_SL = c()
pairs_SL = c()
cn = 0
for (anch in 1:length(anchordrugid))
{
    print(anch)
    anchorgene = unlist(strsplit(relevant_anchorgenes[anch],","))
    for (lib in 1:length(librarydrugid))
    {
        libgene = unlist(strsplit(relevant_librarygenes[lib],","))
        sl1=cbind(rep(anchorgene,each=length(libgene)),rep(libgene,length(anchorgene)))
        sl_drug = sl1[sl1[,1]!=sl1[,2],]
            
        if ((length(anchorgene)<4) & (length(libgene)<4)) { 

        	cn = cn + 1 

            if (length(sl_drug)==2){
                sl_drug = matrix(sl_drug,1,2)
            }
            
            if ((anchordrugid[anch] %in% uaid_keep) & (librarydrugid[lib] %in% ulid_keep)) {

                sl_drug_id = paste(sl_drug[,1],sl_drug[,2],sep="_")
                sum_target_pairs_sl[anch,lib] = sum(sl_drug_id %in% sl.ess.id)
                if( sum(sl_drug_id %in% sl.ess.id) > 0) {
                    anchordrugid_SL = c(anchordrugid_SL, anchordrugid[anch])
                    librarydrugid_SL = c(librarydrugid_SL, librarydrugid[lib])
                    pairs_SL = c(pairs_SL, sum(sl_drug_id %in% sl.ess.id))
                }
            }
        }
    }
}


dt=data.frame(cbind(anchor_ID=anchordrugid_SL, library_ID=librarydrugid_SL, anchor=unq_anchor_name[match(anchordrugid_SL, uaid_keep)], library=unq_library_name[match(librarydrugid_SL, ulid_keep)],  SLpairs = pairs_SL))
#write(colnames(dt), file="drugcombinations_lessthan4targets_with_atleast_one_SL_pair_correcteddata_Aug2019.txt", ncolumns=5, sep="\t", append=FALSE)
#write(t(dt), file="drugcombinations_lessthan4targets_with_atleast_one_SL_pair_correcteddata_Aug2019.txt", ncolumns=5, sep="\t", append=TRUE)

highthresh = quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.75)
lowthresh = quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.25)
cn_small_target = 0
cn_small_target_index = c()
cn_small_target_highsyn = 0
cn_small_target_highsyn_index = c()
cn_small_target_lowsyn = 0
cn_small_target_lowsyn_index = 0
for (i in 1:dim(syn_per_combo_order)[1]){
    anchorgene = unlist(strsplit(relevant_anchorgenes[match(syn_per_combo_order$AnchorID[i], anchordrugid)],","))
    libgene = unlist(strsplit(relevant_librarygenes[match(syn_per_combo_order$LibraryID[i], librarydrugid)],","))
    if ((length(anchorgene)<4) & (length(libgene)<4)) { 
    	cn_small_target = cn_small_target + 1 
    	cn_small_target_index = c(cn_small_target_index, i)
    }
    if ((length(anchorgene)<4) & (length(libgene)<4) & (syn_per_combo_order$percentage_synergistic_celllines[i] >= highthresh)) { 
    	cn_small_target_highsyn = cn_small_target_highsyn + 1 
    	cn_small_target_highsyn_index = c(cn_small_target_highsyn_index, i)
    }    
    if ((length(anchorgene)<4) & (length(libgene)<4) & (syn_per_combo_order$percentage_synergistic_celllines[i] <= lowthresh)) { 
    	cn_small_target_lowsyn = cn_small_target_lowsyn + 1 
    	cn_small_target_lowsyn_index = c(cn_small_target_lowsyn_index, i)
    }  
}


combowithSL = fread("drugcombinations_lessthan4targets_with_atleast_one_SL_pair_correcteddata_Aug2019.txt", header=TRUE)
st_aid_high = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines >= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.75))]
st_lid_high = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines >= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.75))]

st_aid_low = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines <= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.25))]
st_lid_low = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines <= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.25))]

combid_all = paste(syn_per_combo_order$AnchorID, syn_per_combo_order$LibraryID, sep="_")

combid_with_sl = paste(combowithSL$anchor_ID, combowithSL$library_ID, sep="_")
combid_with_highsynergy = paste(st_aid_high, st_lid_high, sep="_")
combid_with_lowsynergy = paste(st_aid_low, st_lid_low, sep="_")

length(intersect(combid_with_sl,combid_all[cn_small_target_highsyn_index]))

int_combos_syn_sl = intersect(combid_with_sl,combid_all[cn_small_target_highsyn_index])
notint_combos_syn_sl = setdiff(combid_all[cn_small_target_highsyn_index],combid_with_sl)


top_syn_per_combo_order_SL = syn_per_combo_order[match(int_combos_syn_sl, paste(syn_per_combo_order$AnchorID, syn_per_combo_order$LibraryID, sep="_")),]
top_syn_per_combo_order_SL$SLpairs = combowithSL$SLpairs[match(int_combos_syn_sl, combid_with_sl)]

# write(colnames(top_syn_per_combo_order_SL), file="top_synergisticdrugcombinations_lessthan4targets_with_atleast_one_SL_pair_corrected_Aug2019.txt", ncolumns=6, sep="\t", append=FALSE)
# write(t(top_syn_per_combo_order_SL), file="top_synergisticdrugcombinations_lessthan4targets_with_atleast_one_SL_pair_corrected_Aug2019.txt", ncolumns=6, sep="\t", append=TRUE)



topsyn_SL = fread("top_synergisticdrugcombinations_lessthan4targets_with_atleast_one_SL_pair_corrected_Aug2019.txt", header=TRUE)    
topsyn_SL$Anchor_Drug = new_anchor_name_June2020$AnchorName[match(topsyn_SL$AnchorID, new_anchor_name_June2020$AnchorID)]
topsyn_SL$Library_Drug = new_library_name_June2020$Drug.Name[match(topsyn_SL$LibraryID, new_library_name_June2020$"Cmt.Cpd#")]
## write(c("Anchor_ID", "Library_ID", "Anchor_Drug_Name", "Library_Drug_Name", "percentage_synergistic_celllines", "no_of_SL_pairs"), file="top_July2020_synergisticdrugcombinations_lessthan4targets_with_atleast_one_SL_pair_corrected_Aug2019.txt", ncolumns=dim(topsyn_SL)[2], sep="\t", append=FALSE)
## write(t(topsyn_SL), file="top_July2020_synergisticdrugcombinations_lessthan4targets_with_atleast_one_SL_pair_corrected_Aug2019.txt", ncolumns=dim(topsyn_SL)[2], sep="\t", append=TRUE) ## this file can be used to create Figure 7G
   

## for figure 7H

exp_SL = fread("experimental_crispr_SLpooled_data_Jan2019.txt", header=FALSE)


anchordrugid = unique(anchordrugfile$"Compound#")
relevant_anchorgenes = anchordrugfile$"targets_latest"[match(anchordrugid,anchordrugfile$"Compound#")]

anchor_drugs_targets_map = data.frame(cbind(relevant_anchorgenes, anchordrugid))


librarydrugid = unique(librarydrugfile$"CMT_CPD#")
relevant_librarygenes = librarydrugfile$"TARGETGENESLatest"[match(librarydrugid,librarydrugfile$"CMT_CPD#")]

library_drugs_targets_map = data.frame(cbind(relevant_librarygenes, librarydrugid))


exp.sl.ess.id = union(paste(exp_SL$V1, exp_SL$V2, sep="_"),paste(exp_SL$V2, exp_SL$V1, sep="_"))

st_aid_high = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)]
st_lid_high = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines >= 10)]

st_aid_low = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)]
st_lid_low = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines <= 0)]

st_id = c(paste(st_aid_high, st_lid_high, sep="_"), paste(st_aid_low, st_lid_low, sep="_")) 
class_id = as.factor(c(rep(1,length(st_aid_high)), rep(-1,length(st_aid_low))))

sum_target_pairs_sl = matrix(NA,length(anchordrugid),length(librarydrugid))
anchordrugid_SL = c()
librarydrugid_SL = c()
pairs_SL = c()
cn = 0
for (anch in 1:length(anchordrugid))
{
    print(anch)
    anchorgene = unlist(strsplit(relevant_anchorgenes[anch],","))
    for (lib in 1:length(librarydrugid))
    {
        libgene = unlist(strsplit(relevant_librarygenes[lib],","))
        sl1=cbind(rep(anchorgene,each=length(libgene)),rep(libgene,length(anchorgene)))
        sl_drug = sl1[sl1[,1]!=sl1[,2],]
            
        if ((length(anchorgene)<4) & (length(libgene)<4)) { 

        	cn = cn + 1 

            if (length(sl_drug)==2){
                sl_drug = matrix(sl_drug,1,2)
            }
            
            if ((anchordrugid[anch] %in% uaid_keep) & (librarydrugid[lib] %in% ulid_keep)) {

                sl_drug_id = paste(sl_drug[,1],sl_drug[,2],sep="_")
                sum_target_pairs_sl[anch,lib] = sum(sl_drug_id %in% exp.sl.ess.id)
                if( sum(sl_drug_id %in% exp.sl.ess.id) > 0) {
                    anchordrugid_SL = c(anchordrugid_SL, anchordrugid[anch])
                    librarydrugid_SL = c(librarydrugid_SL, librarydrugid[lib])
                    pairs_SL = c(pairs_SL, sum(sl_drug_id %in% exp.sl.ess.id))
                }
            }
        }
    }
}


dt=data.frame(cbind(anchor_ID=anchordrugid_SL, library_ID=librarydrugid_SL, anchor=unq_anchor_name[match(anchordrugid_SL, uaid_keep)], library=unq_library_name[match(librarydrugid_SL, ulid_keep)],  SLpairs = pairs_SL))
#write(colnames(dt), file="drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", ncolumns=5, sep="\t", append=FALSE)
#write(t(dt), file="drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", ncolumns=5, sep="\t", append=TRUE)



drugcom_comp = fread("drugcombinations_lessthan4targets_with_atleast_one_SL_pair_correcteddata_Aug2019.txt", header=TRUE)
drugcom_exp = fread("drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", header=TRUE)

drugcommonid = intersect(paste(drugcom_comp$anchor_ID, drugcom_comp$library_ID, sep="_"), paste(drugcom_exp$anchor_ID, drugcom_exp$library_ID, sep="_"))
drugcommonid_mat = matrix(unlist(strsplit(drugcommonid, "_")), ncol=2, byrow=TRUE)

drugcommon_comp_expSL = cbind(unq_anchor_name[match(drugcommonid_mat[,1], uaid_keep)], unq_library_name[match(drugcommonid_mat[,2], ulid_keep)])

#write(c("anchordrug", "librarydrug"), file="drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_AND_COMPUTATIONAL_SL_pair_correcteddata_Aug2019.txt", ncolumns=2, sep="\t", append=FALSE)
#write(t(drugcommon_comp_expSL), file="drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_AND_COMPUTATIONAL_SL_pair_correcteddata_Aug2019.txt", ncolumns=2, sep="\t", append=TRUE)

drugcom_exp = fread("drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", header=TRUE)
drugcom_exp$anchor = new_anchor_name_June2020$AnchorName[match(drugcom_exp$anchor_ID, new_anchor_name_June2020$AnchorID)]
drugcom_exp$library = new_library_name_June2020$Drug.Name[match(drugcom_exp$library_ID, new_library_name_June2020$"Cmt.Cpd#")]
## write(c("Anchor_ID", "Library_ID", "Anchor_Drug_Name", "Library_Drug_Name", "no_of_SL_pairs"), file="drugcombinations_July2020_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", ncolumns=dim(drugcom_exp)[2], sep="\t", append=FALSE)
## write(t(drugcom_exp), file="drugcombinations_July2020_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", ncolumns=dim(drugcom_exp)[2], sep="\t", append=TRUE)
   

combowithSL = fread("drugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_correcteddata_Aug2019.txt", header=TRUE)
st_aid_high = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines >= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.75))]
st_lid_high = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines >= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.75))]

st_aid_low = syn_per_combo_order$AnchorID[which(syn_per_combo_order$percentage_synergistic_celllines <= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.25))]
st_lid_low = syn_per_combo_order$LibraryID[which(syn_per_combo_order$percentage_synergistic_celllines <= quantile(syn_per_combo_order$percentage_synergistic_celllines, prob=0.25))]

combid_all = paste(syn_per_combo_order$AnchorID, syn_per_combo_order$LibraryID, sep="_")

combid_with_sl = paste(combowithSL$anchor_ID, combowithSL$library_ID, sep="_")
combid_with_highsynergy = paste(st_aid_high, st_lid_high, sep="_")
combid_with_lowsynergy = paste(st_aid_low, st_lid_low, sep="_")


int_combos_syn_sl = intersect(combid_with_sl,combid_all[cn_small_target_highsyn_index])
notint_combos_syn_sl = setdiff(combid_all[cn_small_target_highsyn_index],combid_with_sl)


top_syn_per_combo_order_SL = syn_per_combo_order[match(int_combos_syn_sl, paste(syn_per_combo_order$AnchorID, syn_per_combo_order$LibraryID, sep="_")),]
top_syn_per_combo_order_SL$SLpairs = combowithSL$SLpairs[match(int_combos_syn_sl, combid_with_sl)]

# write(colnames(top_syn_per_combo_order_SL), file="top_synergisticdrugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_corrected_Aug2019.txt", ncolumns=6, sep="\t", append=FALSE)
# write(t(top_syn_per_combo_order_SL), file="top_synergisticdrugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_corrected_Aug2019.txt", ncolumns=6, sep="\t", append=TRUE)


    
topsyn_SL = fread("top_synergisticdrugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_corrected_Aug2019.txt", header=TRUE)
topsyn_SL$Anchor_Drug = new_anchor_name_June2020$AnchorName[match(topsyn_SL$AnchorID, new_anchor_name_June2020$AnchorID)]
topsyn_SL$Library_Drug = new_library_name_June2020$Drug.Name[match(topsyn_SL$LibraryID, new_library_name_June2020$"Cmt.Cpd#")]
## write(c("Anchor_ID", "Library_ID", "Anchor_Drug_Name", "Library_Drug_Name", "percentage_synergistic_celllines", "no_of_SL_pairs"), file="top_July2020_synergisticdrugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_corrected_Aug2019.txt", ncolumns=dim(topsyn_SL)[2], sep="\t", append=FALSE)
## write(t(topsyn_SL), file="top_July2020_synergisticdrugcombinations_lessthan4targets_with_atleast_one_EXPERIMENTAL_SL_pair_corrected_Aug2019.txt", ncolumns=dim(topsyn_SL)[2], sep="\t", append=TRUE) ## this file can be used to generate Figure 7H
   




