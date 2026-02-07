
###################### Subglacial precipitates figure making code

# Bianca De Sanctis
# bddesanctis@gmail.com
# Sep 2025
# This R code makes all of the main plots and most of the supplement plots for the paper,
#   and takes as input only text and excel files available on the github. 

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork) 
library(scales)
library(RColorBrewer)
library(colorspace)
library(cowplot)
library(vegan)
library(grid)
library(reshape2)
library(grid)
library(gridExtra)
library(stats)
library(ggrepel)
library(data.table)
library(RColorBrewer)
library(data.table)

##################### Get all the input ######################

# set this to the metadata sheet in the supplement of the paper
metadata = read_excel("metadata.xlsx") 

# this is available on the github here https://github.com/bdesanctis/subglacial-precipitates/blob/main/data/all_min50reads.tsv
new = read.table("all_min50reads.tsv",header = TRUE, sep="\t")

# this is available on the github here https://github.com/bdesanctis/subglacial-precipitates/blob/main/data/alllinecounts.txt
decontamination_df = read.table("alllinecounts.txt",header=TRUE)

# https://github.com/bdesanctis/subglacial-precipitates/blob/main/data/oxymetag_results.csv
oxymetag = read.table("oxymetag_results.csv",sep=",",header=TRUE)

# https://github.com/bdesanctis/subglacial-precipitates/blob/main/data/SLM-SLW-BFB-OxTol.xlsx
existing_oxygen_tolerances = read_excel("SLM-SLW-BFB-OxTol.xlsx")

# https://github.com/bdesanctis/subglacial-precipitates/blob/main/data/Literature_FeMn_data_summary.xlsx
geochem_file = "Literature_FeMn_data_summary.xlsx"
all_geochem_sheets <- lapply(excel_sheets(geochem_file), function(s) {
	read_excel(geochem_file, sheet = s)
})
names(all_geochem_sheets) <- excel_sheets(geochem_file)

# https://github.com/bdesanctis/subglacial-precipitates/blob/main/data/BR_interpolated.txt
barnaby_grab = read.table("from-terry/BR_interpolated.txt",header=TRUE,sep="\t")

###################### Set up  ######################

# process metadata; assign colours

small_metadata <- na.omit(data.frame(
	Sample = metadata$`Paper code`,
	Age = as.numeric(metadata$`Approx age (ka)`)
))
small_metadata$color <- rep("black", nrow(small_metadata))
positive_control_indices <- which(small_metadata$Sample == "PC1")
small_metadata$color[positive_control_indices] <- "cadetblue3"
sunset_palette <- hcl.colors(5, "Sunset")  
sunset_palette[5] <- "#E8C766"  # Replace the 5th color
non_control_indices <- which(small_metadata$Sample != "PC1")
age_bins <- cut(small_metadata$Age[non_control_indices], 
								breaks = c(0, 25, 99, 199, 300, 600), 
								labels = c("16-25", "26-99", "100-199", "200-300", "300-600"),
								include.lowest = TRUE, right = TRUE)
small_metadata$color[non_control_indices] <- sunset_palette[as.numeric(age_bins)]

# extract node tax level
new$taxlevel = unname(sapply( sapply(new$taxpath, function(x) strsplit(x,split=";")[[1]][1]) , function(x) strsplit(x,split=":")[[1]][3]))
# extract phyla names
new$phylum = gsub("p__","",unname(sapply( sapply(new$taxpath, function(x) rev(strsplit(x,split=";")[[1]])[3]) , function(x) strsplit(x,split=":")[[1]][2])))
# extract read cols
readcols = grep("_TotalReads",colnames(new))
sample_names = unname(sapply(colnames(new[,readcols]), function(x) gsub("_TotalReads","",x)))
neg_control_cols = grep("NC",colnames(new))

total_assigned_reads_per_sample = colSums(new[grep("p__",new$TaxName),grep("_TotalReads",colnames(new))])
	

###################### Fig 1: Maps ######################

# blood falls: Latitude: -77° 42' 59.99" S  , Longitude: 162° 15' 60.00" E
# whillans Coordinates	84°15′S 153°30′W
# mercer Coordinates	84.661°S 149.677°W
# pull out the lat long of everything with a numeric lat long and with a filename
md_samples = data.frame(metadata[which(!is.na(metadata$Filename)),
																 which(colnames(metadata) %in% c("Approx age (ka)","Filename","Lat","Lon","GeoCluster","Paper code"))])
# set up data 
controls = c(""NC1","NC2","NC3","NC4","NC5","NC6","NC7")
samples_to_map = md_samples[! md_samples$Filename %in% controls,]
colnames(samples_to_map)[1] = "Age"
all_ages = as.numeric(samples_to_map$Age[!is.na(samples_to_map$Age)])

# Add the three special points
special_points <- data.frame(
	Lat = c(-77.716664, -84.25, -84.661),
	Lon = c(162.266667, -153.5, -149.677),
	Label = c("BFB", "SLW", "SLM")
)

# define the map
world <- map_data("world")
samples_to_map <- samples_to_map %>%
	mutate(PaperGroup = str_extract(Paper.code, "^[A-Za-z]+")) 
labels_df <- samples_to_map %>%
	group_by(PaperGroup) %>%
	slice(1) %>%
	ungroup()
antmap_grey <- ggplot(world, aes(x = long, y = lat, group = group)) + 
	geom_polygon(fill = "gray92", color = "#4D4D4D") +
	geom_point(data = samples_to_map, aes(x = Lon, y = Lat),
						 size = 3, shape = 21, fill = "grey36", color = "black", stroke = 0.5, inherit.aes = FALSE) +
	geom_text(data = labels_df, aes(x = Lon, y = Lat, label = PaperGroup),
						inherit.aes = FALSE, size = 3, vjust = -1) +  # adjust vjust to shift labels upward
	geom_point(data = special_points, aes(x = Lon, y = Lat),
						 size = 3, shape = 21, fill = "pink", color = "black", stroke = 0.5, inherit.aes = FALSE) +
	geom_text(data = special_points, aes(x = Lon, y = Lat, label = Label),
						inherit.aes = FALSE, size = 3, vjust = -1) +
	coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90, -60),
						xlim = c(-50, 180)) +
	theme_minimal() +
	theme(
		panel.background = element_rect(fill = "lightblue", color = NA),
		panel.border = element_rect(colour = "black", fill = NA),
		panel.grid = element_blank(),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
		axis.title = element_blank()
	) +
	labs(title = "Antarctica")
baffin_map_grey <- ggplot(world, aes(x = long, y = lat, group = group)) + 
	geom_polygon(fill = "gray92", color = "#4D4D4D") + # Light grey land with dark grey outlines
	geom_point(data = samples_to_map, aes(x = Lon, y = Lat), 
						 size = 3, shape = 21, fill = "grey36", color = "black", stroke = 0.5, inherit.aes = FALSE) +
	geom_text(data = labels_df, aes(x = Lon, y = Lat, label = PaperGroup),
						inherit.aes = FALSE, size = 3, vjust = -1) +  # adjust vjust to shift labels upward
	coord_map("ortho", orientation = c(75, -70, 0), 
						ylim = c(60, 80), xlim = c(-92, -55)) + # Zoom in to the island region
	theme_minimal() +
	theme(
		panel.background = element_rect(fill = "lightblue", color = NA), # Soft blue for water
		panel.border = element_rect(colour = "black", fill = NA), # Add a subtle black border
		panel.grid = element_blank(), # No gridlines
		axis.text = element_blank(), # Remove axis text
		axis.ticks = element_blank(), # Remove axis ticks
		axis.title = element_blank(),
		legend.position = "none"
	) +
	labs(title = "Baffin Island") # Add title and legend label
# COMBINE THE MAPS
combined_map_grey <- baffin_map_grey + antmap_grey + 
	plot_layout(ncol = 2) # Arrange plots in 2 columns
# DISPLAY THE COMBINED MAP
print(combined_map_grey)







###################### Supplementary figure: What is the effect of removing control taxa? ######################

taxa_data <- decontamination_df %>%
	select(Sample, TaxaWithReadsOver50KeepControlTax, TaxaWithReadsOver50) %>%
	rename(`Before filtering` = TaxaWithReadsOver50KeepControlTax,
				 `After filtering` = TaxaWithReadsOver50) %>%
	pivot_longer(cols = c(`Before filtering`, `After filtering`),
							 names_to = "Status",
							 values_to = "Count") %>%
	mutate(Status = factor(Status, levels = c("Before filtering", "After filtering")))

lca_data <- decontamination_df %>%
	select(Sample, SmallLcaLineCountKeepControlTax, SmallLcaLineCount) %>%
	rename(`Before filtering` = SmallLcaLineCountKeepControlTax,
				 `After filtering` = SmallLcaLineCount) %>%
	pivot_longer(cols = c(`Before filtering`, `After filtering`),
							 names_to = "Status",
							 values_to = "Count") %>%
	mutate(Status = factor(Status, levels = c("Before filtering", "After filtering")))

p1 <- ggplot(taxa_data, aes(x = Sample, y = Count, fill = Status)) +
	geom_bar(stat = "identity", position = "dodge", width = 0.7) +
	scale_fill_manual(values = c("Before filtering" = "salmon", "After filtering" = "#3498DB")) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
				legend.position = "none",
				panel.grid.major.x = element_blank()) +
	labs(title = "Number of taxa with at least 50 assigned reads",
			 x = "Sample",
			 y = "Number of taxa")

p2 <- ggplot(lca_data, aes(x = Sample, y = Count, fill = Status)) +
	geom_bar(stat = "identity", position = "dodge", width = 0.7) +
	scale_fill_manual(values = c("Before filtering" = "salmon", "After filtering" = "#3498DB")) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
				legend.position = "none",
				panel.grid.major.x = element_blank()) +
	labs(title = "Number of assigned reads",
			 x = "Sample",
			 y = "Number of reads")

p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "right", legend.title = element_blank())

###################### Fig 2: Most abundant genera x median damage #########################
# figure showing most abundant genera damage and abundances
# exclude your positive control from the damage mean

tax_level = "genus"
pad_tax_path = TRUE
readcutoff = 40000
minreads= 200
minreads_pc1 = 50
min_samples = 3

relrows = grep(tax_level, new$taxlevel)
newp = new[relrows,-neg_control_cols]
smallss = newp[newp$TotalReads > readcutoff,]

long_data <- smallss[,-3] %>%
	pivot_longer(
		cols = any_of(paste0(rep(sample_names, each = 2), c("_TotalReads", "_Damage.1"))),
		names_to = c("sample", ".value"),
		names_pattern = "^(.*?)_(TotalReads|Damage.1)$"
	) %>%
	rename(reads = TotalReads, damage = Damage.1) %>%
	filter((sample == "PC1" & reads > minreads_pc1) | (sample != "PC1" & reads > minreads)) %>%
	group_by(TaxName) %>%
	filter(n_distinct(sample) >= min_samples) %>%
	ungroup()

# order the x-axis by mean damage (excluding PC1 and NCs)
long_data <- long_data %>%
	group_by(TaxName) %>%
	mutate(mean_damage = median(damage[sample != "PC1" & !grepl("^NC", sample)], na.rm = TRUE)) %>%
	ungroup() %>%
	mutate(TaxName = reorder(TaxName, mean_damage))

long_data <- long_data %>%
	group_by(TaxName) %>%
	mutate(has_PC1 = any(sample == "PC1")) %>%
	ungroup()

# here the box plots are [.25,.75] quartiles and show the median by default, and the black line shows the mean 
# Damage plot - remove the labels
damageplot = ggplot() + 
	geom_vline(
		xintercept = seq_along(levels(long_data$TaxName)),
		color = "gray85",
		linewidth = 0.3,
		alpha = 0.5
	) +
	geom_boxplot(
		data = filter(long_data, sample != "PC1"),
		aes(x = TaxName, y = damage, fill = has_PC1),
		alpha = 0.5,
		outlier.shape = NA, 
		color = "gray70", 
		linewidth = 0.4,
		width = 0.7
	) +
	geom_point(
		data = filter(long_data, sample != "PC1"),
		aes(x = TaxName, y = damage),
		color = "steelblue", 
		alpha = 0.5,  
		size = 2,
		position = position_jitter(width = 0.25, seed = 42)
	) +
	geom_point(
		data = filter(long_data, sample == "PC1"),
		aes(x = TaxName, y = damage),
		color = "black",  
		size = 3.5,
		alpha = 0.8,
		shape = 21, 
		fill = "firebrick2",
		stroke = .5
	) +
	geom_line(
		data = long_data, 
		aes(x = TaxName, y = mean_damage, group = 1), 
		color = "#2C3E50", 
		linetype = "solid", 
		alpha = 0.8, 
		linewidth = .5
	) +
	scale_fill_manual(
		values = c("TRUE" = "#F8F9FA", "FALSE" = "#F8F9FA")
	) +
	scale_x_discrete(
		labels = NULL  # Remove labels from top plot
	) +
	theme_minimal() +
	theme(
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.y = element_text(size = 10, color = "gray30"),
		axis.title = element_text(size = 12, color = "gray20", face = "bold"),
		axis.title.x = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
		panel.background = element_rect(fill = "white", color = NA),
		plot.background = element_rect(fill = "white", color = NA),
		legend.position = "none"
	) +
	labs(
		x = NULL,
		y = "DNA Damage",
		title = "Most abundant genera, ordered by increasing median damage")


# damageplot

keep_taxa <- levels(long_data$TaxName)
original_data <- new %>%
	filter(TaxName %in% keep_taxa)
sample_order <- c("PC1","AIS1","EM1","EM2","EM3","EM4","MBL1","MBL2","MV1","MV2","RM1","RM2",
									"BV1","BV2","BV3","MA1","MA2","PM1","PM2",
									"LG1","LG2","LG3","LG4","LG5","LG6","LG7")
read_cols <- paste0(sample_order, "_TotalReads")
heat_data <- original_data %>%
	select(TaxName, all_of(read_cols)) %>%
	pivot_longer(
		cols = ends_with("_TotalReads"),
		names_to = "sample",
		values_to = "reads"
	) %>%
	mutate(sample = str_remove(sample, "_TotalReads"))
heat_data <- heat_data %>%
	mutate(
		reads = ifelse(reads < 50, NA, reads),
		reads_log10 = log10(reads + 1)
	)
heat_data <- heat_data %>%
	mutate(
		sample = factor(sample, levels = rev(sample_order)),
		TaxName = factor(TaxName, levels = keep_taxa)
	)
# create custom color mapping for PC1 (reds) vs others (blues)
blues <- colorRampPalette(brewer.pal(9, "Blues"))(100)
reds <- colorRampPalette(c("white", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"))(100)
# create separate data for blue and red tiles to get separate legends
heat_data_blue <- heat_data %>%
	filter(sample != "PC1") %>%
	mutate(reads_log10_blue = reads_log10)
heat_data_red <- heat_data %>%
	filter(sample == "PC1") %>%
	mutate(reads_log10_red = reads_log10)

p1 <- ggplot() +
	geom_tile(data = heat_data_blue, 
						aes(x = TaxName, y = sample, fill = reads_log10_blue), 
						color = "white", linewidth = 0.2) +
	geom_tile(data = heat_data_red, 
						aes(x = TaxName, y = sample, fill = reads_log10_red), 
						color = "white", linewidth = 0.2) +
	scale_fill_distiller(
		name = "Reads", 
		palette = "Blues", 
		direction = 1, 
		na.value = "white", 
		aesthetics = "fill",
		breaks = log10(c(50,100, 1000, 10000, 100000, 500000)),
		labels = c("50","100", "1,000", "10,000", "100,000", "500,000")
	) +
	ggnewscale::new_scale_fill() +
	geom_tile(data = heat_data_red, 
						aes(x = TaxName, y = sample, fill = reads_log10_red), 
						color = "white", linewidth = 0.2) +
	scale_fill_distiller(
		name = "Reads", 
		palette = "Reds", 
		direction = 1, 
		na.value = "white",
		breaks = log10(c(50,100, 1000, 5000, 30000)),
		labels = c("50","100", "1,000", "5,000", "30,000")
	) +
	scale_x_discrete(
		limits = keep_taxa,
		labels = function(x) {
			taxname <- gsub("^g__", "", x)
			phylum <- long_data$phylum[match(x, long_data$TaxName)]
			paste(phylum, taxname, sep = " : ")
		}
	) +
	scale_y_discrete(limits = rev(sample_order)) +
	theme_bw(base_size = 11) +
	theme(
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
		axis.text.y = element_text(size = 8),
		axis.title = element_text(size = 12, face = "bold"),
		panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
		panel.grid.minor = element_blank(),
		plot.background = element_rect(fill = "white", color = NA)
	) +
	labs(x = "Genus", y = "Sample", title = "Per-sample read counts")

damageplot / p1



###################### Supplementary figure: Most abundant (other tax levels) x median damage ###################### 

# turn it into a function
plot_increasing_mean_damage <- function(tax_level, pad_tax_path, readcutoff, minreads, minreads_pc1, min_samples){
	relrows = grep(tax_level, new$taxlevel)
	newp = new[relrows,-neg_control_cols]
	smallss = newp[newp$TotalReads > readcutoff,]
	long_data <- smallss[,-3] %>%
		pivot_longer(
			cols = any_of(paste0(rep(sample_names, each = 2), c("_TotalReads", "_Damage.1"))),
			names_to = c("sample", ".value"),
			names_pattern = "^(.*?)_(TotalReads|Damage.1)$"
		) %>%
		rename(reads = TotalReads, damage = Damage.1) %>%
		filter((sample == "PC1" & reads > minreads_pc1) | (sample != "PC1" & reads > minreads)) %>%
		group_by(TaxName) %>%
		filter(n_distinct(sample) >= min_samples) %>%
		ungroup()
	
	# order the x-axis by mean damage (excluding PC1 and NCs)
	long_data <- long_data %>%
		group_by(TaxName) %>%
		mutate(mean_damage = median(damage[sample != "PC1" & !grepl("^NC", sample)], na.rm = TRUE)) %>%
		ungroup() %>%
		mutate(TaxName = reorder(TaxName, mean_damage))
	
	long_data <- long_data %>%
		group_by(TaxName) %>%
		mutate(has_PC1 = any(sample == "PC1")) %>%
		ungroup()
	
	# here the box plots are [.25,.75] quartiles and show the median by default, and the black line shows the mean 
	ggplot() + 
		geom_vline(
			xintercept = seq_along(levels(long_data$TaxName)),
			color = "gray85",
			linewidth = 0.3,
			alpha = 0.5
		) +
		geom_boxplot(
			data = filter(long_data, sample != "PC1"),
			aes(x = TaxName, y = damage, fill = has_PC1),
			alpha = 0.5,
			outlier.shape = NA, 
			color = "gray70", 
			linewidth = 0.4,
			width = 0.7
		) +
		geom_point(
			data = filter(long_data, sample != "PC1"),
			aes(x = TaxName, y = damage),
			color = "steelblue", 
			alpha = 0.5,  
			size = 2,
			position = position_jitter(width = 0.25, seed = 42)
		) +
		geom_point(
			data = filter(long_data, sample == "PC1"),
			aes(x = TaxName, y = damage),
			color = "black",  
			size = 3.5,
			alpha = 0.8,
			shape = 21, 
			fill = "firebrick2",
			stroke = .5
		) +
		geom_line(
			data = long_data, 
			aes(x = TaxName, y = mean_damage, group = 1), 
			color = "#2C3E50", 
			linetype = "solid", 
			alpha = 0.8, 
			linewidth = .5
		) +
		scale_fill_manual(
			values = c("TRUE" = "#F8F9FA", "FALSE" = "#F8F9FA")
		) +
		scale_x_discrete(
			labels = function(x) {
				taxname <- gsub("^g__", "", x)
				phylum <- long_data$phylum[match(x, long_data$TaxName)]
				paste(phylum, taxname, sep = " : ")
			}
		) +
		theme_minimal() +
		theme(
			axis.text.x = element_text(
				angle = 90, 
				hjust = 1,
				vjust = 0.5,
				size = 10,
				color = "gray30"
			),
			axis.text.y = element_text(size = 10, color = "gray30"),
			axis.title = element_text(size = 12, color = "gray20", face = "bold"),
			panel.grid.major.x = element_blank(),
			panel.grid.minor = element_blank(),
			panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
			panel.background = element_rect(fill = "white", color = NA),
			plot.background = element_rect(fill = "white", color = NA),
			legend.position = "none"
		) +
		labs(
			x = "Tax Name",
			y = "DNA Damage",
			title = paste("Most abundant", tax_level, "ordered by increasing median damage"))
}

# what we just did above
tax_level = "genus"
pad_tax_path = TRUE
readcutoff = 40000
minreads = 200
minreads_pc1 = 50
min_samples = 3
plot_increasing_mean_damage(tax_level, pad_tax_path, readcutoff, minreads, minreads_pc1, min_samples)

# supplement figs:

# class level
tax_level = "class"
pad_tax_path = TRUE
readcutoff = 30000
minreads = 200
minreads_pc1 = 200
min_samples = 3
plot_increasing_mean_damage(tax_level, pad_tax_path, readcutoff, minreads, minreads_pc1, min_samples)

# looser read cut off at class level
tax_level = "class"
pad_tax_path = TRUE
readcutoff = 3000
minreads = 200
minreads_pc1 = 200
min_samples = 3
plot_increasing_mean_damage(tax_level, pad_tax_path, readcutoff, minreads, minreads_pc1, min_samples)

# looser cut off for genus level
tax_level = "genus"
pad_tax_path = TRUE
readcutoff = 10000
minreads = 200
minreads_pc1 = 200
min_samples = 3
plot_increasing_mean_damage(tax_level, pad_tax_path, readcutoff, minreads, minreads_pc1, min_samples)



###################### Fig 4 part 1. Clustering using top classes ###################### 

# explicitly write out everything to the right of the positive control on the class x damage figure
top_classes = c("c__Anaerolineae","c__JS1","c__Minisyncoccia","c__Phycisphaerae","c__Dehalococcoidia","c__Humimicrobiia",
								"c__UBA8468","c__Patescibacteriia","c__Thermodesulfovibrionia","c__Nanobdellia","c__Nitrospiria","c__Methanosarcinia",
								"c__Koll11","c__Desulfobacteria","c__Nitrososphaeria","c__Cloacimonadia","c__Syntrophia","c__AB16","c__UBA9942")

subg = new[grep(paste0(top_classes,collapse="|"),new$taxpath),] 
subg = subg[,-grep("NC",colnames(subg))]
numsubg = subg[,grep("_Unaggregated",colnames(subg))]
numsubg[is.na(numsubg)] = 0
rownames(numsubg) = subg$TaxName
smallsubg = subg[,(colSums(numsubg) > 200)] # keep this so you can keep the phyla handy
numsubg = numsubg[,colSums(numsubg) > 200] # ditches all NCs, PC1, EM1 and PM2. 22 columns (samples) left
propsubg = sweep(numsubg, 2, colSums(numsubg), `/`) # now it's proportions

# do the nmds
nmds_result <- metaMDS(t(propsubg), distance = "bray", k = 2, trymax = 100)
nmds_df <- as.data.frame(nmds_result$points)
nmds_df <- nmds_df %>%
	mutate(Label = gsub("_UnaggregatedReads", "", rownames(nmds_df)))
nmds_df <- nmds_df %>%
	mutate(SamplePrefix = substr(Label, 1, 2))  # extract first two letters for colours
green_colour <- "#8FBE5C"
purple_colour <- "#7766FF"

my_colors <- c(
	"AI" = purple_colour,
	"MV" = purple_colour,
	"MB" = purple_colour,
	"RM" = purple_colour,
	"EM" = purple_colour,
	"BV" = green_colour,
	"MA" = green_colour,
	"PM" = green_colour,
	"LG" = green_colour
)
site_nmds_ancient <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, label = Label, fill = SamplePrefix, color = SamplePrefix)) +
	geom_point(shape = 21, size = 3, stroke = 0.8) +
	geom_text(
		size = 2,
		nudge_y = .1,
		color = "black"
	) +
	scale_fill_manual(values = my_colors) +
	scale_color_manual(values = my_colors) +
	theme_bw() +
	theme(
		axis.title = element_text(size = 12),
		axis.text = element_text(size = 10),
		legend.position = "none"
	) +
	labs(x = "MDS1", y = "MDS2", title = "(A) Subglacial taxonomic abundance NMDS")
site_nmds_ancient

# plot the loadings for fun
species_scores <- as.data.frame(scores(nmds_result, display = "species"))
species_scores$Taxa <- rownames(species_scores)
species_scores <- species_scores %>%
	left_join(
		subg %>%
			select(TaxName, taxpath, taxlevel),
		by = c("Taxa" = "TaxName")
	) %>%
	mutate(
		Class = sub(".*(c__[^;]+).*", "\\1", taxpath)  # grep class from taxpath
	)
species_scores <- species_scores %>% filter(!is.nan(NMDS1))
class_order <- species_scores %>%
	group_by(Class) %>%
	summarise(mean_NMDS1 = mean(NMDS1, na.rm = TRUE)) %>%
	arrange(mean_NMDS1) %>%
	pull(Class)
species_scores$Class <- factor(species_scores$Class, levels = class_order)

# filter: remove subspecies and taxa with < 100 total reads just to clean it all up 
species_scores = species_scores[which(species_scores$taxlevel != "subspecies"),]
total_reads_lookup <- subg[, c("TaxName", "TotalReads")]  # column 3 is TotalReads
colnames(total_reads_lookup)[2] <- "TotalReads"
species_scores <- species_scores %>%
	left_join(total_reads_lookup, by = c("Taxa" = "TaxName")) %>%
	filter(TotalReads >= 100)

taxa_nmds <- ggplot(species_scores, aes(x = NMDS1, y = Class)) +
	geom_jitter(width = 0, height = 0.2, alpha = 0.5,color="steelblue",aes(size=TotalReads)) +
	theme_bw() +
	labs(x = "NMDS1 score", y = "Class", color = "Taxonomic level",
			 title = "NMDS taxa loadings by class") +
	theme(axis.text.y = element_text(size = 8))

taxa_nmds


# let's do abundances now

top_classes_without_lower_orders = c("c__Anaerolineae","c__JS1","c__Humimicrobiia",
								"c__UBA8468","c__Thermodesulfovibrionia","c__Nanobdellia","c__Nitrospiria","c__Methanosarcinia",
								"c__Koll11","c__Desulfobacteria","c__Nitrososphaeria","c__Cloacimonadia","c__Syntrophia","c__AB16","c__UBA9942")

# get the class ones first
classsubg = subg[subg$taxlevel=="class",]
rownames(classsubg) = classsubg$TaxName 
csubg = classsubg[,grep("_TotalReads",colnames(classsubg))] 
csubg[csubg < 25] <- 0

cpropsubg = sweep(csubg, 2, colSums(csubg), `/`) # now it's proportions
cpropsubg <- cpropsubg %>%
	mutate(Taxa = rownames(cpropsubg))
cpropsubg_long <- cpropsubg %>%
	pivot_longer(
		cols = -Taxa,
		names_to = "Sample",
		values_to = "Abundance"
	) %>%
	mutate(Sample = gsub("_TotalReads", "", Sample)) %>%  # clean names
	filter(!is.na(Sample))  # remove any NAs

# hardcoded sample order , separator is an empty sample to separate columns
sample_order <- c("AIS1","MV1","MV2","MBL1","MBL2","RM1","RM2","EM2","EM3","EM4",
									" ","BV1","BV2","BV3","MA1","MA2","PM1",
									"LG1","LG2","LG3","LG4","LG5","LG6","LG7")

# Keep only samples in your order
cpropsubg_long <- cpropsubg_long %>%
	filter(Sample %in% sample_order)
dummy <- data.frame(Taxa = " ", Sample = " ", Abundance = 0)
cpropsubg_long <- bind_rows(cpropsubg_long, dummy)
cpropsubg_long$Sample <- factor(cpropsubg_long$Sample, levels = sample_order)
taxa_phylum_map <- subg %>% 
	select(TaxName, phylum)
cpropsubg_long <- cpropsubg_long %>%
	left_join(taxa_phylum_map, by = c("Taxa" = "TaxName")) %>%
	mutate(Taxa = ifelse(Taxa == " ", "",
											 paste0(phylum, " : ", gsub("^c__", "", Taxa)))) %>%
	select(-phylum)


palette <- c("white","darkgreen","#85C285","#8E0F9C","#BF3EFF","pink","#018571",
						 "darkslategrey","#CC5500","hotpink","#CDDC39",
						 "#FFB230","deepskyblue","red","yellow","#D7BFE3","lightblue",
						 "#4B1D91","salmon","#1F78B4","deepskyblue","white")


props = ggplot(cpropsubg_long, aes(x = Sample, y = Abundance, fill = Taxa)) +
	geom_bar(stat = "identity", color = "black", size = 0.2) +  # thin black outlines
	scale_y_continuous(labels = scales::percent) +
	scale_fill_manual(values = palette) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1),
				legend.position = "bottom",
				legend.box = "vertical",
				legend.text = element_text(size = 8)) +
	guides(fill = guide_legend(ncol = 3)) +
	labs(y = "Relative abundance", x = "Sample", fill = "Class", title="(C) Subglacial class abundances")
props

site_nmds_ancient / props   + plot_layout(heights = c(1, 2.5))

# now plot cliff's stuff
green_colour <- "#8FBE5C"
purple_colour <- "#7766FF"
oxymetag_clean <- oxymetag %>%
	mutate(Label = sampleID)
plot_df <- inner_join(nmds_df, oxymetag_clean, by = "Label")

oxyplot = ggplot(plot_df, aes(x = MDS1, y = Per_aerobe, label = Label, color = Cluster, fill = Cluster)) +
	geom_point(shape = 21, size = 2.5, stroke = 0.8) +
	scale_color_manual(values = c("Left" = purple_colour, "Right" = green_colour)) +
	scale_fill_manual(values = c("Left" = purple_colour, "Right" = green_colour)) +
#	geom_text_repel(size = 3, color = "black", max.overlaps = 100) +
	theme_bw() +
	labs(x = "MDS1", y = "Predicted aerobes (%)", title="(B) Predicted percent of aerobes")















###################### Fig 4 part 2: Placing Whillans, Mercer and Blood Falls into clusters with DNA ######################

# prep: define top taxa from other papers

# MERCER TOP TAX
# using the single cell paper
# https://www.researchsquare.com/article/rs-4392950/v1
# has in fig 1 the abundant genera and phyla
# here they are, luckily all matching between the gtdb 214 taxonomy they used and the 226 i'm using
mercer_abundant_genera <- list(
	"g__SYFI01",
	"g__F1-60-MAGs163",
	"g__F1-60-MAGs149",
	"g__UBA4592",
	"g__UBA10799",
	"g__Planktophila",
	"g__UBA3006",
	"g__RBG-16-66-20",
	"g__PALSA-1004",
	"g__Polaromonas",
	"g__Nitrotoga",
	"g__39-52-133",
	"g__SURF-13",
	"g__12-FULL-67-14b",
	"g__SPCO01",
	"g__C7867-001",
	"g__UBA1550",
	"g__Nitrosarchaeum"
)
# also these above genus level:
#	"f__Nanopelagicaceae",  
#	"f__Burkholderiaceae",  
#	"f__GW2011-AR1"

#mercer_abundant_phyla <- list(
#	"p__Actinobacteriota",
#	"p__Proteobacteria",
#	"p__Acidobacteriota",
#	"p__Bacteroidota",
#	"p__Chloroflexota",
#	"p__Patescibacteria"
#)
mercer_taxa = c(mercer_abundant_genera)





# WHILLANS TOP TAXA
# https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2016.01457/full fig 2 here

# the whillans taxa with more than 90% identity are:
#	"Sideroxydans lithotrophicus",
#	"Feriphaselus amnicola",
#	"Albidiferax ferrireducens",
#	"Thiobacillus denitrificans",
#	"Acidiferrobacter thiooxydans",
#	"Candidatus Nitrotoga arctica",
#	"Jettenia asiatica",
#	"Nitrosoarchaeum koreensis",
#	"Nitrosospira multiformis",
#	"Methylobacter tundripaludum",
#	"Methylobacillus glycogenes",
#	"Methyloversatilis thermotolerans",
#	"Polaromonas glacialis",
#	"Aggregicoccus edonensis",
#	"Solitalea koreensis",
#	"Ohtaekwangia koreensis",
#	"Acinetobacter lwoffii",
#	"Smithella propionica",
#	"Candidatus Planktophila limnetica",
#	"Ignavibacterium album", -> "f__Ignavibacteriaceae", 
#	"Ilumatobacter fluminis",
#	"Dietzia alimentaria",
#	"Thermomarininilacea lacunofontalis"

# i want to move everything up to genus level. if there is no exact analog in gtdb, here are my chosen mappings:
# Sideroxydans -> Sideroxyarcus
# Feriphaselus is not in GTDB
# Albidiferax is in the genus Rhodoferax
# Nitrososarchaeum -> Nitrosarchaeum
# Jettenia is a real gtdb thing, I just don't see it in our data
# Aggregicoccus is a real gtdb thing, I just don't see it
# Solitalea is real and I don't see it
# Ohtaekwangia likewise
# Ignavibacterium -> Ignavibacteriaceae, but it's a family in gtdb
# 
whillans_taxa <- c(
	"g__Sideroxyarcus", "g__Rhodoferax", "g__Thiobacillus", "g__Acidiferrobacter",
	"g__Nitrotoga", "g__Jettenia", "g__Nitrosarchaeum", "g__Nitrosospira",
	"g__Methylobacter", "g__Methylobacillus", "g__Methyloversatilis",
	"g__Polaromonas", "g__Aggregicoccus", "g__Solitalea", "g__Ohtaekwangia",
	"g__Acinetobacter", "g__Smithella", "g__Planktophila", "g__Ilumatobacter",
	"g__Dietzia", "g__Ignavibacterium"
)




# BLOOD FALLS TOP TAXA
# https://enviromicro-journals.onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.14607 
# fig 4
# this is much harder because somehow it's to a very different db and the names are very different.
# but i'll do my best to map them.

# use the ones in the EMB specifically. these are the guys determined to be a biomarker for EMB (end member output). 
# they have
# Xanthomonadaceae_unc # > f__Xanthomonadaceae
# Thiomicrospira # > f__Thiomicrospiraceae 
# NB1-n_ge # > c__Candidatus Izimiplasma
# Bacteria_unc # SKIP
# DHVE6_ge  # > o__Woesearchaeales and o__Pacearchaeales , both within c__Nanobdellia
# Lutibacter # > g__Lutibacter
# Deltaproteobacteria_unc # SKIP. Deltaproteobacteria is a class in NCBI that maps to different things in gtdb... not even to one single phylum. i'm gonna skip it.
# Desulfobacteraceae_unc # > f__Desulfobacteraceae 
# Geopsychrobacter # > f__Geopsychrobacteraceae
# Desulfobulbaceae_unc # > f__Desulfocapsaceae
# Bacteroidetes_unc # > p__Bacteroidota in gtdb; wikipedia says "The phylum Bacteroidota (synonym Bacteroidetes)"
	# ditch this, it's in everything 
# Sva0996_marine_group_ge # > s__marine group bacterium Sva0996 bin134 . bump to gtdb taxonomy genus level for this which is g__Poriferisocius
# Maritimimonas # > o__Flavobacteriales (no analog in gtdb, wiki says it's a "Flavobacteriaceae", which in ncbi maps to many things in gtdb, all which seem to be in this order)
# Atribacteria_ge # > p__Atribacterota
# f__Izemoplasmataceae

# logic:

# Was hard to figure out what NB1-n_ge is. The paper says "The Tenericutes were represented by OTU45, which classified in the NB1-n order of the Mollicutes class."
# and then cites this https://academic.oup.com/ismej/article-abstract/10/11/2679/7538138
# which calls them Izimiplasma, which DOES have an analog in gtdb 226. f__Izemoplasmataceae

# Blood falls paper says ""The   SILVA (v128) reference taxonomy contains sequences from   the DHVE6 archaeal group classiﬁed as Woesearchaeota   but also containing sequences closer to Pacearchaeota""
# so i will map them to both I guess

bloodfalls_taxa <- c("g__Lutibacter","g__Poriferisocius","f__Xanthomonadaceae","f__Thiomicrospiraceae",
										 "f__Izemoplasmataceae","f__Desulfobacteraceae", "f__Geopsychrobacteraceae",
										 "f__Desulfocapsaceae",  "o__Woesearchaeales","o__Pacearchaeales",
										 "o__Flavobacteriales", "p__Atribacterota")


extract_taxa_abundance <- function(data, taxa_list, total_reads) {
	sample_cols <- grep("_TotalReads$", colnames(data), value = TRUE)
	sample_names <- gsub("_TotalReads$", "", sample_cols)
	
	result_matrix <- matrix(0, nrow = length(taxa_list), ncol = length(sample_names))
	rownames(result_matrix) <- taxa_list
	colnames(result_matrix) <- sample_names
	
	for (i in seq_along(taxa_list)) {
		taxon <- taxa_list[i]
		matching_rows <- grep(taxon, data$TaxName, fixed = TRUE)
		
		if (length(matching_rows) > 0) {
			for (j in seq_along(sample_cols)) {
				sample_col <- sample_cols[j]
				sample_name <- sample_names[j]
				
				total_reads_taxon <- sum(data[matching_rows, sample_col], na.rm = TRUE)
				
				if (sample_name %in% names(total_reads)) {
					result_matrix[i, j] <- total_reads_taxon / total_reads[sample_name]
				}
			}
		}
	}
	
	return(result_matrix)
}

# colour palettes for each group
reds <- c("white", "#EE8E8E", "#D95C5C", "#8B0000", "#5A0000", "black")
purple_blues <- c("white", "#F0E6FF", "#B399FF", "#7766FF", "#5544CC", "#3D2D99")
greens <- c("white",   "#8FBE5C", "#6B9E3E", "#1B5E20","#1C5F2F","black")

create_heatmap_with_separator <- function(abundance_matrix, title, sample_order) {
	plot_data <- as.data.frame(abundance_matrix)
	plot_data$Taxon <- rownames(plot_data)
	
	plot_long <- pivot_longer(plot_data, 
														cols = -Taxon, 
														names_to = "Sample", 
														values_to = "Proportion")
	
	# Remove NA samples
	plot_long <- plot_long[!is.na(plot_long$Sample), ]
	
	# Define the three groups
	group1_samples <- c("PC1")
	group2_samples <- c("AIS1","MV1","MV2","MBL1","MBL2","RM1","RM2","EM2","EM3","EM4")
	
	plot_long$Group <- ifelse(plot_long$Sample %in% group1_samples, "Group1",
														ifelse(plot_long$Sample %in% group2_samples, "Group2", "Group3"))
	
	sample_order_clean <- sample_order[sample_order != " " & !is.na(sample_order)]
	
	# Only keep samples that exist in sample_order_clean
	plot_long <- plot_long[plot_long$Sample %in% sample_order_clean, ]
	
	plot_long$Sample <- factor(plot_long$Sample, levels = sample_order_clean)
	
	# Create separate data for each group
	plot_g1 <- plot_long[plot_long$Group == "Group1", ]
	plot_g2 <- plot_long[plot_long$Group == "Group2", ]
	plot_g3 <- plot_long[plot_long$Group == "Group3", ]
	
	# Set hard max at 7.5%
	scale_max <- 0.1
	
	p <- ggplot() +
		geom_tile(data = plot_g1, aes(x = Sample, y = Taxon, fill = Proportion), color = "white", linewidth = 0.5) +
		scale_fill_gradientn(colors = reds,
												 limits = c(0, scale_max),
												 values = c(0, 0.005, 0.01, 0.05, 0.1, 1),
												 name = "Proportion",
												 labels = scales::percent,
												 aesthetics = "fill") +
		ggnewscale::new_scale_fill() +
		geom_tile(data = plot_g2, aes(x = Sample, y = Taxon, fill = Proportion), color = "white", linewidth = 0.5) +
		scale_fill_gradientn(colors = purple_blues,
												 limits = c(0, scale_max),
												 values = c(0, 0.005, 0.01, 0.05, 0.1, 1),
												 name = "Proportion",
												 labels = scales::percent) +
		ggnewscale::new_scale_fill() +
		geom_tile(data = plot_g3, aes(x = Sample, y = Taxon, fill = Proportion), color = "white", linewidth = 0.5) +
		scale_fill_gradientn(colors = greens,
												 limits = c(0, scale_max),
												 values = c(0, 0.005, 0.01, 0.05, 0.1, 1),
												 name = "Proportion",
												 labels = scales::percent) +
		theme_linedraw() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
					axis.text.y = element_text(size = 7),
					panel.grid = element_blank(),
					plot.title = element_text(hjust = 0.5, face = "bold")) +
		labs(title = title, x = "Sample", y = "Taxon") +
		geom_vline(xintercept = 1.5, color = "black", linewidth = .7, linetype = "solid") +
		geom_vline(xintercept = 11.5, color = "black", linewidth = .7, linetype = "solid")
	
	return(p)
}

sample_order <- c("PC1"," ","AIS1","MV1","MV2","MBL1","MBL2","RM1","RM2","EM2","EM3","EM4",
									" ","BV1","BV2","BV3","MA1","MA2","PM1",
									"LG1","LG2","LG3","LG4","LG5","LG6","LG7")

sample_cols <- grep("_TotalReads$", colnames(new), value = TRUE)
sample_names <- gsub("_TotalReads$", "", sample_cols)
total_assigned_reads_per_sample <- sapply(sample_cols, function(col) sum(new[[col]], na.rm = TRUE))
names(total_assigned_reads_per_sample) <- sample_names

sample_order_clean <- sample_order[sample_order != " "]
missing_samples <- setdiff(sample_order_clean, sample_names)
if (length(missing_samples) > 0) {
	warning("Some samples in sample_order not found in data: ", paste(missing_samples, collapse = ", "))
}
available_samples <- intersect(sample_order_clean, sample_names)
cat("Processing", length(available_samples), "samples\n")

mercer_abundance <- extract_taxa_abundance(new, mercer_taxa, total_assigned_reads_per_sample)
whillans_abundance <- extract_taxa_abundance(new, whillans_taxa, total_assigned_reads_per_sample)
bloodfalls_abundance <- extract_taxa_abundance(new, bloodfalls_taxa, total_assigned_reads_per_sample)

mercer_heatmap <- create_heatmap_with_separator(mercer_abundance, 
																								"Relative abundances of top Subglacial Lake Mercer taxa in our samples", 
																								sample_order)
whillans_heatmap <- create_heatmap_with_separator(whillans_abundance, 
																									"Relative abundances of top Subglacial Lake Whillans taxa in our samples", 
																									sample_order)
bloodfalls_heatmap <- create_heatmap_with_separator(bloodfalls_abundance, 
																										"Relative abundances of top Blood Falls taxa in our samples", 
																										sample_order)
heatmaps <- mercer_heatmap / whillans_heatmap / bloodfalls_heatmap +
	plot_layout(heights = c(2, 2, 1.2))

nmdsabundance <- site_nmds_ancient / oxyplot / props +
	plot_layout(heights = c(1, 1, 2))   # NMDS half the size of props

final_plot <- nmdsabundance | heatmaps +
	plot_layout(widths = c(1, 1.5))  # adjust column widths if desired

final_plot





###################### Fig 5: Geochemistry ######################

subglacial_sheet = all_geochem_sheets$`Antarctic precipitate data`
# colour by $`0= oxidized / 1 =  reduced`
eh_sheet = all_geochem_sheets$`Fe_Mn_Eh_Barnaby_&_Rimstidt`
riggs_sheet = all_geochem_sheets$Riggs_etal_2022
green_sheet = all_geochem_sheets$Green_etal_2015
sinha_sheet = all_geochem_sheets$`Sinha et al. 2018`
sirona_sheet = all_geochem_sheets$Sironi_etal_2023
drake_sheet = all_geochem_sheets$Drake_etal_2017
kusturica_sheet = all_geochem_sheets$Kusturica_etal_2022
yuguchi_sheet = all_geochem_sheets$Yuguchi_etal_2022
wang_sheet = all_geochem_sheets$Wang_etal_2020
earthchem_sheet = all_geochem_sheets$Earthchem_FeMnOnly


riggs_df <- riggs_sheet %>%
	select(Fe = 1, Mn = 2) %>%
	slice(-1) %>%  
	mutate(Fe = as.numeric(Fe), Mn = as.numeric(Mn)) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

green_df <- green_sheet %>%
	select(Fe, Mn) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

sinha_df <- sinha_sheet %>%
	select(Fe, Mn) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

sirona_df <- sirona_sheet %>%
	select(Fe = Fe_ppm, Mn = Mn_ppm) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

drake_df <- drake_sheet %>%
	select(Fe, Mn) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

kusturica_df <- kusturica_sheet %>%
	select(Fe = `Fe (ppm)`, Mn = `Mn (ppm)`) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

yuguchi_df <- yuguchi_sheet %>%
	select(Fe, Mn) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

wang_df <- wang_sheet %>%
	select(Fe = `Fe (ppm)`, Mn = `Mn (ppm)`) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

earthchem_df <- earthchem_sheet %>%
	select(Fe, Mn) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0)

background_df <- bind_rows(
	riggs_df %>% mutate(Source = "Riggs et al. 2022"),
	green_df %>% mutate(Source = "Green et al. 2015"),
	sinha_df %>% mutate(Source = "Sinha et al. 2018"),
	sirona_df %>% mutate(Source = "Sironi et al. 2023"),
	drake_df %>% mutate(Source = "Drake et al. 2017"),
	kusturica_df %>% mutate(Source = "Kusturica et al. 2022"),
	yuguchi_df %>% mutate(Source = "Yuguchi et al. 2022"),
	wang_df %>% mutate(Source = "Wang et al. 2020"),
	earthchem_df %>% mutate(Source = "Earthchem")
)

subglacial_individual <- subglacial_sheet %>%
	select(sample_name = ...1,
				 Fe = `Fe (ppm)`, 
				 Mn = `Mn (ppm)`, 
				 oxidation_state = `0= oxidized / 1 =  reduced`) %>%
	filter(!is.na(Fe) & !is.na(Mn) & Fe > 0 & Mn > 0) %>%
	filter(sample_name != "EM1")

subglacial_df <- subglacial_individual %>%
	group_by(sample_name, oxidation_state) %>%
	summarize(
		Fe = mean(Fe, na.rm = TRUE),
		Mn = mean(Mn, na.rm = TRUE),
		.groups = "drop"
	)

eh_df = barnaby_grab
colnames(eh_df) = c("Fe","Mn","Eh")

source_mn_avg <- background_df %>%
	group_by(Source) %>%
	summarize(avg_Mn = mean(Mn, na.rm = TRUE)) %>%
	mutate(color_group = ifelse(avg_Mn < 100, "cool", "warm"))

green_colour <- "#8FBE5C"
purple_colour <- "#7766FF"

source_colors <- source_mn_avg %>%
	arrange(avg_Mn) %>%
	mutate(
		source_color = case_when(
			color_group == "cool" & row_number() == 1 ~ "#C89080",  # Riggs et al. 2022
			color_group == "cool" & row_number() == 2 ~ "#E8C4A0",  # Sinha et al. 2018
			color_group == "cool" & row_number() == 3 ~ "#C5B8D4",  # Sironi et al. 2023
			color_group == "cool" & row_number() == 4 ~ "#E8B8C8",  # Wang et al. 2020
			color_group == "cool" & row_number() == 5 ~ "#E8A87C",  # Kusturica et al. 2022
			color_group == "warm" & row_number() == 6 ~ "#A0B5C8",  # Green et al. 2015
			color_group == "warm" & row_number() == 7 ~ "#A0E6E1",  # Drake et al. 2017
			color_group == "warm" & row_number() == 8 ~ "#7FCBD0",  # Yuguchi et al. 2022
			color_group == "warm" & row_number() == 9 ~ "#74C69D", # Earthchem
			TRUE ~ "#CCCCCC"
		)
		,
		bins = case_when(
			Source == "Riggs et al. 2022" ~ 3,
			Source == "Sinha et al. 2018" ~ 3,
			Source == "Sironi et al. 2023" ~ 3,
			Source == "Wang et al. 2020" ~ 3,
			Source == "Kusturica et al. 2022" ~ 3,
			Source == "Green et al. 2015" ~ 4,
			Source == "Drake et al. 2017" ~ 4,
			Source == "Yuguchi et al. 2022" ~ 4,
			Source == "Earthchem" ~ 5,
			TRUE ~ 4
		)
	)

background_df <- background_df %>%
	left_join(source_colors, by = "Source")

background_large <- background_df %>%
	group_by(Source) %>%
	filter(n() >= 10) %>%
	ungroup()

background_small <- background_df %>%
	group_by(Source) %>%
	filter(n() < 10) %>%
	ungroup()

source_list <- split(background_large, background_large$Source)

# fit bfb and slw
eh_interp_Mn <- approxfun(eh_df$Eh, eh_df$Mn)
eh_interp_Fe <- approxfun(eh_df$Eh, eh_df$Fe)
eh_382_Mn <- eh_interp_Mn(382)
eh_382_Fe <- eh_interp_Fe(382)
eh_90_Mn <- eh_interp_Mn(90)
eh_90_Fe <- eh_interp_Fe(90)
special_points <- data.frame(
	Eh = c(382, 90),
	Mn = c(eh_382_Mn, eh_90_Mn),
	Fe = c(eh_382_Fe, eh_90_Fe),
	color = c("white", "black")
)




p <- ggplot()

for (src in names(source_list)) {
	src_data <- source_list[[src]]
	src_color <- unique(src_data$source_color)[1]
	p <- p + geom_point(data = src_data,
											aes(x = Mn, y = Fe),
											color = src_color,
											alpha = 0.2, size = 1)
}

for (src in names(source_list)) {
	src_data <- source_list[[src]]
	src_color <- unique(src_data$source_color)[1]
	src_bins <- source_colors %>% filter(Source == src) %>% pull(bins)
	
	p <- p + stat_density_2d(data = src_data,
													 aes(x = Mn, y = Fe),
													 geom = "polygon",
													 fill = src_color,
													 alpha = 0.15,
													 color = NA,
													 bins = src_bins,
													 contour_var = "ndensity")
	
	p <- p + stat_density_2d(data = src_data,
													 aes(x = Mn, y = Fe, alpha = after_stat(level)),
													 geom = "density_2d",
													 color = src_color,
													 linewidth = 0.3,
													 bins = src_bins,
													 contour_var = "ndensity",
													 show.legend = FALSE)
}

p <- p +
	geom_point(data = background_small,
						 aes(x = Mn, y = Fe),
						 color = "grey60", alpha = 0.5, size = 1) +
	geom_path(data = eh_df, 
						aes(x = Mn, y = Fe, color = Eh),
						linewidth = 3, alpha = 1) +
	geom_point(data = eh_df, 
						 aes(x = Mn, y = Fe, color = Eh),
						 size = 1.5, alpha = 0.9) +
	geom_point(data = subglacial_individual,
						 aes(x = Mn, y = Fe, fill = factor(oxidation_state)),
						 shape = 21, size = 3, alpha = 0.4, color = "black", stroke = 0.6) +
	geom_point(data = subglacial_df, 
						 aes(x = Mn, y = Fe, fill = factor(oxidation_state)),
						 shape = 21, size = 4, alpha = 1, color = "black", stroke = 0.6) +
	geom_label(data = subglacial_df,
						 aes(x = Mn, y = Fe, label = sample_name),
						 size = 2.2, 
						 nudge_y = 0.1,
						 label.padding = unit(0.15, "lines"),
						 label.size = 0,
						 fill = "white",
						 alpha = 0.8,
						 color = "grey10", fontface = "bold") +
	geom_point(data = special_points,
						 aes(x = Mn, y = Fe),
						 shape = 4,
						 size = 2,
						 stroke = 1,
						 color = special_points$color)+
	scale_x_log10(
		limits = c(0.3, 35000),
		breaks = 10^(0:4),
		labels = trans_format("log10", math_format(10^.x)),
		expand = expansion(mult = c(0.02, 0.05))
	) +
	scale_y_log10(
		limits = c(1, 40000),
		breaks = 10^(0:4),
		labels = trans_format("log10", math_format(10^.x)),
		expand = expansion(mult = c(0.02, 0.05))
	) +
	scale_color_viridis_c(
		name = "Eh (mV)",
		option = "viridis",
		limits = c(40, 400),
		direction= -1
	) +
	scale_fill_manual(
		name = "Subglacial Precipitates",
		values = c("0" = purple_colour, "1" = green_colour),
		labels = c("0" = "Oxidized", "1" = "Reduced")
	) +
	labs(
		x = "Mn (ppm)",
		y = "Fe (ppm)",
		title = "Fe-Mn concentrations"
	) +
	theme_bw(base_size = 13) +
	theme(
		legend.position = "right",
		legend.box = "vertical",
		panel.grid.major = element_line(color = "grey90", size = 0.3),
		panel.grid.minor = element_blank(),
		panel.border = element_rect(color = "grey60", size = 0.8),
		plot.title = element_text(size = 14, face = "bold"),
		axis.title = element_text(size = 12),
		legend.title = element_text(size = 11, face = "bold"),
		legend.text = element_text(size = 10),
		legend.key.height = unit(1.2, "cm")
	) +
	annotation_logticks()

print(p)



 


