# -------------------------------------------------------------------------
# Astra PheWAS
# 12Jan2-7Feb22 / wjst
# Exon variants associated with asthma and allergy
# -------------------------------------------------------------------------

usehere <- c("Gviz","biomaRt","dplyr","DBI","rtracklayer","BSgenome.Hsapiens.UCSC.hg38","ggVennDiagram","GenomicRanges","ggplot2","ggrepel","grDevices","IRanges","RSQLite","stringr","tibble","tidyr","utils")
lapply(usehere, require, character.only = TRUE)
work <- c("~/Desktop") # or getwd()

# not_run!
# if there are errors in magittr chains, use the latest dplyr version (deinstall tidyr) and prepend the package dplyr::command
# *** 11Feb22 open issue: the allelic dataset from azphenas is missing some variants 
# *** found in the dominant or recessive dataset
# *** until resolved use scramble2()

# functions ---------------------------------------------------------------

# data cleaning
scramble <- function(df) {
	print( dim(df) )
	words <- paste(c("_variant", "_region_variant","start_codon_","_transcript_exon"), collapse = "|")
	res <- df %>%
		separate(Variant,sep="-",into = c("chromosome_name","pos","alleleA","alleleB")) %>%
		left_join(genepos, by = c("chromosome_name" = "chromosome_name")) %>%
		mutate(Consequence = gsub(words, "\\1", Consequence.type) ) %>%
		mutate(id=paste(chromosome_name,pos,Gene,Consequence,sep="_")) %>%
		mutate(chr = as.integer( case_when(chromosome_name=="X"~"23",chromosome_name=="Y"~"24",TRUE~chromosome_name) ) ) %>%
		filter(Model=="allelic") %>%
		mutate(pos =as.numeric(pos)) %>%
		mutate(cumpos = pos+cumpos ) %>%
		mutate(p = -log10(p.value)) %>%    # P≤2×10-9 in Wang = 0.000000002; log10 is 8.69897
		mutate(p = case_when(p>30 ~ 30, TRUE ~ p)) %>%
		mutate_all(function(x) ifelse(is.infinite(x),30, x)) %>%
		arrange(chr,pos)
	print( dim(res) )
  return(res)
}

scramble2 <- function(df) {
	print( dim(df) )
	words <- paste(c("_variant", "_region_variant","start_codon_","_transcript_exon"), collapse = "|")
	res <- df %>%
		separate(Variant,sep="-",into = c("chromosome_name","pos","alleleA","alleleB")) %>%
		left_join(genepos, by = c("chromosome_name" = "chromosome_name")) %>%
		mutate(Consequence = gsub(words, "\\1", Consequence.type) ) %>%
		mutate(id=paste(chromosome_name,pos,Gene,Consequence,sep="_")) %>%
		mutate(chr = as.integer( case_when(chromosome_name=="X"~"23",chromosome_name=="Y"~"24",TRUE~chromosome_name) ) ) %>%
		mutate(pos =as.numeric(pos)) %>%
		mutate(cumpos = pos+cumpos ) %>%
		mutate(p = -log10(p.value)) %>%    # P≤2×10-9 in Wang = 0.000000002; log10 is 8.69897
		mutate(p = case_when(p>30 ~ 30, TRUE ~ p))
	rebuild <- res %>%
		distinct(id) %>%
		left_join(res %>% dplyr::select(-Model,-p)  %>% group_by(id) %>% filter(row_number()==1), by="id") %>%
		left_join(res %>% filter(Model=="allelic")  %>% dplyr::select(id,p) %>% dplyr::rename(p_all=p), by="id") %>%
		left_join(res %>% filter(Model=="dominant") %>% dplyr::select(id,p) %>% dplyr::rename(p_dom=p), by="id") %>%
		left_join(res %>% filter(Model=="recessive")%>% dplyr::select(id,p) %>% dplyr::rename(p_rec=p), by="id") %>%
		mutate(p= purrr::pmap_dbl(list(p_all, p_dom, p_rec, na.rm=TRUE), max)) %>%
		arrange(chr,pos)
	print( dim(rebuild) )
	return(rebuild)
}

# show genomic positions
spot_on <- function(chr,from,to,gene) {
	cchr = paste0("chr",chr)
	genome = "hg38"
	lift <- liftOver( GRanges(seqnames=cchr,ranges=IRanges(from,to),strand=c("*")), chain)@unlistData@ranges
	ff = start(lift)[1]
	tt = end(lift)[1]
	
	sql <- paste0("SELECT RS_ID, CHR, BP, BETA, PVALUE, EFFECT_ALLELE FROM res where CHR is '",chr,"' AND BP>",ff," AND BP<",tt,";")
	db1 <- DBI::dbGetQuery(conGWAS,sql) %>%
		mutate(p = -log10(PVALUE)) %>%
		mutate(p = case_when(p>30 ~ 30, TRUE ~ p)) %>%
		mutate_all(function(x) ifelse(is.infinite(x),30, x))
		
	sql <- paste0("SELECT * FROM ld where CHR is '",chr,"' AND BP>",ff," AND BP<",tt,";")
	db2 <- DBI::dbGetQuery(conLD,sql)

	#  GRCh37=hg19; current version is GRCh38.p13, always attn to offset
	grch38.snp = useMart(biomart="ENSEMBL_MART_SNP",host="https://dec2021.archive.ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")
	if ( dim(db1)[1]>0 ) {
		ens <- getBM(attributes = c("refsnp_id", "chrom_start"), filters = "snp_filter", values = db1$RS_ID, mart = grch38.snp)
		db1 <- db1 %>%
			left_join(ens, by = c("RS_ID"="refsnp_id")) %>%
			filter( complete.cases(.) )
	}
	if ( dim(db2)[1]>0 ) {
		ens <- getBM(attributes = c("refsnp_id", "chrom_start"), filters = "snp_filter", values = db2$SNP, mart = grch38.snp)
		db2 <- db2 %>%
			left_join(ens, by = c("SNP"="refsnp_id")) %>%
			filter( complete.cases(.) )
	}
	
  itrack <- IdeogramTrack(genome = genome, chromosome = cchr)
	
	gtrack <- GenomeAxisTrack(genome = genome, chromosome = cchr, fontsize=25)
	
	# strack <- SequenceTrack(Hsapiens)
	
	btrack <- BiomartGeneRegionTrack(genome = genome, chromosome = cchr, start = from, end = to, name = "transcripts", fontsize=25, cex.group=1.1)

	dtrack1ar <-  AnnotationTrack(data=1, start=1, width=1, name="A", id=1, genome=genome, chromosome=cchr, fontsize=25)
	dtrack1as <-  AnnotationTrack(data=1, start=1, width=1, name="AR", id=1, genome=genome, chromosome=cchr, fontsize=25)
	dtrack1ad <-  AnnotationTrack(data=1, start=1, width=1, name="AD", id=1, genome=genome, chromosome=cchr, fontsize=25)
	tmp <- ar %>%
		filter(p>8.7) %>%
		filter(chr==!!chr) %>%
		mutate(color = case_when(Odds.ratio<1 ~"green",Odds.ratio>=1 ~"red"))
	if (dim(tmp)[1] >0) {
		dtrack1ar <-  AnnotationTrack(start = tmp[,"pos"], width=50, name="A", group = paste(tmp[,"cDNA.change"],tmp[,"Consequence.type"],sep="\n"), 
								  stacking="squish", id=tmp[,"pos"], fill=tmp[,"color"], col=tmp[,"color"], shape = "box", genome=genome, chromosome=cchr, cex=0,
								  fontcolor.group = "gray", lineheight.group = 0.8, fontsize=25, cex.group=1.1)
		              #tmp[,"NO..cases"],"/",[,"NO..cases"]," ",round(tmp[,"Case.MAF"],4),"/",round(tmp[,"Control.MAF"],4)
	} 
	tmp <- as %>%
		filter(p>8.7) %>%
		filter(chr==!!chr) %>%
		mutate(color = case_when(Odds.ratio<1 ~"green",Odds.ratio>=1 ~"red"))
	if (dim(tmp)[1] >0) {
		dtrack1as <-  AnnotationTrack(start = tmp[,"pos"], width=50, name="AR", group = paste(tmp[,"cDNA.change"],tmp[,"Consequence.type"],sep="\n"),
								  stacking="squish", id=tmp[,"pos"], fill=tmp[,"color"], col=tmp[,"color"], shape = "box", genome=genome, chromosome=cchr, cex=0,
								  fontcolor.group = "gray", lineheight.group = 0.8, fontsize=25, cex.group=1.1)
	}
	tmp <- ad %>%
		filter(p>8.7) %>%
		filter(chr==!!chr) %>%
		mutate(color = case_when(Odds.ratio<1 ~"green",Odds.ratio>=1 ~"red"))
	if (dim(tmp)[1] >0) {
		dtrack1ad <-  AnnotationTrack(start = tmp[,"pos"], width=50, name="AD", group = paste(tmp[,"cDNA.change"],tmp[,"Consequence.type"],sep="\n"),
								  stacking="squish", id=tmp[,"pos"], fill=tmp[,"color"], col=tmp[,"color"], shape = "box", genome=genome, chromosome=cchr, cex=0,
								  fontcolor.group = "gray", lineheight.group = 0.8, fontsize=2, cex.group=1.1)
	}
	
	dtrack2 <-  AnnotationTrack(data=1, start=1, width=1, name="GWAS", id=1, genome=genome, chromosome=cchr, fontsize=25)
	if ( dim(db1)[1]>0 ) {
		dtrack2 <- DataTrack(data = db1[,"BETA"], start = db1[,"chrom_start"], width=5, name="GWAS", identifier=db1[,"RS_ID"], genome=genome, chromosome=cchr, type="histogram", fill="gray", fontsize=25)
	}
	ldtrack <- DataTrack(data=1, start=1, width=10, name="LD", genome=genome, chromosome=cchr, type="l", fill="red")
	if ( dim(db2)[1]>0 ) {
		ldtrack <- DataTrack(data = db2[,"L2"], start=db2[,"chrom_start"], width=5, name="LD", id=db2[,"SNP"], genome=genome, chromosome=cchr, type=c("p","l"), col="gray", fontsize=25)
	}
	
	fn <- paste0(work,"/",gene,".png")
	png(file=fn, width=1024*2, height=610*2)
	plot.new()
	plotTracks(list(itrack,gtrack,btrack,dtrack1as,dtrack1ar,dtrack2,ldtrack) , from=from, to=to, showId = TRUE, sizes= c(1,1,6,3,3,2,1) )
	dev.off()
}

# data --------------------------------------------------------------------

# calculate genomic positions
chain <- import.chain(paste0(work,"/hg38ToHg19.over.chain"))

# chromosome positions Manhattan plot
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
BMgene <- getBM(attributes = c("external_gene_name","chromosome_name", "start_position","end_position","description"),mart = ensembl)
genepos <- BMgene %>%
	filter(chromosome_name %in% c(1:22,"X","Y")) %>%
	group_by(chromosome_name) %>%
	filter(end_position == max(end_position)) %>%
	dplyr::select(end_position,chromosome_name) %>%
	mutate(chr=as.numeric( case_when(chromosome_name=="X"~"23",chromosome_name=="Y"~"24",TRUE~chromosome_name)) ) %>%
	arrange(chr) %>%
	ungroup() %>%
	mutate(cumpos = cumsum(as.numeric(end_position))) %>%
	dplyr::select(cumpos,chromosome_name)

# GWAS https://genepi.qimr.edu.au/staff/manuelf/gwas_results/intro.html
conGWAS <- dbConnect(RSQLite::SQLite(), paste0(work,"/2/share.sqlite3" )

# LD https://alkesgroup.broadinstitute.org/LDSCORE/
conLD <- dbConnect(RSQLite::SQLite(), "/Users/wjst/Documents/Daten/Projekte/Il2/ldscore.sqlite3")

# allergic rhinitis https://azphewas.com/phenotypeView/7e2a7fab-97f0-45f7-9297-f976f7e667c8/eb827281-c7bc-4852-aaca-07fa69cab580/glr
fn <- paste0(work,"/20002_1387_hayfever_allergic_rhinitis_filtered_table_AZ_PheWAS_Portal.csv")
ar <- read.csv(fn) %>%
	scramble()

# asthma https://azphewas.com/phenotypeView/7e2a7fab-97f0-45f7-9297-f976f7e667c8/2bc2f761-4074-48f5-9106-938ed34db496/glr
fn <- paste0(work,"/Union_J45_J45_Asthma_filtered_table_AZ_PheWAS_Portal.csv")
as <- read.csv(fn) %>%
	scramble()

# atopic dermatitis https://azphewas.com/phenotypeView/7e2a7fab-97f0-45f7-9297-f976f7e667c8/f7885511-e269-434b-acea-f7591bc0a4d1/glr
fn <- paste0(work,"/Union_L20_L20_Atopic_dermatitis_filtered_table_AZ_PheWAS_Portal.csv")
ad <- read.csv(fn) %>%
	scramble()

# Fig 1 genomewide --------------------------------------------------------

p <- ggplot() +
	geom_point(data=as,aes(y=p,x=cumpos),shape=21, fill = "red",size=5) +
	geom_text_repel(data=as,aes(x=cumpos,y=p,label=Gene), nudge_x=2, color="red",max.overlaps=12, nudge_y=1, size=4, box.padding = .5, max.iter=100000) +	

	geom_point(data=ar,aes(y=p,x=cumpos), shape=21, fill = "black",size=5) +
	geom_text_repel(data=ar,aes(x=cumpos,y=p,label=Gene), nudge_x=2, color="black",max.overlaps=12, size=4, box.padding = .5, max.iter=100000) +	
	
	geom_point(data=ad,aes(y=p,x=cumpos), shape=21, fill = "blue",size=5) +	
	geom_text_repel(data=ad,aes(x=cumpos,y=p,label=Gene), nudge_x=2, nudge_y=-.5, color="blue", size=4) +	
	
	labs(y=bquote(-log[10](P)), x="chromosomes") +

	scale_y_continuous( limits=c(8,30) ) +
	
	scale_x_continuous(position = "top") +
	
	geom_vline(xintercept = genepos$cumpos, linetype="dotted", color="gray") +
	annotate("text", x=genepos$cumpos, y=30, label=c(1:20," ","22","X","Y"), size=7 ) +
	
	clean_theme( axis.text.x=element_blank(), axis.text=element_text(size=20)), axis.ticks.x=element_blank() ) 

fn = paste0(work,"/fig1_genome.png")
png(file =  fn, width = 1600, height = 900, units = 'px')
print(p)
dev.off()

# Fig 2 Odds ratios -------------------------------------------------------

rbind(
	ar %>%
		filter(p>8.7) %>%
		mutate(id=paste(str_pad(chr,2,pad="0"),str_pad(pos,9,pad="0"),sep="_")) %>%
		rename(value1=Odds.ratio,group=Gene) %>%
		mutate(value2=0) %>%
		dplyr::select(id,group,value1,value2),
	as %>%
		filter(p>8.7) %>%
		mutate(id=paste(str_pad(chr,2,pad="0"),str_pad(pos,9,pad="0"),sep="_")) %>%
		mutate(value1=0) %>%
		rename(value2=Odds.ratio,group=Gene) %>%
		dplyr::select(id,group,value1,value2)) %>%
	arrange(id) %>%
	rename(individual=id) %>%
	gather(key = "observation", value="value", -c(1,2)) %>%
	mutate(observation = case_when(observation=="value1"~"as",TRUE~"ar") ) %>%
	arrange(individual) %>%
	rowid_to_column(var = "r")	%>%
	mutate(angle = 90 - 360/max(r)*(r-0.5)) %>%
	mutate(hjust = case_when(angle<(-90)~1,TRUE~0)) %>%
	mutate(angle = case_when(angle<(-90)~angle+180,TRUE~angle)) %>%
  (function(x){
	  data <<- left_join(x,	  
	  	x %>%
	  	  group_by(group) %>%
	  	  mutate( first = first(r), last = last(r) ) %>%
  		  filter(row_number()==ceiling(n()/2)) %>%
  		  rename(label=group) %>%
  		  dplyr::select(r,label,first,last), by="r")
  })

p <- ggplot(data, aes(fill=observation, y=value, x=r ) )+ 
	geom_bar(position="dodge", stat="identity") +
	scale_fill_manual(values=c("red","black")) +
	coord_cartesian(ylim=c(-.5,1.6)) +
	geom_hline(yintercept = 1, linetype="dotted", color="black") +
	geom_text_repel( aes(x=r,y=1.5,label=label,angle=angle), size=4, segment.color = NA, box.padding = .5, max.iter=100000) + 
	geom_segment(aes(x=first, y = 1.3, xend=last, yend = 1.3), colour = "grey", alpha=.8, size=1.8 ) +
	theme_minimal() +
	theme(
		legend.position = "none",
		axis.text = element_blank(),
		axis.title = element_blank(),
		panel.grid = element_blank()
	) +
	coord_polar()

fn = paste0(work,"/fig2_polar.png")
png(file =  fn, width = 1600, height = 900, units = 'px')
print(p)
dev.off()

# Fig 3 shared genes ------------------------------------------------------

variants <- rbind(
	ar %>% mutate(trait="ar"),
	as %>% mutate(trait="as")) %>%
	filter(p>8.7) %>% 
	arrange(chromosome_name,pos)

variants.gene.list <- list(
	rhinitis.plus =    variants %>% filter(trait=="ar" & Odds.ratio>1) %>% dplyr::select(Gene) %>% distinct(Gene) %>% pull(),
	asthma.plus =      variants %>% filter(trait=="as" & Odds.ratio>1) %>% dplyr::select(Gene) %>% distinct(Gene) %>% pull(),
	rhinitis.minus =   variants %>% filter(trait=="ar" & Odds.ratio<1) %>% dplyr::select(Gene) %>% distinct(Gene) %>% pull(),
	asthma.minus =     variants %>% filter(trait=="as" & Odds.ratio<1) %>% dplyr::select(Gene) %>% distinct(Gene) %>% pull())
venn = Venn(variants.gene.list)
data = process_data(venn)
cont <- venn_region(data)
out$name <- c("rhinitis\nOR>1","asthma\nOR>1","rhinitis\nOR<1","asthma\nOR<1")

p <- ggplot() +
	scale_fill_gradient(low="grey", high="red") +
	geom_sf(aes(fill=count), data = venn_region(data)) +
	geom_sf(data = venn_setedge(data)) +
	geom_sf_text(aes(label = name), data = out, size=5) +
	geom_sf_label(aes(label = count), data = venn_region(data), size=5) +
	theme_void(base_size = 20)
fn = paste0(work,"/fig3_venn.png")
png(file =  fn, width = 1600, height = 900, units = 'px')
print(p)
dev.off()

# description -------------------------------------------------------------

p <- variants %>%
	mutate(risk= case_when(Odds.ratio<1~"OR decreased",TRUE~"OR increased")) %>%
	ggplot( aes(x=Odds.ratio, fill=risk)) +
	geom_histogram( alpha=0.6, position = 'identity') +
	xlab("OR") +
	clean_theme(base_size = 20,legend.title=element_blank() )
fn = paste0(work,"/figs3_or.png")
png(file =  fn, width = 1600, height = 900, units = 'px')
print(p)
dev.off()

p <- variants %>%
	mutate(risk= case_when(Odds.ratio<1~"OR decreased",TRUE~"OR increased")) %>%
	ggplot( aes(x=Case.MAF, fill=risk)) +
	geom_histogram( alpha=0.6, position = 'identity', bins=10) +
	xlab("cases MAF") +
	clean_theme(base_size = 20,legend.title=element_blank() )
fn = paste0(work,"/figs2_distr.png")
png(file =  fn, width = 1600, height = 900, units = 'px')
print(p)
dev.off()

p <- variants %>%
	mutate(risk= case_when(Odds.ratio<1~"OR decreased",TRUE~"OR increased")) %>%
	dplyr::select(risk,Consequence) %>%
	group_by(risk,Consequence) %>%
	mutate(count = n()) %>%
	ggplot(aes(x = risk, y = count, fill = factor(Consequence))) +
	geom_bar(stat = "identity", position = position_dodge()) +
	labs(x="MAF cases") +
	clean_theme(base_size = 20,legend.title=element_blank() )

for (i in seq(1,14,by=2)) {
  words <- names(table(variants$Consequence))[i:(i+1)] #14
  p <- variants %>%
  	mutate(risk= case_when(Odds.ratio<1~"OR decreased",TRUE~"OR increased")) %>%
  	filter(Consequence %in% words) %>%
  	ggplot( aes(y=Odds.ratio, x=Case.MAF, colour=trait, label=Gene) ) +
  	scale_colour_manual(values=c("red", "grey")) +
  	geom_point(size=8) +
  	geom_label_repel(aes(color=trait), size=7) +
  	labs(x="MAF cases",y="OR") +
  	geom_hline(yintercept = 1, linetype="dotted", color="black") +
  	geom_vline(xintercept = .25, linetype="dotted", color="black") +
  	facet_wrap(vars(Consequence)) +
  	clean_theme(base_size = 20,legend.title=element_blank(), panel.border = element_rect(color = "black", fill = NA, size = .5) )

  fn = paste0(work,"/figs_combined",i,".png")
  png(file =  fn, width = 1600, height = 900, units = 'px')
  print(p)
  dev.off()
}

# genes -------------------------------------------------------------------

cand <- variants %>%
	dplyr::select(Gene) %>%
	distinct(Gene) %>%
	left_join( BMgene %>%
	  filter(chromosome_name %in% c(1:22,"X","Y")) %>%
		group_by(chromosome_name,external_gene_name) %>%
		filter(start_position == min(start_position)) %>%
		filter(end_position == max(end_position)), by=c("Gene"="external_gene_name") ) %>%
	filter(complete.cases(.)) %>%
	rename(from=start_position,to=end_position,chr=chromosome_name) 

for (i in c(1)) { #115
  	chr = as.numeric(cand[i,"chr"])
  	from = as.numeric(cand[i,"from"]-50000)
  	to = as.numeric(cand[i,"to"]+50000)
  	gene = cand[i,"Gene"]
  	if (from==to) {
  		from=from-25000
  		to=to+25000
  	}
  	tryCatch(
  	  expr= spot_on(chr,from,to,gene),
      error = function(e){ print(e) }
	)
  cat(i," ")
}
