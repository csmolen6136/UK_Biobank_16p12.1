# Create logistic models to examine the effects of case-control status on each ICD10 code
chroms <- head(read.csv('analysis_files/chromosomes.list', header=F, col.names=c('col1'))$col1, -1)
controls <- c('NoCNV_Control', 'LargeRare_Control', 'NEJM_Control')

output <- data.frame(chrom=rep(NA, 34000), gene=rep(NA, 34000), control=rep(NA, 34000), case_n=rep(NA, 34000), control_n=rep(NA, 34000),
	cc_est=rep(NA, 34000), cc_err=rep(NA, 34000), cc_p=rep(NA, 34000),
	sex_est=rep(NA, 34000), sex_err=rep(NA, 34000), sex_p=rep(NA, 34000), stringsAsFactors=FALSE)

i <- 1
for (chr in chroms) {
	print(chr)
	file=paste('analysis_files/4_final_secondhit_tables/4_', chr, '_secondhit_table.csv', sep='')
	df <- read.csv(file)

	for (c in controls) {
		df2 <- df[(df$Case_Control==c | df$Case_Control=='Case'), ]
		# Make inputs numeric
		df2$case_input <- 0
		df2[df2$Case_Control=='Case', 'case_input'] <- 1

		case_n <- sum(df2$case_input)
		control_n <- length(rownames(df2))-case_n

		df2$Case_Control <- 1
		# Remove invariant genes
		df2 <- df2[, colSums(df2)>0]
		# Make inputs into factors
		df2$Case_Control <- as.factor(df2$case_input)
		df2$Sex <- as.factor(df2$Sex)

		genes <- colnames(df2)
		genes <- genes[! genes %in% c('Case_Control', 'Sex', 'Sample', 'case_input')]
		
		for (g in genes) {
			df2$output <- as.factor(df2[, g])
			model <- glm(output ~ Case_Control + Sex, family='binomial', data=df2)
			cc_est <- unlist(summary(model)$coefficients[2, 1])
			cc_err <- unlist(summary(model)$coefficients[2, 2])
			cc_p <- unlist(summary(model)$coefficients[2, 4])

			sex_est <- unlist(summary(model)$coefficients[3, 1])
			sex_err <- unlist(summary(model)$coefficients[3, 2])
			sex_p <- unlist(summary(model)$coefficients[3, 4])

			output[i, ] <- list(chr, g, c, case_n, control_n, cc_est, cc_err, cc_p, sex_est, sex_err, sex_p)
			i <- i+1
		}
	}
}

# Remove any empty rows
output <- output[!is.na(output$chr), ]

# Save to file
write.csv(output, 'result_tables/5_logistic_models.csv', row.names=F)
