# Create logistic models to examine the effects of case-control status on each ICD10 code
files <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 93, 98)
controls <- c('NoCNV_Control', 'LargeRare_Control', 'NEJM_Control')

output <- data.frame(GO_id=rep(NA, 40000), control=rep(NA, 40000), case_w=rep(NA, 40000), case_wo=rep(NA, 40000), control_w=rep(NA, 40000), control_wo=rep(NA, 40000),
	cc_est=rep(NA, 40000), cc_err=rep(NA, 40000), cc_p=rep(NA, 40000),
	sex_est=rep(NA, 40000), sex_err=rep(NA, 40000), sex_p=rep(NA, 40000), stringsAsFactors=FALSE)

i <- 1
for (suff in files) {
	print(suff)
	file=paste('analysis_files/12_go_matrix/', as.character(suff), '.csv', sep='')
	df <- read.csv(file)

	for (c in controls) {
		df2 <- df[(df$Case_Control==c | df$Case_Control=='Case'), ]
		# Make inputs numeric
		df2$case_input <- 0
		df2[df2$Case_Control=='Case', 'case_input'] <- 1

		df2$Case_Control <- 1
		# Remove invariant GO terms
		df2 <- df2[, colSums(df2)>0]
		# Make inputs into factors
		df2$Case_Control <- as.factor(df2$case_input)
		df2$Sex <- as.factor(df2$Sex)

		terms <- colnames(df2)
		terms <- terms[! terms %in% c('Case_Control', 'Sex', 'Sample', 'case_input')]
		
		for (t in terms) {
			df2$output <- as.factor(df2[, t])
			model <- glm(output ~ Case_Control + Sex, family='binomial', data=df2)
			cc_est <- unlist(summary(model)$coefficients[2, 1])
			cc_err <- unlist(summary(model)$coefficients[2, 2])
			cc_p <- unlist(summary(model)$coefficients[2, 4])

			sex_est <- unlist(summary(model)$coefficients[3, 1])
			sex_err <- unlist(summary(model)$coefficients[3, 2])
			sex_p <- unlist(summary(model)$coefficients[3, 4])

			case_w <- length(rownames(df2[((df2[, t]==1) & (df2$case_input==1)), ]))
			case_wo <- length(rownames(df2[((df2[, t]==0) & (df2$case_input==1)), ]))
			cont_w <- length(rownames(df2[((df2[, t]==1) & (df2$case_input==0)), ]))
			cont_wo <- length(rownames(df2[((df2[, t]==0) & (df2$case_input==0)), ]))

			output[i, ] <- list(t, c, case_w, case_wo, cont_w, cont_wo, cc_est, cc_err, cc_p, sex_est, sex_err, sex_p)
			i <- i+1
		}
	}
}

# Remove any empty rows
output <- output[!is.na(output$GO_id), ]

# Save to file
write.csv(output, 'result_tables/13_go_logistic_models.csv', row.names=F)
