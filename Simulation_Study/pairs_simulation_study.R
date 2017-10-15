library('gtools');
sim_Markov_data <- function(CPs, transition_matrices_list, ndata, numberOfStates) {
# simulates data with specified changepoints and specified matrix structure
	simdata <- rep(0, ndata);
	simdata[1] <- sample(numberOfStates, 1);
	matrixon <- 1;
	for(i in 2:ndata - 1) {
		if(i %in% CPs) {
			matrixon <- matrixon + 1
		}
		simdata[i + 1] <- sample(1:numberOfStates, 1, prob = transition_matrices_list[[matrixon]][simdata[i],]);
	}
	simdata - 1;
}


generate_transition_vector_list <- function(rho, markov_alpha, number_of_states) {
	#for each process, generate a set of multinomial probability vectors
	#each element of the list represents a process, a matrix of (number of regimes) x (number of states)
	row_list <- list()
	for (i in 1:number_of_states) {
		number_of_regimes <- rgeom(1, rho) + 1
		row_list[[i]] <- rdirichlet(number_of_regimes, rep(markov_alpha, number_of_states))
	}
	row_list
}

# for each process, create a transition matrix between regimes
generate_regime_transition_matrices_list = function(transition_vector_list, regime_markov_alpha, number_of_states) {
	regime_transition_matrices_list = list()
	m = number_of_states
	for (i in 1:m) {
		number_of_regimes <- dim(transition_vector_list[[i]])[1]
		regime_transition_matrices_list[[i]] <- rdirichlet(number_of_regimes, rep(regime_markov_alpha, number_of_regimes))
	}
	regime_transition_matrices_list
}

simulate_regime_vectors <- function(cps, regime_markov_alpha, transition_vector_list, regime_transition_matrices_list, number_of_states, wildcard_first_regime) {
	m <- number_of_states
	regime_vectors <- list()
	# each element of regime vectors is a vector saying which regimes each process adopts at each changepoint
	if (wildcard_first_regime) {
		number_of_regimes = sapply(TVL, function(x) {dim(x)[1]})
		regimes = numeric(number_of_states)
		for (i in 1:number_of_states) { regimes[i] = sample(number_of_regimes[i], 1) }
	} else {
		regimes <- rep(1, m)
	}
	regime_vectors[[1]] = regimes
	if (length(cps) > 0) {
		for (cp_idx in 1:length(cps)) {
			# create the TPM that begins operation at changepoint cp_idx
			old_regimes <- regimes
			# generate the next regime for each process
			regimes <- numeric(m)
			for (proc in 1:m) {
				number_of_regimes <- dim(transition_vector_list[[proc]])[1]
				reg_prob <- regime_transition_matrices_list[[proc]][old_regimes[proc],]
				regimes[proc] <- sample(number_of_regimes, 1, prob = reg_prob)
			}
			regime_vectors[[cp_idx + 1]] <- regimes
		}
	}
	regime_vectors
}

generate_TPMs_list <- function(transition_vector_list, regime_vectors, number_of_states) {
	# generate the set of TPMs for the Markov chain based on the transition vector list and the regime transition matrices list
	m <- number_of_states
	transition_matrices_list <- list()
	# generate the first TPM
	regimes = regime_vectors[[1]]
	tpm <- matrix(0, nrow = m, ncol = m)
	for (proc in 1:m) {
		tpm[proc,] <- transition_vector_list[[proc]][regimes[proc],] # using 1 here because this is the first TPM
	}
	transition_matrices_list[[1]] <- tpm
	number_of_transitions <- length(regime_vectors) - 1
	if (number_of_transitions > 0) {
		for (cp_idx in 1:number_of_transitions) {
			# create the TPM that begins operation at changepoint cp_idx
			regimes <- regime_vectors[[cp_idx + 1]]
			tpm <- matrix(0, nrow = m, ncol = m)
			for (proc in 1:m) {
				tpm[proc,] <- transition_vector_list[[proc]][regimes[proc],]
			}
			transition_matrices_list[[cp_idx + 1]] <- tpm
		}
	}
	transition_matrices_list
}

generate_changepoints <- function(n, p) {
	which(rbinom(n - 1, 1, p) == 1)
}

calculate_number_of_observations = function(process, data, indices, cps_with_0_end) {
	number_of_observations = 0
	for (index in indices) {
		# miraculously this works even with 0's e.g. > (1:4)[0:2] outputs [1] 1 2
		left_position = cps_with_0_end[index]
		right_position = cps_with_0_end[index + 1] - 1
		segment_data = data[left_position:right_position]
		number_of_observations = number_of_observations + sum(segment_data == process)
	}
	return(number_of_observations)
}

calculate_min_prop_similarity = function(first_RVs, second_RVs, first_data, second_data, first_cps_with_0_end, second_cps_with_0_end) {
	# work out the amount of observations spent in common regimes using the regime vectors for both traces
	first_common_regime_obs = 0
	second_common_regime_obs = 0
	number_of_states = length(first_RVs[[1]])
	for (proc in 1:number_of_states) {
		first_observed_regimes = unique(sapply(first_RVs, function(x) { x[proc] }))
		second_observed_regimes = unique(sapply(second_RVs, function(x) { x[proc] }))
		common_regimes = intersect(first_observed_regimes, second_observed_regimes)
		for (reg in common_regimes) {
			# work out which CPs represent this regime for this process
			first_regimes = sapply(first_RVs, function(x) { x[proc] })
			first_reg_indices = which(first_regimes == reg)
			first_number_of_observations = calculate_number_of_observations(proc, first_data, first_reg_indices, first_cps_with_0_end)
			first_common_regime_obs = first_common_regime_obs + first_number_of_observations
			second_regimes = sapply(second_RVs, function(x) { x[proc] })
			second_reg_indices = which(second_regimes == reg)
			second_number_of_observations = calculate_number_of_observations(proc, second_data, second_reg_indices, second_cps_with_0_end)
			second_common_regime_obs = second_common_regime_obs + second_number_of_observations
		}
	}
	similarity = min(first_common_regime_obs, second_common_regime_obs) / (length(first_data) - 1)
	similarity
}

calculate_similarity = function(first_RVs, second_RVs, first_data, second_data, first_cps_with_0_end, second_cps_with_0_end) {
	common_regime_obs = 0
	number_of_states = length(first_RVs[[1]])
	for (proc in 1:number_of_states) {
		first_observed_regimes = unique(sapply(first_RVs, function(x) { x[proc] }))
		second_observed_regimes = unique(sapply(second_RVs, function(x) { x[proc] }))
		common_regimes = intersect(first_observed_regimes, second_observed_regimes)
		for (reg in common_regimes) {
			# work out which CPs represent this regime for this process
			first_regimes = sapply(first_RVs, function(x) { x[proc] })
			first_reg_indices = which(first_regimes == reg)
			first_number_of_observations = calculate_number_of_observations(proc, first_data, first_reg_indices, first_cps_with_0_end)
			common_regime_obs = common_regime_obs + first_number_of_observations
			second_regimes = sapply(second_RVs, function(x) { x[proc] })
			second_reg_indices = which(second_regimes == reg)
			second_number_of_observations = calculate_number_of_observations(proc, second_data, second_reg_indices, second_cps_with_0_end)
			common_regime_obs = common_regime_obs + second_number_of_observations
		}
	}
	similarity = common_regime_obs / (length(first_data) - 1 + length(second_data) - 1)
	similarity
}

calculate_effective_number_of_CPs <- function(regime_vectors) {
	if (length(regime_vectors) == 1) {
		return(0)
	}
	prev_regime_vec <- regime_vectors[[1]]
	effective_number_of_CPs <- 0
	k <- length(regime_vectors)
	for (i in 2:k) {
		new_regime_vec <- regime_vectors[[i]]
		if (any(new_regime_vec != prev_regime_vec)) {
			effective_number_of_CPs <- effective_number_of_CPs + 1
		}
		prev_regime_vec <- new_regime_vec
	}
	return(effective_number_of_CPs)
}


args <- commandArgs(TRUE)
basename <- args[1]
set.seed(as.integer(basename))
num_states <- 4
ndata <- 10000
markov_alpha <- 1
dirichlet_alpha <- 1
rho <- 0.2
p <- 5 / ndata
# concatenate two traces drawn from the same regime transition vector list (which defines regimes) and the same regime_transition_matrices_list (which does probs transition 'twixt regimes), and writes the trace to file as well as writing the actual min proportion similarity of the two traces to file.
# generate the regimes and regime transition probabilities for each process
TVL = generate_transition_vector_list(rho, markov_alpha, num_states)
RTM = generate_regime_transition_matrices_list(TVL, dirichlet_alpha, num_states)

# generate the change points and TPMs list for the first trace
first_cps = generate_changepoints(ndata, p)
first_cps_with_0_end = c(0, first_cps, ndata)
first_RVs = simulate_regime_vectors(first_cps, regime_markov_alpha, TVL, RTM, num_states, FALSE)
first_TPMs = generate_TPMs_list(TVL, first_RVs, num_states)

# generate the change points and TPMs list for the second trace
second_cps = generate_changepoints(ndata, p)
second_cps_with_0_end = c(0, second_cps, ndata)
second_RVs = simulate_regime_vectors(second_cps, regime_markov_alpha, TVL, RTM, num_states, TRUE)
second_TPMs = generate_TPMs_list(TVL, second_RVs, num_states)

# generate the two traces
first_data = sim_Markov_data(first_cps, first_TPMs, ndata, num_states)
second_data = sim_Markov_data(second_cps, second_TPMs, ndata, num_states)

min_prop_similarity = calculate_min_prop_similarity(first_RVs, second_RVs, first_data + 1, second_data + 1, first_cps_with_0_end, second_cps_with_0_end)
similarity = calculate_similarity(first_RVs, second_RVs, first_data + 1, second_data + 1, first_cps_with_0_end, second_cps_with_0_end)
x <- c(first_data, second_data)
separators = c(ndata, 2 * ndata)
for(state in 1:num_states - 1) {
	idx <- which(x == state)
	idx <- idx[!is.element(idx, separators)]
	trans = x[idx + 1]
	name1 <- paste(basename, "_", "data_state_", state, ".txt", sep = "")
	name2 <- paste(basename, "_", "times_state_", state, ".txt", sep = "")
	write(trans, file = name1, ncolumns = length(trans));
	write(idx, file = name2, ncolumns = length(idx))
}
write(x, file = paste(basename, "_data.txt", sep = ""), ncol = 1)
write(min_prop_similarity, paste(basename, "_actual_min_proportion_similarity.txt", sep = ""))
write(similarity, paste(basename, "_actual_similarity.txt", sep = ""))
write(length(first_cps) + length(second_cps), file = paste(basename, "_actual_number_of_CPs.txt", sep = ""))
write(calculate_effective_number_of_CPs(first_RVs) + calculate_effective_number_of_CPs(second_RVs), file = paste(basename, "_actual_effective_number_of_CPs.txt", sep = ""))













