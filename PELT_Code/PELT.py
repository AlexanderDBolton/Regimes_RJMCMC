import numpy as np

class probability_model:
    def __init__(self, x, end_time, beta, likelihood_model):
        self.x = x
        self.end_time = end_time
        self.beta = beta
        self.likelihood_model = likelihood_model
    
    def cost(self, t0, t1):
        if self.likelihood_model == "gaussian":
            n = t1 - t0
            # check if there is only one observation in the interval
            if n == 1:
                # the likelihood for the observation will be 1, so log(likelihood) = 0
                return 0
            
            mu_hat = np.mean(self.x[t0:t1])
            sigma_squared_hat = np.std(self.x[t0:t1])
    
            if sigma_squared_hat < 0:
                raise ValueError("Gaussian variance estimated to be below 0")
    
            # check if the variance is 0
            if sigma_squared_hat < 0.0000001:
                # the likelihood for the observation will be 1, so log(likelihood) = 0
                return 0
    
            maximum_log_likelihood = -(np.log(2*np.pi) + np.log(sigma_squared_hat) + 1) * n / 2
            return -maximum_log_likelihood * 2.0
        else:
            raise ValueError("Model not implemented")


def calculate_F_tau_star_and_tau_1(F, R_tau_star, beta, tau_1, tau_star, pm):
    tau_1[0] = R_tau_star[0]
    min_cost = F[tau_1[0]] + pm.cost(tau_1[0], tau_star) + beta
    for tau in R_tau_star[1:]:
        cost = F[tau] + pm.cost(tau, tau_star) + beta
        if cost < min_cost:
            min_cost = cost
            tau_1[0] = tau
    F[tau_star] = min_cost


def set_R_tau_star_plus_1(R_tau_star_plus_1, R_tau_star, F, tau_star, K, pm):
    for tau in R_tau_star:
        if F[tau] + pm.cost(tau, tau_star) + K <= F[tau_star]:
            R_tau_star_plus_1.append(tau)
    R_tau_star_plus_1.append(tau_star)


#def print_asterisk_if_needed(x, n, number_already_printed, num_asterisks):
#    while x > n * number_already_printed[0] / num_asterisks:
#        print("*")
#        number_already_printed[0] += 1


def PELT_cp_n(x, end_time, likelihood_model = "gaussian", K = 0, penalty = "BIC"):
    if penalty == "BIC":
        beta = 2.0 * np.log(end_time + 1)
    elif penalty == "AIC":
        beta = 2.0 * 2.0
    
    pm = probability_model(x, end_time, beta, likelihood_model)

    F = np.zeros(end_time + 2)
    R = [[] for _ in range(end_time + 3)]
    cp = [[] for _ in range(end_time + 2)]

    F[0] = -beta
    R[1].append(0)
    number_of_asterisks_already_printed = [0]

    for tau_star in range(1, end_time + 2):
        tau_1 = [0]
        calculate_F_tau_star_and_tau_1(F, R[tau_star], beta, tau_1, tau_star, pm)
        cp[tau_star] = cp[tau_1[0]][:]
        cp[tau_star].append(tau_1[0])
        set_R_tau_star_plus_1(R[tau_star + 1], R[tau_star], F, tau_star, K, pm)
        #print_asterisk_if_needed(tau_star, end_time, number_of_asterisks_already_printed, 10)

    return cp