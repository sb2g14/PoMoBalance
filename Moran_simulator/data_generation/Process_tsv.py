# Script to convert RB output in .tsv format to PoMo input for inference
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

def rate_matrix_Rui(N, mu, phi, beta, B):
    num_states = 2 * (3 * N - 1)
    m = np.zeros((num_states, num_states))
    m_test = np.zeros((num_states, num_states))
    # diagonal = np.zeros(num_states)
    m[0][4] = mu[0]
    m[1][N + 2] = mu[1]

    # AG
    m[0][N + 3] = mu[2]
    m[2][2 * N + 1] = mu[3]

    # AT
    m[0][2 * N + 2] = mu[4]
    m[3][3 * N] = mu[5]

    # CG
    m[1][3 * N + 1] = mu[6]
    m[2][4 * N - 1] = mu[7]

    # CT
    m[1][4 * N] = mu[8]
    m[3][5 * N - 2] = mu[9]

    # GT
    m[2][5 * N - 1] = mu[10]
    m[3][6 * N - 3] = mu[11]

    rN = 1.0 / (N)

    # fixations
    # AC
    m[4][0] = (N - 1.0) * (phi[0]) * rN
    m[N + 2][1] = (N - 1.0) * (phi[1]) * rN

    # AG
    m[N + 3][0] = (N - 1.0) * (phi[0]) * rN
    m[2 * N + 1][2] = (N - 1.0) * (phi[2]) * rN

    # AT
    m[2 * N + 2][0] = (N - 1.0) * (phi[0]) * rN
    m[3 * N][3] = (N - 1.0) * (phi[3]) * rN

    # CG
    m[3 * N + 1][1] = (N - 1.0) * (phi[1]) * rN
    m[4 * N - 1][2] = (N - 1.0) * (phi[2]) * rN

    # CT
    m[4 * N][1] = (N - 1.0) * (phi[1]) * rN
    m[5 * N - 2][3] = (N - 1.0) * (phi[3]) * rN

    # GT
    m[5 * N - 1][2] = (N - 1.0) * (phi[2]) * rN
    m[6 * N - 3][3] = (N - 1.0) * (phi[3]) * rN


    #fixations and mutations if N = 2
    if N > 2:

        # frequency shifts from singletons
        # AC
        m[4][5]     = (N-1.0) * (phi[1]) * pow(beta[0], 0.5 * (abs(N-1-B[0])-abs(N-2-B[0])+1.0)) * rN
        m[N+2][N+1]   = (N-1.0) * (phi[0]) * pow(beta[0], 0.5 * (abs(  1-B[0])-abs(  2-B[0])+1.0)) * rN

        m_test[4][5] += 1
        m_test[N + 2][N + 1] += 1
        # AG
        m[N+3][N+4]   = (N-1.0) * (phi[2]) * pow(beta[1], 0.5 * (abs(N-1-B[1])-abs(N-2-B[1])+1.0)) * rN
        m[2 * N+1][2 * N]   = (N-1.0) * (phi[0]) * pow(beta[1], 0.5 * (abs(  1-B[1])-abs(  2-B[1])+1.0)) * rN

        m_test[N+3][N+4] += 1
        m_test[2 * N+1][2 * N] += 1
        # AT
        m[2 * N+2][2 * N+3] = (N-1.0) * (phi[3]) * pow(beta[2], 0.5 * (abs(N-1-B[2])-abs(N-2-B[2])+1.0)) * rN
        m[3 * N][3 * N-1] = (N-1.0) * (phi[0]) * pow(beta[2], 0.5 * (abs(  1-B[2])-abs(  2-B[2])+1.0)) * rN

        m_test[2 * N+2][2 * N+3] += 1
        m_test[3 * N][3 * N-1] += 1
        # CG
        m[3 * N+1][3 * N+2] = (N-1.0) * (phi[2]) * pow(beta[3], 0.5 * (abs(N-1-B[3])-abs(N-2-B[3])+1.0)) * rN
        m[4 * N-1][4 * N-2] = (N-1.0) * (phi[1]) * pow(beta[3], 0.5 * (abs(  1-B[3])-abs(  2-B[3])+1.0)) * rN

        m_test[3 * N+1][3 * N+2] += 1
        m_test[4 * N-1][4 * N-2] += 1
        # CT
        m[4 * N][4 * N+1] = (N-1.0) * (phi[3]) * pow(beta[4], 0.5 * (abs(N-1-B[4])-abs(N-2-B[4])+1.0)) * rN
        m[5 * N-2][5 * N-3] = (N-1.0) * (phi[1]) * pow(beta[4], 0.5 * (abs(  1-B[4])-abs(  2-B[4])+1.0)) * rN

        m_test[4 * N][4 * N+1] += 1
        m_test[5 * N-2][5 * N-3] += 1
        # GT
        m[5 * N-1][5 * N]   = (N-1.0) * (phi[3]) * pow(beta[5], 0.5 * (abs(N-1-B[5])-abs(N-2-B[5])+1.0)) * rN
        m[6 * N-3][6 * N-4] = (N-1.0) * (phi[2]) * pow(beta[5], 0.5 * (abs(  1-B[5])-abs(  2-B[5])+1.0)) * rN

        m_test[5 * N-1][5 * N] += 1
        m_test[6 * N-3][6 * N-4] += 1

        # frequency shifts for all the other polymorphic states
        if N > 3:

            # polymorphic states are populated in two fronts, thus the need for the middle frequency
            S = round(N / 2)+1

            for n in range(2, S):

                # populates the first half of the polymorphic edges
                # AC
                m[n+3][n+4]      = n * (N-n) * (phi[1]) * pow(beta[0], 0.5 * (abs(N-n-B[0])-abs(N-n-1-B[0])+1.0)) * rN
                m[n+3][n+2]      = n * (N-n) * (phi[0]) * pow(beta[0], 0.5 * (abs(N-n-B[0])-abs(N-n+1-B[0])+1.0)) * rN
                m_test[n+3][n+4] += 1
                m_test[n+3][n+2] += 1
                # AG
                m[N+n+2][N+n+3]    = n * (N-n) * (phi[2]) * pow(beta[1], 0.5 * (abs(N-n-B[1])-abs(N-n-1-B[1])+1.0)) * rN
                m[N+n+2][N+n+1]    = n * (N-n) * (phi[0]) * pow(beta[1], 0.5 * (abs(N-n-B[1])-abs(N-n+1-B[1])+1.0)) * rN
                m_test[N+n+2][N+n+3] += 1
                m_test[N+n+2][N+n+1] += 1
                # AT
                m[2 * N+n+1][2 * N+n+2]  = n * (N-n) * (phi[3]) * pow(beta[2], 0.5 * (abs(N-n-B[2])-abs(N-n-1-B[2])+1.0)) * rN
                m[2 * N+n+1][2 * N+n]    = n * (N-n) * (phi[0]) * pow(beta[2], 0.5 * (abs(N-n-B[2])-abs(N-n+1-B[2])+1.0)) * rN
                m_test[2 * N+n+1][2 * N+n+2] += 1
                m_test[2 * N+n+1][2 * N+n] += 1
                # CG
                m[3 * N+n][3 * N+n+1]  = n * (N-n) * (phi[2]) * pow(beta[3], 0.5 * (abs(N-n-B[3])-abs(N-n-1-B[3])+1.0)) * rN
                m[3 * N+n][3 * N+n-1]  = n * (N-n) * (phi[1]) * pow(beta[3], 0.5 * (abs(N-n-B[3])-abs(N-n+1-B[3])+1.0)) * rN
                m_test[3 * N+n][3 * N+n+1] += 1
                m_test[3 * N+n][3 * N+n-1] += 1
                # CT
                m[4 * N+n-1][4 * N+n]    = n * (N-n) * (phi[3]) * pow(beta[4], 0.5 * (abs(N-n-B[4])-abs(N-n-1-B[4])+1.0)) * rN
                m[4 * N+n-1][4 * N+n-2]  = n * (N-n) * (phi[1]) * pow(beta[4], 0.5 * (abs(N-n-B[4])-abs(N-n+1-B[4])+1.0)) * rN
                m_test[4 * N+n-1][4 * N+n] += 1
                m_test[4 * N+n-1][4 * N+n-2] += 1
                # GT
                m[5 * N+n-2][5 * N+n-1]  = n * (N-n) * (phi[3]) * pow(beta[5], 0.5 * (abs(N-n-B[5])-abs(N-n-1-B[5])+1.0)) * rN
                m[5 * N+n-2][5 * N+n-3]  = n * (N-n) * (phi[2]) * pow(beta[5], 0.5 * (abs(N-n-B[5])-abs(N-n+1-B[5])+1.0)) * rN
                m_test[5 * N+n-2][5 * N+n-1] += 1
                m_test[5 * N+n-2][5 * N+n-3] += 1
                # populates the second half of the polymorphic edges
                # AC
                m[N-n+3][N-n+2]     = (N-n) * n * (phi[0]) * pow(beta[0], 0.5 * (abs(n-B[0])-abs(n+1-B[0])+1.0)) * rN
                m[N-n+3][N-n+4]     = (N-n) * n * (phi[1]) * pow(beta[0], 0.5 * (abs(n-B[0])-abs(n-1-B[0])+1.0)) * rN

                # AG
                m[2 * N-n+2][2 * N-n+1] = (N-n) * n * (phi[0]) * pow(beta[1], 0.5 * (abs(n-B[1])-abs(n+1-B[1])+1.0)) * rN
                m[2 * N-n+2][2 * N-n+3] = (N-n) * n * (phi[2]) * pow(beta[1], 0.5 * (abs(n-B[1])-abs(n-1-B[1])+1.0)) * rN

                # AT
                m[3 * N-n+1][3 * N-n]   = (N-n) * n * (phi[0]) * pow(beta[2], 0.5 * (abs(n-B[2])-abs(n+1-B[2])+1.0)) * rN
                m[3 * N-n+1][3 * N-n+2] = (N-n) * n * (phi[3]) * pow(beta[2], 0.5 * (abs(n-B[2])-abs(n-1-B[2])+1.0)) * rN

                # CG
                m[4 * N-n][4 * N-n-1]   = (N-n) * n * (phi[1]) * pow(beta[3], 0.5 * (abs(n-B[3])-abs(n+1-B[3])+1.0)) * rN
                m[4 * N-n][4 * N-n+1]   = (N-n) * n * (phi[2]) * pow(beta[3], 0.5 * (abs(n-B[3])-abs(n-1-B[3])+1.0)) * rN

                # CT
                m[5 * N-n-1][5 * N-n-2] = (N-n) * n * (phi[1]) * pow(beta[4], 0.5 * (abs(n-B[4])-abs(n+1-B[4])+1.0)) * rN
                m[5 * N-n-1][5 * N-n]   = (N-n) * n * (phi[3]) * pow(beta[4], 0.5 * (abs(n-B[4])-abs(n-1-B[4])+1.0)) * rN

                # GT
                m[6 * N-n-2][6 * N-n-3] = (N-n) * n * (phi[2]) * pow(beta[5], 0.5 * (abs(n-B[5])-abs(n+1-B[5])+1.0)) * rN
                m[6 * N-n-2][6 * N-n-1] = (N-n) * n * (phi[3]) * pow(beta[5], 0.5 * (abs(n-B[5])-abs(n-1-B[5])+1.0)) * rN

    for i in range(len(m[0])):
        m[i, i] = -sum(m[i, :])
    return m

if __name__ == '__main__':

    N = 20
    case_list = ["neutral", "sel", "bal", "sel+bal"]
    num = 0
    case = case_list[num]
    n_sam = 100000
    # Theoretical values
    pi = [0.25, 0.25, 0.25, 0.25]
    rho = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    if num == 0 or num == 2:
        sigma = 0
    elif num == 1 or num == 3:
        sigma = 0.1

    phi = [1., 1.+sigma, 1.+sigma, 1.]
    if num == 0 or num == 1:
        beta_val = 1.0
    elif num == 2 or num == 3:
        beta_val = 2.0

    beta = np.multiply([1, 1, 1, 1, 1, 1], beta_val)
    B_val = 2
    B = np.multiply([1, 1, 1, 1, 1, 1], B_val)

    path_to_data = '/Users/Documents/Tools/PoMoBalance/Moran_simulator/data_generation/{0}/output/Validation_Sim_0/'.format(case, n_sam)
    poMoStates = np.genfromtxt(path_to_data+'sequences.tsv', dtype=int)
    n_pop = len(poMoStates[:, 0])
    n_states = 4 + 6 * (N - 1)
    n_nucl = len(poMoStates[0, :]) - 1
    p_read = np.zeros((n_pop, n_nucl))
    p = np.flip(p_read, axis=0)
    print(p[p > n_states - 1])
    # Output file
    f = open(path_to_data + 'sequences_{0}.txt'.format(case), 'w')
    for sel in range(n_pop):
        p[sel, :] = poMoStates[sel, 1:].astype(int)
        counts = np.zeros(n_states)
        for i in range(0, n_states):
            counts[i] = list(p[sel, :]).count(i)

        freq_sim = counts/n_nucl

        # Non-reversible
        mu = [pi[1] * rho[0], pi[0] * rho[0], pi[2] * rho[1], pi[0] * rho[1], pi[3] * rho[2], pi[0] * rho[2],
              pi[2] * rho[3], pi[1] * rho[3], pi[3] * rho[4], pi[1] * rho[4], pi[3] * rho[5], pi[2] * rho[5]]

        rm = rate_matrix_Rui(N, mu, phi, beta, B)

        t = 1e6
        Pnr = expm(rm * t)
        freq = Pnr[0, :]
        x = np.arange(0, N+1)
        y_sim = np.zeros((6, N+1))
        y = np.zeros((6, N+1))


        comb = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        for i in range(6):
            # From simulations
            y_sim[i, 0] = freq_sim[comb[i][0]]
            y_sim[i, 1:N] = freq_sim[4 + (N-1)*i:4 + (N - 1)*(i+1)]
            y_sim[i, N] = freq_sim[comb[i][1]]
            # Theory
            y[i, 0] = freq[comb[i][0]]
            y[i, 1:N] = freq[4 + (N-1)*i:4 + (N - 1)*(i+1)]
            y[i, N] = freq[comb[i][1]]
        markers = ['*', 'd', 'X', 'o', 'h', '<', 'x', '1', 's', '+', 'P', '*', 'd', '.', 'o']
        colors = ['teal', 'orangered', 'forestgreen', 'purple', 'gold', 'magenta', 'cyan', 'coral', 'b']
        font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 15}
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif', size=25)
        fig = plt.figure(figsize=(28, 18), facecolor='w', edgecolor='k')
        lab = ["AC", "AG", "AT", "CG", "CT", "GT"]
        path_to_figures = path_to_data
        for i in range(0, 6):
            fig.add_subplot(2, 3, i + 1)
            plt.plot(x, y_sim[i, :], linestyle="None", markerfacecolor='none', markeredgewidth=3, color=colors[1],
                     marker=markers[1],
                     markersize=30, alpha=0.6, label="Data")
            plt.plot(x, y[i, :], linestyle="None", markerfacecolor='none', markeredgewidth=3, color=colors[0],
                     marker=markers[0],
                     markersize=30, alpha=0.6, label="Theory")
            plt.xlabel('N')
            plt.ylabel('St freq {0}'.format(lab[i]))
        # Set axes:
        axes = plt.gca()
        plt.legend()
        plt.tight_layout()
        # plt.savefig(path_to_figures + "Stationary_frequencies_{0}.pdf".format(sel+1))
        plt.show()
        f.write("pop_{0}".format(sel))
        for state in p[sel, :]:
            if (state).is_integer() is False:
                print("Whoops {0} is not integer".format(state))
            f.write(" {0}".format(int(state)))
        f.write("\n")
    f.close()
