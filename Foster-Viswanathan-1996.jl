# Compute the solution for the difference equations of the market equilibrium in
# Foster & Viswanathan (1996, JF)
#
# by Wataru Nozawa
# June 21, 2019

# Pkg.add("Roots")
using Roots

# Pkg.add("Plots")
using Plots
gr()

#--- Parameters

N = 4  # number of auctions
M = 3  # number of informed

sv = 1  # variance of v, fundamental value of the asset

S_0 = 1  # conditional variance of [estimate of v given all signals]

chi = 0.2

theta = 3

L_0 = S_0/theta^2 + chi*(M-1)/M  # variance of signal
O_0 = S_0/theta^2 - chi/M  # covariance of signal


#--- boundary condition for a_N

a_N = 0
ps_N = 0
m_N = 0
d_N = 0



# # variables for storing solution
# a = zeros(1,N)
# ps = zeros(1,N)
# m = zeros(1,N)
# d = zeros(1,N)
# l = zeros(1,N)
# g = zeros(1,N)
# b = zeros(1,N)
# d = zeros(1,N)
# eta = zeros(1,N)
# ph = zeros(1,N)
# L = zeros(1,N)
# O = zeros(1,N)
# S = zeros(1,N)

#--- Intial guess for L_N

L[N] = S_0/theta^2 + chi*(M-1)/M

O[N] = L[N] - chi
S[N] = (L[N] + (M-1)*O[N])*theta^2/M


# solve

#--- functions for solution

function eq_l_n(l_n, S_n, ps_n, L_n)
    fourth = l_n^4*(theta*sv^4*(M-1)*chi/(M^2*S_n^2))
    third = l_n^3*sv^2*ps_n*L_n/S_n
    second = -l_n^2*(theta*sv^2/M*S_n)*((M+1)*L_n - (M-1)*chi)
    first = -l_n*L_n*ps_n
    constant = (theta/M)*(M*L_n - (M-1)*chi)
    return fourth + third + second + first + constant  # equation (36)
end

get_beta(l_n, S_n) = l_n*theta*sv^2/(M*S_n)  # third equation in Proposition 1

get_phi(O_n, l_n, b_n, L_n) = (O_n + (l_n*b_n/theta)*chi)/(L_n - (l_n*b_n/theta)*(M-1)*chi)  # equation (35)

get_gamma(m_n, l_n, b_n) = (1-2*m_n*l_n)*(1-(l_n*b_n*(M-1)/theta))/(2*l_n*(1-m_n*l_n))  # second equation in Proposition 1
get_eta(phi_n) = theta/M*(1+(M-1)*phi_n)

get_alpha(eta_n, l_n, b_n, a_n, ph_n) = (eta_n - l_n*b_n*(1+(M-1)*ph_n))*b_n + a_n*(1 - (l_n*b_n/theta)*(1+ph_n*(M-1)))^2
get_psi(eta_n, l_n, b_n, a_n, ph_n, gamma_n, psi_n) = (eta_n - l_n*b_n*(1+(M-1)*ph_n))*gamma_n - l_n*gamma_n*b_n + b_n*(1-l_n*b_n*(M-1)/theta) + psi_n*(1-l_n*b_n*(1+(M-1)*ph_n)/theta)*(1-l_n*b_n*(M-1)/theta - l_n*gamma_n)
get_mu(l_n, gamma_n, b_n, m_n) = -l_n*gamma_n^2 + gamma_n*(1-l_n*b_n*(M-1)/theta) + m_n*(1-l_n*b_n*(M-1)/theta - l_n*gamma_n)^2

get_L(L_n, l_n, b_n) = (L_n - (l_n*b_n/theta)*(M-1)*chi)/(1-M*l_n*b_n/theta)
get_O(O_n, l_n, b_n) = (O_n + (l_n*b_n/theta)*chi)/(1-(M*l_n*b_n/theta))

function get_delta(d_n, a_n, l_n, b_n, L_n_1, O_n_1)
    Var = (L_n_1^2 - O_n_1^2)/L_n_1
    Cov = O_n_1*(L_n_1 - O_n_1)/L_n_1
    return d_n + a_n*l_n^2*sv^2/theta^2 + a_n*(l_n^2*b_n^2/theta^2)*(M-1)*(Var + (M-2)*Cov)
end


#--- solve

function solve_system(L_N)
    # variables for storing solution
    a = zeros(1,N)
    ps = zeros(1,N)
    m = zeros(1,N)
    d = zeros(1,N)
    l = zeros(1,N)
    g = zeros(1,N)
    b = zeros(1,N)
    d = zeros(1,N)
    eta = zeros(1,N)
    ph = zeros(1,N)
    L = zeros(1,N)
    O = zeros(1,N)
    S = zeros(1,N)

    L[N] = L_N
    O[N] = L_N - chi
    S[N] = (L_N + (M-1)*O[N])*theta^2/M

    for n in [N:-1:2;]
        # compute lam_n from S_n & a_n
        eq_l(l) = eq_l_n(l, S[n], ps[n], L[n])
        l[n] = find_zero(eq_l, (0, S[n]^(1/2)/sv))

        b[n] = get_beta(l[n], S[n])
        ph[n] = get_phi(O[n], l[n], b[n], L[n])
        g[n] = get_gamma(m[n], l[n], b[n])
        eta[n] = get_eta(ph[n])

        a[n-1] = get_alpha(eta[n], l[n], b[n], a[n], ph[n])
        ps[n-1] = get_psi(eta[n], l[n], b[n], a[n], ph[n], g[n], ps[n])
        m[n-1] = get_mu(l[n], g[n], b[n], m[n])

        L[n-1] = get_L(L[n], l[n], b[n])
        O[n-1] = get_O(O[n], l[n], b[n])
        S[n-1] = (L[n-1] + (M-1)*O[n-1])*theta^2/M
        d[n-1] = get_delta(d[n], a[n], l[n], b[n], L[n-1], O[n-1])
    end
    return L[1] - L_0
end

# solve
L_N_eqm = find_zero(solve_system, (0, 5))
## draw figures

plot(l[2:N]')

plot(S')
