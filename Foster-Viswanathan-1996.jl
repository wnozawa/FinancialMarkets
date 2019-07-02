# Compute the solution for the difference equations of the market equilibrium in
# Foster & Viswanathan (1996, JF)
#
# by Wataru Nozawa
# June 21, 2019
#
# July 2, 2019; Now reproduces Panel A in Figure 1

#--- Packages
# Pkg.add("Roots")
using Roots

# Pkg.add("LaTeXStrings")
using LaTeXStrings

# Pkg.add("Plots")
# using Plots; pgfplots()
using Plots; gr(dpi=400)
# using Plots; plotlyjs()
# using Plotly


#--- Parameters

N = 4  # number of auctions
M = 3  # number of informed

su = 2  # variance of v, fundamental value of the asset

S_0 = 1  # conditional variance of [estimate of v given all signals]

theta = 3

chis = [0.000004 0.2 0.333333 1.0]
Nchi = length(chis)

#--- boundary condition for a_N

a_N = 0
ps_N = 0
m_N = 0
d_N = 0

L_N_eqm = zeros(1,Nchi)
l = zeros(Nchi,N+1)
L = zeros(Nchi,N+1)
S = zeros(Nchi,N+1)
O = zeros(Nchi,N+1)


for idx_chi in 1:Nchi
    chi = chis[idx_chi]
    println(chi)

    L_0 = S_0/theta^2 + chi*(M-1)/M  # variance of signal
    O_0 = S_0/theta^2 - chi/M  # covariance of signal

    #--- functions for solution

    function eq_l_n(l_n, S_n, ps_n, L_n)
        fourth = l_n^4*(theta*su^4*(M-1)*chi/(M^2*S_n^2))
        third = l_n^3*su^2*ps_n*L_n/S_n
        second = -l_n^2*(theta*su^2/(M*S_n))*((M+1)*L_n - (M-1)*chi)
        first = -l_n*L_n*ps_n
        constant = (theta/M)*(M*L_n - (M-1)*chi)
        return fourth + third + second + first + constant  # equation (36)
    end

    get_beta(l_n, S_n) = l_n*theta*su^2/(M*S_n)  # third equation in Proposition 1

    get_phi(O_n, l_n, b_n, L_n) = (O_n + (l_n*b_n/theta)*chi)/(L_n - (l_n*b_n/theta)*(M-1)*chi)
    # equation (35)

    get_gamma(m_n, l_n, b_n) = (1-2*m_n*l_n)*(1-(l_n*b_n*(M-1)/theta))/(2*l_n*(1-m_n*l_n))
    # second equation in Proposition 1
    get_eta(phi_n) = theta/M*(1+(M-1)*phi_n)  # 9th equation in Proposition 1

    get_alpha(eta_n, l_n, b_n, a_n, ph_n) = (eta_n - l_n*b_n*(1+(M-1)*ph_n))*b_n + a_n*(1 - (l_n*b_n/theta)*(1+ph_n*(M-1)))^2  # 4th equation in Proposition 1

    get_psi(eta_n, l_n, b_n, a_n, ph_n, gamma_n, psi_n) = (eta_n - l_n*b_n*(1+(M-1)*ph_n))*gamma_n - l_n*gamma_n*b_n + b_n*(1-l_n*b_n*(M-1)/theta) + psi_n*(1-l_n*b_n*(1+(M-1)*ph_n)/theta)*(1-l_n*b_n*(M-1)/theta - l_n*gamma_n)  # 5th equation in Proposition 1

    get_mu(l_n, gamma_n, b_n, m_n) = -l_n*gamma_n^2 + gamma_n*(1-l_n*b_n*(M-1)/theta) + m_n*(1-l_n*b_n*(M-1)/theta - l_n*gamma_n)^2  # 6th equation in Proposition 1

    get_L(L_n, l_n, b_n) = (L_n - (l_n*b_n/theta)*(M-1)*chi)/(1-M*l_n*b_n/theta)  # equation (33)
    get_O(O_n, l_n, b_n) = (O_n + (l_n*b_n/theta)*chi)/(1-(M*l_n*b_n/theta))  # equation (34)

    function get_delta(d_n, a_n, l_n, b_n, L_n_1, O_n_1)
        Var = (L_n_1^2 - O_n_1^2)/L_n_1
        Cov = O_n_1*(L_n_1 - O_n_1)/L_n_1
        return d_n + a_n*l_n^2*su^2/theta^2 + a_n*(l_n^2*b_n^2/theta^2)*(M-1)*(Var + (M-2)*Cov)
    end


    #--- solve

    function solve_system(L_N)
        # variables for storing solution
        a = zeros(1,N+1)
        ps = zeros(1,N+1)
        m = zeros(1,N+1)
        d = zeros(1,N+1)
        l = zeros(1,N+1)
        g = zeros(1,N+1)
        b = zeros(1,N+1)
        d = zeros(1,N+1)
        eta = zeros(1,N+1)
        ph = zeros(1,N+1)
        LI = zeros(1,N+1)
        OI = zeros(1,N+1)
        SI = zeros(1,N+1)

        LI[N+1] = L_N
        OI[N+1] = L_N - chi
        SI[N+1] = (L_N + (M-1)*OI[N+1])*theta^2/M

        for n in [N+1:-1:2;]
            # compute lam_n from S_n & a_n
            eq_l(l) = eq_l_n(l, SI[n], ps[n], LI[n])
            l[n] = find_zero(eq_l, (0, SI[n]^(1/2)/su))

            b[n] = get_beta(l[n], SI[n])
            ph[n] = get_phi(OI[n], l[n], b[n], LI[n])
            g[n] = get_gamma(m[n], l[n], b[n])
            eta[n] = get_eta(ph[n])

            a[n-1] = get_alpha(eta[n], l[n], b[n], a[n], ph[n])
            ps[n-1] = get_psi(eta[n], l[n], b[n], a[n], ph[n], g[n], ps[n])
            m[n-1] = get_mu(l[n], g[n], b[n], m[n])

            LI[n-1] = get_L(LI[n], l[n], b[n])
            OI[n-1] = get_O(OI[n], l[n], b[n])
            SI[n-1] = (LI[n-1] + (M-1)*OI[n-1])*theta^2/M  # equation (8)
            d[n-1] = get_delta(d[n], a[n], l[n], b[n], LI[n-1], OI[n-1])
        end
        return LI[1] - L_0
    end

    #--- solve
    L_N_eqm[idx_chi] = find_zero(solve_system, ((M-1)/M*chi+1e-7, 10))


    #--- get solution
    a = zeros(1,N+1)
    ps = zeros(1,N+1)
    m = zeros(1,N+1)
    d = zeros(1,N+1)
    g = zeros(1,N+1)
    b = zeros(1,N+1)
    d = zeros(1,N+1)
    eta = zeros(1,N+1)
    ph = zeros(1,N+1)

    L[idx_chi,N+1] = L_N_eqm[idx_chi]
    O[idx_chi,N+1] = L_N_eqm[idx_chi] - chi
    S[idx_chi,N+1] = (L_N_eqm[idx_chi] + (M-1)*O[idx_chi,N+1])*theta^2/M

    for n in [N+1:-1:2;]
        # compute lam_n from S_n & a_n
        eq_l(l) = eq_l_n(l, S[idx_chi,n], ps[n], L[idx_chi,n])
        l[idx_chi,n] = find_zero(eq_l, (0, S[idx_chi,n]^(1/2)/su))

        b[n] = get_beta(l[idx_chi,n], S[idx_chi,n])
        ph[n] = get_phi(O[idx_chi,n], l[idx_chi,n], b[n], L[idx_chi,n])
        g[n] = get_gamma(m[n], l[idx_chi,n], b[n])
        eta[n] = get_eta(ph[n])

        a[n-1] = get_alpha(eta[n], l[idx_chi,n], b[n], a[n], ph[n])
        ps[n-1] = get_psi(eta[n], l[idx_chi,n], b[n], a[n], ph[n], g[n], ps[n])
        m[n-1] = get_mu(l[idx_chi,n], g[n], b[n], m[n])

        L[idx_chi,n-1] = get_L(L[idx_chi,n], l[idx_chi,n], b[n])
        O[idx_chi,n-1] = get_O(O[idx_chi,n], l[idx_chi,n], b[n])
        S[idx_chi,n-1] = (L[idx_chi,n-1] + (M-1)*O[idx_chi,n-1])*theta^2/M  # equation (8)
        d[n-1] = get_delta(d[n], a[n], l[idx_chi,n], b[n], L[idx_chi,n-1], O[idx_chi,n-1])
    end

end
#--- draw figures

# plot(l[2:N])

Corr_label = [latexstring("\\mathrm{Corr}=", trunc(O[idx_chi,1]/L[idx_chi,1],digits=4)) for idx_chi = 1:Nchi]

plot(S', label=Corr_label)
# title!("Information decay")
xaxis!(L"\mathrm{Trading \ time}")
yaxis!(L"\Sigma\mathrm{: \ Information \ decay}")

savefig("test.pdf")
# savefig("test.tex")
