
# Compute the solution for the difference equations of the market equilibrium in Holden & Subrahmanyam (1994, EL)
# Check the effect of liquidity order variance (su2) on price discovery (S)
#
# by Wataru Nozawa
# June 18, 2019

using Roots

using Plots
gr()

# Parameters

N = 20  # number of auctions
dt = 1/N  # time interval
A = 0  # risk aversion parameter
su2s = [0.25 0.5 1 2 4]  # variance of liquidity trader order size
S_0 = 1  # asset value variance
M = 1  # number of informed traders

# boundary condition for a_N

a_N = 0


# Intial guess for S_N

S_N = 0.02

# variables for storing solution
N_su2 = length(su2s)
a = zeros(N_su2,N)
S = zeros(N_su2,N)
l = zeros(N_su2,N)


# solve

for idx_su2 in 1:N_su2
    # functions for solution
    su2 = su2s[idx_su2]

    for_a_n_1(a_n, l_n) = (1-a_n*l_n + A*l_n*su2*dt/2)/(l_n*(M*(1-2*a_n*l_n) + 1 + A*l_n*su2*dt)^2)  # equation (7)

    for_l_n(l_n, a_n, S_n) = l_n^3*su2*(2*M*a_n-A*su2*dt)*dt - l_n^2*su2*(M+1)*dt - 2*M*a_n*S_n*l_n + M*S_n  # equation (23)

    for_S_n_1(l_n, S_n) = S_n/(1-l_n^2*su2*dt/S_n)  # equation (22 & 21)



    a[idx_su2,N] = a_N

    function diff_system(S_N0)
        S[idx_su2,N] = S_N0

        # repeat the following for n=N, N-1, ..., 2
        for n in [N:-1:2;]
            # compute lam_n from S_n & a_n
            for_l(l) = for_l_n(l, a[idx_su2,n],S[idx_su2,n])
            l[idx_su2,n] = find_zero(for_l, 0.5)

            # compute S_n-1 from lam_n & S_n
            S[idx_su2,n-1] = for_S_n_1(l[idx_su2,n], S[idx_su2,n])

            # compute a_n-1 from lam_n & a_n
            a[idx_su2,n-1] = for_a_n_1(a[idx_su2,n], l[idx_su2,n])
        end
        return S[idx_su2,1] - S_0
    end
    # solve

    S_N_eqm = find_zero(diff_system, 0.02)
    S[N] = S_N_eqm

    # repeat the following for n=N, N-1, ..., 2
    for n in [N:-1:2;]
        # compute lam_n from S_n & a_n
        for_l(l) = for_l_n(l, a[idx_su2,n],S[idx_su2,n])
        l[idx_su2,n] = find_zero(for_l, 0.5)

        # compute S_n-1 from lam_n & S_n
        S[idx_su2,n-1] = for_S_n_1(l[idx_su2,n], S[idx_su2,n])

        # compute a_n-1 from lam_n & a_n
        a[idx_su2,n-1] = for_a_n_1(a[idx_su2,n], l[idx_su2,n])
    end
end


# draw figures

plot(l[:,2:N]')

plot(S')
