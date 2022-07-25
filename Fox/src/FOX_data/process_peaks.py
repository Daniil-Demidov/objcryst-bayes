import numpy as np
from scipy.special import erf
from scipy import integrate
import math
import sys
# import blur_rho
# import converter

# changed so that F includes visible peaks only
# zero divisions not checked - be careful if see warning

# class Pattern:
#     def __init__(self, peaks, sigmas, tt1, tt2, blur_sigma, lam):
#         self.peaks = peaks
#         self.sigmas = sigmas
#         self.tt1 = tt1
#         self.tt2 = tt2
#         self.blur_sigma = blur_sigma
#         self.lam = lam

def tt2d(tt, lam=1.540562):
    return lam/2/math.sin(math.radians(tt/2))

def d2tt(d, lam=1.540562):
    return 2*math.degrees(math.asin(lam/2/d))

def x2tt(x, lam=1.540562):
    return d2tt(1/x, lam)

def tt2x(tt, lam=1.540562):
    return 1/tt2d(tt, lam)

def tt2q(tt, lam=1.540562):
    return tt2x(tt, lam)**2

def q2tt(q, lam=1.540562):
    return x2tt(math.sqrt(q), lam)

def mod_erf(x, s):
    return erf(x/math.sqrt(2)/s)

def jeffreys(x, mu, s1, s2):
    x -= mu
    if abs(x) < 0.05*s1:
        return 1/math.sqrt(2*math.pi)/math.log(s2/s1)*(1/s1 - 1/s2)
    return 1/2/math.log(s2/s1)/x*(mod_erf(x, s1) - mod_erf(x, s2))

def gaussian(x, mu, sig):
    return 1/(math.sqrt(2*math.pi)*sig)*math.exp(-((x - mu)/sig)**2/2)

def mirrored_gaussian(tt, mu, s, tt1, tt2):
    mu1 = 2*tt1 - mu
    mu2 = 2*tt2 - mu
    if tt1 <= tt <= tt2:
        ans = sum([gaussian(tt, mui, s) for mui in (mu1, mu, mu2)])
    else:
        ans = 0
    return ans

def calc_rho(tt):
    return sum([mirrored_gaussian(tt, peaks[i], blur_sigma, tt1, tt2)
                for i in range(peaks.size)])/(1 - ext_prob)

def calc_f(tt):
##    f = calc_rho(tt)*ext_prob
    f = 0
    f += sum([jeffreys(tt, peaks[i], sigmas[i], 10*sigmas[i])
                       for i in range(peaks.size)])
    return f

def make_FRho():
    N = 50000
    # N = 5000
    qmax = tt2q(tt2)
    qs = np.linspace(0, qmax, N)
    qs[0] = qs[1]*0.99
    tts = np.array([q2tt(q) for q in qs])
    tt_rhos = np.array([calc_rho(tt) for tt in tts])
    q_rhos = tt_rhos*lam**2/2/np.sin(np.radians(tts))*180/np.pi
    tt_fs = np.array([calc_f(tt) for tt in tts])
    q_fs = tt_fs*lam**2/2/np.sin(np.radians(tts))*180/np.pi
    min_index = np.searchsorted(tts, tt1)
    q_rhos[:min_index] = np.full(min_index, q_rhos[min_index])
    q_fs[:min_index] = np.full(min_index, 0)
    qs[0] = 0

    q_rho_exts = q_rhos.copy()
    q_rho_exts[min_index:] *= ext_prob

    Rhos = integrate.cumtrapz(q_rhos, qs, initial=0)
    Fs = integrate.cumtrapz(q_fs, qs, initial=0)
    Rho_exts = integrate.cumtrapz(q_rho_exts, qs, initial=0)
    np.savetxt("FOX_data/integrals.txt", np.vstack((qs, Rhos, Fs, Rho_exts)).T, fmt='%.10f', delimiter='\t', header="q\tRho\t\tF\tRho-ext")

def calc_fi(tt, peak, sigma):
    return jeffreys(tt, peak, sigma, 10*sigma)

def make_K_bounds():
    d_sigmas = []
    for peak, sigma in zip(peaks, sigmas):
        for sign in [-1, 1]:
            bound = peak + sign*30*sigma
            # assert calc_fi(bound, peak, sigma)/calc_rho(bound) < (1 - ext_prob)
        n = 100
        delta = 60*sigma/n
        bounds = [0, 0]
        for sign, index in zip([-1, 1], [0, 1]):
            x = peak
            while True:
                # false also when left-hand side is nan
                if calc_fi(x, peak, sigma)/calc_rho(x) < (1 - ext_prob):
                    break
                # shouldn't be here
                if abs(x - peak) > 30*sigma:
                    break
                x += sign*delta
            bounds[index] = x
    ##        print(peak, x)
    ##        print(calc_fi(x, peak, sigma)/calc_rho(x))
    ##        assert 0
        center, left, right = list(
            map(tt2d, [peak, bounds[0], bounds[1]]))
        sigma = (left - right) / 2
        d_sigmas.append([center, sigma])
    ##    exp_intervals.append([peak, bounds[0], bounds[1]])
    ##    exp_intervals.append(d_values)

    with open("FOX_data/exp_intervals.txt", 'w') as u:
        u.write("#\td\tsigma\n")
        for interval in d_sigmas:
            u.write("\t".join(map(str, interval)) + '\n')


file = sys.argv[1]
tt1, tt2, lam = list(map(float, sys.argv[2:]))

# python process_peaks.py input_peaks.txt 10 32 1.540562
# ./process_peaks input_peaks.txt 10 32 1.540562

# file = "input_peaks.txt"
# tt1 = 10
# tt2 = 32
# lam = 1.540562

vpeak_sigma = np.loadtxt(file)
vpeak_sigma = vpeak_sigma[vpeak_sigma[:, 0].argsort()]
peaks = vpeak_sigma[:, 0]
sigmas = vpeak_sigma[:, 1]

# peaks, sigmas = np.loadtxt(file, unpack=True)

lines = np.insert(peaks, 0, tt1)
blur_sigma = 1.5034880923743086*np.diff(lines).max()
# blur_sigma = 4

ext_prob = 0.5


make_FRho()
make_K_bounds()