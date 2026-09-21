# SELF_BUDGET: 900
# Tolerancia Doppler de cantos de murciélago: HFM vs LFM vs CF-FM.
# V232 · preregistro en PREREG_v232_sonar.md (apuestas #2172-#2178).
#
# El Doppler de un blanco que se acerca NO es un desplazamiento de frecuencia:
# es una COMPRESIÓN del tiempo por η = (c+v)/(c-v). Con c=343 y v=5 m/s, η=1.0296.
# Aquí se mide, no se aproxima: el eco se evalúa exactamente en el argumento escalado.
#
# 900 s de presupuesto: la parte caras es el CF-FM (51.5 ms a 2 MHz = 103k muestras)
# y los barridos (2 formas x 40 velocidades x 3 rejillas). Piloto medido: ~0.4 s por
# correlación de 4k muestras, ~0.2 s por una de 103k con FFT. Cota alta: ~200 s.

import json
import numpy as np
from scipy.signal import correlate

C = 343.0          # m/s, aire a ~25 C
OUT = {}


# ─────────────────────────── envolventes y leyes de fase ───────────────────────

def env_tukey(n, alpha=0.2):
    """Envolvente Tukey sobre n muestras. alpha=0 -> rectangular."""
    if alpha <= 0:
        return np.ones(n)
    w = np.ones(n)
    m = int(np.floor(alpha * (n - 1) / 2.0)) + 1
    if m < 2:
        return w
    t = np.arange(m)
    ramp = 0.5 * (1 - np.cos(np.pi * t / (m - 1)))
    w[:m] = ramp
    w[-m:] = ramp[::-1]
    return w


def env_cos(n, fs, rise_s=5e-4):
    """Rampa coseno de duración fija (para el CF-FM, donde Tukey se comería la cola)."""
    w = np.ones(n)
    m = min(int(rise_s * fs), n // 2)
    if m < 2:
        return w
    t = np.arange(m)
    ramp = 0.5 * (1 - np.cos(np.pi * t / (m - 1)))
    w[:m] = ramp
    w[-m:] = ramp[::-1]
    return w


class LFM:
    """f(t) = f1 + k t."""
    name = "LFM"

    def __init__(self, f1, f2, T):
        self.f1, self.f2, self.T = f1, f2, T
        self.k = (f2 - f1) / T

    def phase(self, t):
        return 2 * np.pi * (self.f1 * t + 0.5 * self.k * t * t)


class HFM:
    """f(t) = 1/(a + b t)  (periodo lineal en el tiempo)."""
    name = "HFM"

    def __init__(self, f1, f2, T):
        self.f1, self.f2, self.T = f1, f2, T
        self.a = 1.0 / f1
        self.b = (1.0 / f2 - 1.0 / f1) / T

    def phase(self, t):
        return (2 * np.pi / self.b) * np.log1p((self.b / self.a) * t)

    def delta_exact(self, eta):
        """El retardo al que el escalado por eta manda la ley: HFM(eta t)=HFM(t+D)."""
        return (self.a / self.b) * (1.0 - eta) / eta


class CFFM:
    """CF de duración Tcf a fcf, seguida de FM de Tfm que baja fcf -> f_end."""
    name = "CF-FM"

    def __init__(self, fcf, Tcf, f_end, Tfm):
        self.fcf, self.Tcf, self.f_end, self.Tfm = fcf, Tcf, f_end, Tfm
        self.T = Tcf + Tfm
        self.k = (f_end - fcf) / Tfm

    def phase(self, t):
        t = np.asarray(t, dtype=float)
        ph = np.where(t <= self.Tcf,
                      2 * np.pi * self.fcf * t,
                      2 * np.pi * (self.fcf * self.Tcf
                                   + self.fcf * (t - self.Tcf)
                                   + 0.5 * self.k * (t - self.Tcf) ** 2))
        return ph


class CF:
    """Tono puro (control positivo de catástrofe)."""
    name = "CF"

    def __init__(self, f0, T):
        self.f0, self.T = f0, T

    def phase(self, t):
        return 2 * np.pi * self.f0 * t


# ─────────────────────────── síntesis y medida ─────────────────────────────────

def emit(sig, fs, envfun):
    """El canto tal cual se emite, muestreado en [0, T)."""
    n = int(round(sig.T * fs))
    t = np.arange(n) / fs
    return envfun(n) * np.exp(1j * sig.phase(t))


def echo_wideband(sig, fs, eta, envfun):
    """Eco exacto: r(t) = s(eta t). Soporte [0, T/eta): se comprime de verdad."""
    n = int(round(sig.T / eta * fs))
    t = np.arange(n) / fs
    return envfun(n) * np.exp(1j * sig.phase(eta * t))


def echo_narrowband(sig, fs, nu, envfun):
    """Eco de la APROXIMACIÓN: mismo soporte, frecuencia desplazada nu."""
    n = int(round(sig.T * fs))
    t = np.arange(n) / fs
    return envfun(n) * np.exp(1j * (sig.phase(t) + 2 * np.pi * nu * t))


def matched(echo, replica, fs):
    """Correlación normalizada. Devuelve (pico en [0,1], retardo del pico en us)."""
    r = correlate(echo, replica, mode='full')
    denom = np.linalg.norm(echo) * np.linalg.norm(replica)
    mag = np.abs(r) / denom
    kmax = int(np.argmax(mag))
    # interpolación parabólica sobre el módulo (que ya es la envolvente)
    if 0 < kmax < len(mag) - 1:
        y0, y1, y2 = mag[kmax - 1], mag[kmax], mag[kmax + 1]
        den = (y0 - 2 * y1 + y2)
        dk = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        peak = y1 - 0.25 * (y0 - y2) * dk
    else:
        dk, peak = 0.0, mag[kmax]
    lag = (kmax + dk) - (len(replica) - 1)
    return float(min(peak, 1.0)), float(lag / fs * 1e6)


def loss_db(peak):
    return float(-20 * np.log10(max(peak, 1e-300)))


def width_3db(echo, replica, fs):
    """Ancho -3 dB (en amplitud: 1/sqrt(2)) del pico comprimido, en us."""
    r = np.abs(correlate(echo, replica, mode='full'))
    r = r / r.max()
    k = int(np.argmax(r))
    thr = 1 / np.sqrt(2)

    def cross(idx, step):
        i = idx
        while 0 < i < len(r) - 1 and r[i] > thr:
            i += step
        # interpolación lineal en el cruce
        a, b = r[i], r[i - step]
        if b == a:
            return float(i)
        return float(i - step + step * (b - thr) / (b - a))

    left, right = cross(k, -1), cross(k, +1)
    return float((right - left) / fs * 1e6)


def eta_of(v):
    """Factor de escala de ida y vuelta para velocidad de acercamiento v."""
    return (C + v) / (C - v)


# ─────────────────────────── parámetros de las señales ─────────────────────────

FS = 2_000_000.0
F1, F2, T_CALL = 55e3, 25e3, 2.0e-3          # Eptesicus fuscus, aproximación
B = F1 - F2                                   # 30 kHz
FC_ARITH = 0.5 * (F1 + F2)                    # 40 kHz
V_REF = 5.0                                   # m/s, vuelo normal
TUK = 0.2

lfm = LFM(F1, F2, T_CALL)
hfm = HFM(F1, F2, T_CALL)
cffm = CFFM(83e3, 50.0e-3, 68e3, 1.5e-3)      # Rhinolophus ferrumequinum
cf2ms = CF(40e3, 2.0e-3)

ev_t = lambda n: env_tukey(n, TUK)
ev_r = lambda n: np.ones(n)

print("η(5 m/s) = %.6f · (2v/c)·T·B = %.3f · TB = %.0f"
      % (eta_of(V_REF), (2 * V_REF / C) * T_CALL * B, T_CALL * B))
