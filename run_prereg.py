# SELF_BUDGET: 900
# Corre el preregistro V232 sobre las definiciones de chirp.py.
# (chirp.py se quedó solo con las leyes de fase y la medida, para que
#  posthoc.py pueda importarlas sin relanzar el experimento entero.)
import json
import numpy as np
import sys
sys.path.insert(0, '/root/agent-personality/projects/sonar_murcielago')
from chirp import *  # noqa: F403
from chirp import (LFM, HFM, CFFM, CF, emit, echo_wideband, echo_narrowband,
                   matched, loss_db, width_3db, eta_of, env_tukey, env_cos, C,
                   FS, F1, F2, T_CALL, B, FC_ARITH, V_REF, TUK,
                   lfm, hfm, cffm, cf2ms, ev_t, ev_r, OUT)

# ─────────────────────────── control NULO (v = 0) ──────────────────────────────

print("\n── CONTROL NULO (v=0): pérdida y sesgo deben ser 0 ──")
null = {}
for sig, ev, tag in ((lfm, ev_t, 'lfm'), (hfm, ev_t, 'hfm'),
                     (cffm, lambda n: env_cos(n, FS), 'cffm')):
    s = emit(sig, FS, ev)
    r0 = echo_wideband(sig, FS, 1.0, ev)
    pk, lag = matched(r0, s, FS)
    null[tag] = {'loss_db': loss_db(pk), 'bias_us': lag}
    print("  %-5s pérdida %.3e dB · sesgo %.4f us" % (tag, loss_db(pk), lag))
OUT['control_nulo'] = null


# ─────────────────────────── P1 · ancho -3 dB del LFM (rect) ───────────────────

s_rect = emit(lfm, FS, ev_r)
w = width_3db(s_rect, s_rect, FS)
OUT['p1_width_us'] = w
OUT['p1_teorico_us'] = 0.886 / B * 1e6
print("\n── P1 ── ancho -3 dB LFM rect = %.3f us (teórico 0.886/B = %.3f us)"
      % (w, 0.886 / B * 1e6))


# ─────────────────────────── P2 · pendiente de cresta (narrowband) ─────────────

nus = np.linspace(-2000, 2000, 21)
lags = []
for nu in nus:
    r = echo_narrowband(lfm, FS, nu, ev_t)
    _, lag = matched(r, emit(lfm, FS, ev_t), FS)
    lags.append(lag)
lags = np.array(lags)
A = np.vstack([nus, np.ones_like(nus)]).T
coef, res, *_ = np.linalg.lstsq(A, lags, rcond=None)
slope_us_per_hz = coef[0]
pred = A @ coef
r2 = 1 - np.sum((lags - pred) ** 2) / np.sum((lags - lags.mean()) ** 2)
OUT['p2_slope_ns_per_hz'] = float(abs(slope_us_per_hz) * 1000)
OUT['p2_slope_signo'] = float(np.sign(slope_us_per_hz))
OUT['p2_r2'] = float(r2)
OUT['p2_teorico_ns_per_hz'] = float(T_CALL / B * 1e9)
print("── P2 ── |pendiente| = %.3f ns/Hz (signo %+d) · T/B = %.3f ns/Hz · R2 = %.6f"
      % (abs(slope_us_per_hz) * 1000, np.sign(slope_us_per_hz), T_CALL / B * 1e9, r2))


# ─────────────────────────── P3/P4/P5 · a v = 5 m/s ────────────────────────────

eta = eta_of(V_REF)
nu_ref = (2 * V_REF / C) * FC_ARITH

print("\n── P3/P5 ── banda ancha, Tukey(%.1f), v = %.1f m/s (ν equiv = %.1f Hz)"
      % (TUK, V_REF, nu_ref))
wb = {}
for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
    s = emit(sig, FS, ev_t)
    r = echo_wideband(sig, FS, eta, ev_t)
    pk, lag = matched(r, s, FS)
    wb[tag] = {'peak': pk, 'loss_db': loss_db(pk), 'bias_us': lag,
               'bias_cm': lag * 1e-6 * C / 2 * 100}
    print("  %-4s pérdida %6.3f dB · sesgo %+8.2f us (%+.3f cm)"
          % (sig.name, loss_db(pk), lag, lag * 1e-6 * C / 2 * 100))

OUT['p3_loss_lfm_db'] = wb['lfm']['loss_db']
OUT['p3_loss_hfm_db'] = wb['hfm']['loss_db']
OUT['p5_bias_lfm_us'] = abs(wb['lfm']['bias_us'])
OUT['p5_bias_hfm_us'] = abs(wb['hfm']['bias_us'])
OUT['p5_ratio'] = abs(wb['hfm']['bias_us']) / abs(wb['lfm']['bias_us'])
OUT['p5_signo_lfm'] = float(np.sign(wb['lfm']['bias_us']))
OUT['p5_signo_hfm'] = float(np.sign(wb['hfm']['bias_us']))
OUT['p5_bias_cm_lfm'] = wb['lfm']['bias_cm']
OUT['p5_bias_cm_hfm'] = wb['hfm']['bias_cm']
OUT['p5_hfm_delta_analitico_us'] = hfm.delta_exact(eta) * 1e6
OUT['p5_lfm_delta_analitico_us'] = (eta - 1) * FC_ARITH / lfm.k * 1e6
print("  razón |sesgo HFM| / |sesgo LFM| = %.4f" % OUT['p5_ratio'])
print("  analítico: HFM %.2f us · LFM %.2f us"
      % (OUT['p5_hfm_delta_analitico_us'], OUT['p5_lfm_delta_analitico_us']))

# el mismo v con el modelo NARROWBAND
print("\n── P4 ── el mismo Doppler tratado como desplazamiento de frecuencia")
nb = {}
for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
    s = emit(sig, FS, ev_t)
    r = echo_narrowband(sig, FS, nu_ref, ev_t)
    pk, lag = matched(r, s, FS)
    nb[tag] = {'loss_db': loss_db(pk), 'bias_us': lag}
    print("  %-4s pérdida %6.3f dB · sesgo %+8.2f us" % (sig.name, loss_db(pk), lag))
OUT['p4_loss_lfm_nb_db'] = nb['lfm']['loss_db']
OUT['p4_loss_hfm_nb_db'] = nb['hfm']['loss_db']
OUT['p4_bias_lfm_nb_us'] = nb['lfm']['bias_us']
print("  brecha banda-ancha − narrowband (LFM) = %.3f dB"
      % (OUT['p3_loss_lfm_db'] - OUT['p4_loss_lfm_nb_db']))

# control de catástrofe: tono de 2 ms
s_cf = emit(cf2ms, FS, ev_t)
r_cf = echo_wideband(cf2ms, FS, eta, ev_t)
pk_cf, _ = matched(r_cf, s_cf, FS)
OUT['control_cf2ms_loss_db'] = loss_db(pk_cf)
print("\n── CONTROL+ catástrofe ── tono CF 2 ms a 40 kHz: pérdida %.2f dB (exigido >10)"
      % loss_db(pk_cf))

# control de signo: blanco que se ALEJA
print("── CONTROL signo ── blanco alejándose a 5 m/s:")
for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
    s = emit(sig, FS, ev_t)
    r = echo_wideband(sig, FS, eta_of(-V_REF), ev_t)
    pk, lag = matched(r, s, FS)
    OUT['signo_alejandose_' + tag + '_bias_us'] = lag
    OUT['signo_alejandose_' + tag + '_loss_db'] = loss_db(pk)
    print("  %-4s pérdida %6.3f dB · sesgo %+8.2f us" % (sig.name, loss_db(pk), lag))


# ─────────────────────────── P6 · descomposición (réplica continuada) ──────────
# Misma duración (T/eta), misma envolvente, y la LEY DE FASE DEL CANTO continuada
# a partir de un desfase D que se optimiza. Lo que quede es mismatch de LEY, no de
# soporte ni de envolvente.

def corr_ley_continuada(sig, fs, eta, envfun, span_us=400.0, step_us=0.5):
    n = int(round(sig.T / eta * fs))
    t = np.arange(n) / fs
    e = envfun(n)
    r = e * np.exp(1j * sig.phase(eta * t))
    nr = np.linalg.norm(r)
    best, bestD = -1.0, None
    Ds = np.arange(-span_us, span_us + step_us, step_us) * 1e-6
    for D in Ds:
        rep = e * np.exp(1j * sig.phase(t + D))
        c = abs(np.vdot(rep, r)) / (nr * np.linalg.norm(rep))
        if c > best:
            best, bestD = c, D
    # refino alrededor del mejor
    for D in np.arange(bestD - 0.5e-6, bestD + 0.5e-6, 0.01e-6):
        rep = e * np.exp(1j * sig.phase(t + D))
        c = abs(np.vdot(rep, r)) / (nr * np.linalg.norm(rep))
        if c > best:
            best, bestD = c, D
    return float(min(best, 1.0)), float(bestD * 1e6)


print("\n── P6 ── réplica con la ley continuada y duración igualada (T/η)")
# control de la construcción: sin escalar (eta=1) debe dar 1.000 en D=0
for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
    c1, d1 = corr_ley_continuada(sig, FS, 1.0, ev_t, span_us=20, step_us=0.5)
    OUT['p6_control_eta1_' + tag] = c1
    print("  control η=1 %-4s corr = %.9f en D = %+.2f us" % (sig.name, c1, d1))
for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
    c, d = corr_ley_continuada(sig, FS, eta, ev_t)
    OUT['p6_corr_' + tag + '_ext'] = c
    OUT['p6_delta_' + tag + '_us'] = d
    print("  %-4s corr = %.9f (pérdida %.4f dB) en D = %+.2f us"
          % (sig.name, c, loss_db(c), d))


# ─────────────────────────── P7 · CF-FM de Rhinolophus ─────────────────────────

print("\n── P7 ── CF-FM (83 kHz CF 50 ms + FM 1.5 ms 83→68 kHz), v = 5 m/s")
ev_c = lambda n: env_cos(n, FS)
s_c = emit(cffm, FS, ev_c)
r_c = echo_wideband(cffm, FS, eta, ev_c)
pk_c, lag_c = matched(r_c, s_c, FS)
OUT['p7_loss_cfm_db'] = loss_db(pk_c)
OUT['p7_bias_cfm_us'] = lag_c
print("  canto completo: pérdida %.3f dB · sesgo %+.2f us" % (loss_db(pk_c), lag_c))

fm_tail = LFM(83e3, 68e3, 1.5e-3)
s_f = emit(fm_tail, FS, ev_c)
r_f = echo_wideband(fm_tail, FS, eta, ev_c)
pk_f, lag_f = matched(r_f, s_f, FS)
OUT['p7_loss_fmtail_db'] = loss_db(pk_f)
OUT['p7_bias_fmtail_us'] = lag_f
print("  cola FM sola:   pérdida %.3f dB · sesgo %+.2f us" % (loss_db(pk_f), lag_f))
OUT['p7_energia_cola_frac'] = float(1.5 / 51.5)
OUT['p7_loss_esperada_por_energia_db'] = float(-20 * np.log10(np.sqrt(1.5 / 51.5)))
print("  esperado si la CF no aporta nada: %.2f dB"
      % OUT['p7_loss_esperada_por_energia_db'])


# ─────────────────────────── sensibilidades declaradas ─────────────────────────

print("\n── SENSIBILIDAD envolvente (rect vs Tukey 0.2), v = 5 ──")
sens_env = {}
for envname, ev in (('rect', ev_r), ('tukey0.2', ev_t)):
    row = {}
    for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
        s = emit(sig, FS, ev)
        r = echo_wideband(sig, FS, eta, ev)
        pk, lag = matched(r, s, FS)
        row[tag] = {'loss_db': loss_db(pk), 'bias_us': lag}
    sens_env[envname] = row
    print("  %-9s LFM %6.3f dB / %+7.2f us · HFM %6.3f dB / %+7.2f us"
          % (envname, row['lfm']['loss_db'], row['lfm']['bias_us'],
             row['hfm']['loss_db'], row['hfm']['bias_us']))
OUT['sens_envolvente'] = sens_env

print("── SENSIBILIDAD rejilla (fs 1/2/4 MHz), v = 5, Tukey ──")
sens_fs = {}
for fs in (1e6, 2e6, 4e6):
    row = {}
    for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
        s = emit(sig, fs, ev_t)
        r = echo_wideband(sig, fs, eta, ev_t)
        pk, lag = matched(r, s, fs)
        row[tag] = {'loss_db': loss_db(pk), 'bias_us': lag}
    sens_fs['%.0fMHz' % (fs / 1e6)] = row
    print("  %5.1f MHz  LFM %6.3f dB / %+7.2f us · HFM %6.3f dB / %+7.2f us"
          % (fs / 1e6, row['lfm']['loss_db'], row['lfm']['bias_us'],
             row['hfm']['loss_db'], row['hfm']['bias_us']))
OUT['sens_rejilla'] = sens_fs


# ─────────────────────────── barrido en velocidad ──────────────────────────────

print("\n── BARRIDO v = 0..12 m/s ──")
vs = np.round(np.arange(0, 12.01, 0.25), 3)
sweep = {'v': vs.tolist()}
for sig, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
    s = emit(sig, FS, ev_t)
    L, Bi = [], []
    for v in vs:
        r = echo_wideband(sig, FS, eta_of(v), ev_t)
        pk, lag = matched(r, s, FS)
        L.append(loss_db(pk)); Bi.append(lag)
    sweep[tag + '_loss_db'] = L
    sweep[tag + '_bias_us'] = Bi
OUT['barrido_v'] = sweep

# velocidad a la que cada forma cruza 1 dB de pérdida
def cruce(vs, L, thr=1.0):
    L = np.array(L)
    idx = np.where(L >= thr)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(vs[0])
    v0, v1, l0, l1 = vs[i - 1], vs[i], L[i - 1], L[i]
    return float(v0 + (thr - l0) * (v1 - v0) / (l1 - l0))

OUT['v_cruce_1db_lfm'] = cruce(vs, sweep['lfm_loss_db'])
OUT['v_cruce_1db_hfm'] = cruce(vs, sweep['hfm_loss_db'])
OUT['v_criterio_clasico'] = float(C / (2 * T_CALL * B))
print("  LFM cruza 1 dB a v = %s m/s · HFM a %s m/s · criterio c/(2TB) = %.3f m/s"
      % (OUT['v_cruce_1db_lfm'], OUT['v_cruce_1db_hfm'], OUT['v_criterio_clasico']))


# ─────────────────────────── barrido en producto TB ────────────────────────────
# ¿Dónde vive el murciélago en el plano (TB, v)? Mantengo B=30 kHz y muevo T.

print("\n── BARRIDO TB (B fijo 30 kHz, T variable), v = 5 m/s ──")
Ts = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 14.0, 20.0]) * 1e-3
tb = {'T_ms': (Ts * 1e3).tolist(), 'TB': (Ts * B).tolist()}
for cls, tag in ((LFM, 'lfm'), (HFM, 'hfm')):
    L = []
    for T in Ts:
        sg = cls(F1, F2, T)
        s = emit(sg, FS, ev_t)
        r = echo_wideband(sg, FS, eta, ev_t)
        pk, _ = matched(r, s, FS)
        L.append(loss_db(pk))
    tb[tag + '_loss_db'] = L
for i, T in enumerate(Ts):
    print("  T=%5.1f ms (TB=%5.0f, (2v/c)TB=%6.2f)  LFM %7.3f dB · HFM %6.3f dB"
          % (T * 1e3, T * B, (2 * V_REF / C) * T * B,
             tb['lfm_loss_db'][i], tb['hfm_loss_db'][i]))
OUT['barrido_tb'] = tb


# ─────────────────────────── contra la agudeza del murciélago ──────────────────

OUT['agudeza'] = {
    'jitter_us': 1.0,               # Simmons: cambios de retardo de ~1 us o menos
    'jnd_rango_us': [60.0, 120.0],  # 1-2 cm de umbral de discriminación de rango
    'jnd_con_cabeza_us': [30.0, 40.0],
}
OUT['meta'] = {
    'c_m_s': C, 'fs_Hz': FS, 'f1_Hz': F1, 'f2_Hz': F2, 'T_s': T_CALL,
    'TB': T_CALL * B, 'v_ref_m_s': V_REF, 'eta_ref': eta,
    'criterio_2vc_TB': (2 * V_REF / C) * T_CALL * B,
    'envolvente': 'tukey0.2 (rect como sensibilidad)',
    'cffm': {'fcf': 83e3, 'Tcf': 50e-3, 'f_end': 68e3, 'Tfm': 1.5e-3},
}

with open('/root/agent-personality/res_v232.json', 'w') as f:
    json.dump(OUT, f, indent=1)
print("\n✓ res_v232.json escrito")
