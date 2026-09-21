# SELF_BUDGET: 900
# POST-HOC (declarado: NO estaba preregistrado; nace de que P5 falló y la razón
# medida fue 0.985 en vez de 0.61). Tres cosas:
#  A) ¿es LEY que el sesgo de rango sea el MISMO para LFM y HFM? Δ ≈ ε·T·f_fin/B.
#     Se prueba fuera del punto donde nació: subidas y bajadas, 12 combinaciones.
#  B) la inversión de ranking narrowband vs banda ancha, a lo largo de v.
#  C) tolerancia = ceguera: cuánta velocidad hay que cambiar para que el filtro
#     lo note (1 dB). Lo que no pierdes amplitud, no lo puedes medir.

import json
import numpy as np
import sys
sys.path.insert(0, '/root/agent-personality/projects/sonar_murcielago')
from chirp import (LFM, HFM, CFFM, CF, emit, echo_wideband, echo_narrowband,
                   matched, loss_db, eta_of, env_tukey, env_cos, C)

FS = 2_000_000.0
ev = lambda n: env_tukey(n, 0.2)
OUT = {}

# ─── A) la ley del sesgo ───────────────────────────────────────────────────────
print("── A) Δ ≈ ε·T·f_fin/B fuera del punto donde nació ──")
print("  (ε = η-1; f_fin = frecuencia al final del barrido; B = |f1-f2|)")
casos = [
    (55e3, 25e3, 2.0e-3, 'Eptesicus aprox (bajada)'),
    (100e3, 40e3, 1.0e-3, 'bajada ancha corta'),
    (40e3, 20e3, 5.0e-3, 'bajada larga'),
    (25e3, 55e3, 2.0e-3, 'SUBIDA (espejo del 1)'),
    (30e3, 90e3, 3.0e-3, 'subida ancha'),
    (60e3, 50e3, 2.0e-3, 'bajada estrecha (B=10k)'),
]
tabla = []
for v in (2.0, 5.0, 9.0):
    eta = eta_of(v)
    eps = eta - 1.0
    for f1, f2, T, etiq in casos:
        B = abs(f1 - f2)
        pred = eps * T * f2 / B * 1e6          # us, primer orden
        row = {'v': v, 'f1': f1, 'f2': f2, 'T': T, 'etiq': etiq,
               'pred_us': pred}
        for cls, tag in ((LFM, 'lfm'), (HFM, 'hfm')):
            sg = cls(f1, f2, T)
            s = emit(sg, FS, ev)
            r = echo_wideband(sg, FS, eta, ev)
            pk, lag = matched(r, s, FS)
            row[tag + '_bias_us'] = lag
            row[tag + '_loss_db'] = loss_db(pk)
        row['err_lfm'] = abs(abs(row['lfm_bias_us']) - abs(pred)) / abs(pred)
        row['err_hfm'] = abs(abs(row['hfm_bias_us']) - abs(pred)) / abs(pred)
        tabla.append(row)
        print("  v=%4.1f  %-26s pred %8.2f us · LFM %8.2f (%5.1f%%) · HFM %8.2f (%5.1f%%)"
              % (v, etiq, pred, row['lfm_bias_us'], 100 * row['err_lfm'],
                 row['hfm_bias_us'], 100 * row['err_hfm']))
OUT['ley_sesgo'] = tabla
el = np.array([r['err_lfm'] for r in tabla])
eh = np.array([r['err_hfm'] for r in tabla])
raz = np.array([abs(r['hfm_bias_us']) / abs(r['lfm_bias_us']) for r in tabla])
OUT['ley_sesgo_resumen'] = {
    'n': len(tabla),
    'err_rel_lfm_mediana': float(np.median(el)), 'err_rel_lfm_max': float(el.max()),
    'err_rel_hfm_mediana': float(np.median(eh)), 'err_rel_hfm_max': float(eh.max()),
    'razon_hfm_lfm_min': float(raz.min()), 'razon_hfm_lfm_max': float(raz.max()),
    'razon_hfm_lfm_mediana': float(np.median(raz)),
}
print("  resumen: error rel. mediano LFM %.2f%% (máx %.2f%%) · HFM %.2f%% (máx %.2f%%)"
      % (100 * np.median(el), 100 * el.max(), 100 * np.median(eh), 100 * eh.max()))
print("  razón |sesgo HFM|/|sesgo LFM| en los %d casos: %.3f – %.3f (mediana %.3f)"
      % (len(tabla), raz.min(), raz.max(), np.median(raz)))

# ─── B) la inversión de ranking ────────────────────────────────────────────────
print("\n── B) narrowband vs banda ancha: quién gana, según el modelo ──")
lfm = LFM(55e3, 25e3, 2.0e-3)
hfm = HFM(55e3, 25e3, 2.0e-3)
FCA = 40e3
inv = {'v': [], 'wb_lfm': [], 'wb_hfm': [], 'nb_lfm': [], 'nb_hfm': []}
for v in (1.0, 2.0, 3.0, 5.0, 7.0, 10.0):
    eta = eta_of(v); nu = (2 * v / C) * FCA
    row = []
    for sg, tag in ((lfm, 'lfm'), (hfm, 'hfm')):
        s = emit(sg, FS, ev)
        pw, _ = matched(echo_wideband(sg, FS, eta, ev), s, FS)
        pn, _ = matched(echo_narrowband(sg, FS, nu, ev), s, FS)
        inv['wb_' + tag].append(loss_db(pw))
        inv['nb_' + tag].append(loss_db(pn))
        row += [loss_db(pw), loss_db(pn)]
    inv['v'].append(v)
    print("  v=%4.1f  BANDA ANCHA: LFM %6.3f  HFM %6.3f  → gana %s"
          "   |  NARROWBAND: LFM %6.3f  HFM %6.3f  → gana %s"
          % (v, row[0], row[2], 'HFM' if row[2] < row[0] else 'LFM',
             row[1], row[3], 'HFM' if row[3] < row[1] else 'LFM'))
OUT['inversion'] = inv
n_rev = sum(1 for i in range(len(inv['v']))
            if (inv['wb_hfm'][i] < inv['wb_lfm'][i]) != (inv['nb_hfm'][i] < inv['nb_lfm'][i]))
OUT['inversion_n_reversiones'] = n_rev
print("  los dos modelos eligen ganador OPUESTO en %d de %d velocidades"
      % (n_rev, len(inv['v'])))

# ─── C) tolerancia = ceguera ───────────────────────────────────────────────────
print("\n── C) ¿cuánta velocidad hay que cambiar para que el filtro pierda 1 dB? ──")
cffm = CFFM(83e3, 50e-3, 68e3, 1.5e-3)
ev_c = lambda n: env_cos(n, FS)
cf2 = CF(40e3, 2.0e-3)


def v_para_perder(sig, envf, thr_db=1.0, vmax=40.0, paso=0.05):
    s = emit(sig, FS, envf)
    v = 0.0
    while v <= vmax:
        pk, _ = matched(echo_wideband(sig, FS, eta_of(v), envf), s, FS)
        if loss_db(pk) >= thr_db:
            return float(v)
        v += paso
    return None


ceg = {}
for sig, envf, tag in ((lfm, ev, 'LFM 2ms 55-25k'), (hfm, ev, 'HFM 2ms 55-25k'),
                       (cf2, ev, 'CF tono 2ms'), (cffm, ev_c, 'CF-FM Rhinolophus')):
    vv = v_para_perder(sig, envf, 1.0)
    ceg[tag] = vv
    print("  %-20s → %s m/s" % (tag, ('%.2f' % vv) if vv is not None else '>40 (ciega)'))
OUT['ceguera_v_1db'] = ceg

with open('/root/agent-personality/res_v232_posthoc.json', 'w') as f:
    json.dump(OUT, f, indent=1)
print("\n✓ res_v232_posthoc.json escrito")
