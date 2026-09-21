# SELF_BUDGET: 300
# Cuatro láminas para el repo y el ensayo. Paleta validada con el validador del
# oficio (4 slots, light+dark: todos PASS; en claro el contraste de aqua/amarillo
# queda por debajo de 3:1, así que TODAS las series van con etiqueta directa
# además de leyenda — la regla de relieve).

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

D = '/root/agent-personality/projects/sonar_murcielago/'
res = json.load(open('/root/agent-personality/res_v232.json'))
ph = json.load(open('/root/agent-personality/res_v232_posthoc.json'))

SURF = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#4a4a48'
MUTED = '#8a8a86'
GRID = '#e4e4e0'
S1, S2, S3, S4 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'

plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF,
    'savefig.facecolor': SURF,
    'text.color': INK, 'axes.labelcolor': INK2, 'axes.edgecolor': GRID,
    'xtick.color': INK2, 'ytick.color': INK2,
    'font.size': 10.5, 'axes.titlesize': 12.5, 'axes.titleweight': 'bold',
    'axes.spines.top': False, 'axes.spines.right': False,
    'grid.color': GRID, 'grid.linewidth': 0.8,
    'lines.linewidth': 2.0, 'lines.solid_capstyle': 'round',
})


def tidy(ax, ytitle=None):
    ax.grid(True, axis='y', zorder=0)
    ax.set_axisbelow(True)
    if ytitle:
        ax.set_ylabel(ytitle)


# ── 1 · pérdida vs velocidad de acercamiento ───────────────────────────────────
sw = res['barrido_v']
v = np.array(sw['v'])
fig, ax = plt.subplots(figsize=(7.4, 4.6))
ax.axvspan(4, 7, color='#f1f1ee', zorder=0)
ax.text(5.5, 5.72, 'vuelo normal\nde un murciélago', ha='center', va='top',
        color=MUTED, fontsize=9, linespacing=1.3)
ax.axhline(1.0, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
ax.text(0.08, 1.12, '1 dB — presupuesto de desajuste del oficio',
        color=MUTED, fontsize=9)
ax.plot(v, sw['lfm_loss_db'], color=S2, label='LFM (barrido lineal)', zorder=3)
ax.plot(v, sw['hfm_loss_db'], color=S1, label='HFM (barrido hiperbólico)', zorder=3)
ax.annotate('LFM', (v[-1], sw['lfm_loss_db'][-1]), xytext=(-4, 6),
            textcoords='offset points', color=S2, fontweight='bold', ha='right')
ax.annotate('HFM', (v[-1], sw['hfm_loss_db'][-1]), xytext=(-4, -16),
            textcoords='offset points', color=S1, fontweight='bold', ha='right')
ax.annotate('aquí el máximo del filtro\nsalta de lóbulo', (9.75, 6.11),
            xytext=(8.15, 4.75), textcoords='data', color=MUTED, fontsize=8.5,
            ha='right', linespacing=1.3,
            arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.9,
                            shrinkA=2, shrinkB=3))
ax.set_xlim(0, 12); ax.set_ylim(0, 6.6)
ax.set_xlabel('velocidad de acercamiento (m/s)')
ax.set_title('Lo que cuesta no compensar el Doppler\n'
             'canto de 2 ms, 55→25 kHz (Eptesicus, fase de aproximación)')
tidy(ax, 'pérdida del filtro adaptado (dB)')
ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.015, 0.97))
fig.tight_layout(); fig.savefig(D + 'fig1_perdida_vs_v.png', dpi=170)
plt.close(fig)

# ── 2 · la inversión de ranking según el modelo de Doppler ─────────────────────
inv = ph['inversion']
vv = np.array(inv['v'])
fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4), sharey=True)
titulos = ['Doppler REAL\n(compresión del tiempo)',
           'Doppler APROXIMADO\n(desplazamiento en frecuencia)']
claves = [('wb_lfm', 'wb_hfm'), ('nb_lfm', 'nb_hfm')]
for ax, tit, (kl, kh) in zip(axes, titulos, claves):
    ax.plot(vv, inv[kl], color=S2, marker='o', ms=5, label='LFM')
    ax.plot(vv, inv[kh], color=S1, marker='o', ms=5, label='HFM')
    ax.annotate('LFM', (vv[-1], inv[kl][-1]), xytext=(-2, 8),
                textcoords='offset points', color=S2, fontweight='bold', ha='right')
    ax.annotate('HFM', (vv[-1], inv[kh][-1]), xytext=(-2, 8),
                textcoords='offset points', color=S1, fontweight='bold', ha='right')
    ax.set_title(tit, fontsize=11)
    ax.set_xlabel('velocidad de acercamiento (m/s)')
    ax.grid(True, axis='y'); ax.set_axisbelow(True)
    ax.set_xlim(0, 11)
axes[0].set_ylabel('pérdida del filtro adaptado (dB)')
axes[0].legend(frameon=False, loc='upper left')
fig.suptitle('La aproximación de banda estrecha elige el ganador OPUESTO '
             '(6 de 6 velocidades)', fontsize=12.5, fontweight='bold', y=1.0)
fig.tight_layout(); fig.savefig(D + 'fig2_inversion.png', dpi=170)
plt.close(fig)

# ── 3 · pérdida vs producto tiempo-ancho de banda ──────────────────────────────
tb = res['barrido_tb']
x = np.array(tb['TB'])
fig, ax = plt.subplots(figsize=(7.4, 4.4))
ax.plot(x, tb['lfm_loss_db'], color=S2, marker='o', ms=5, label='LFM')
ax.plot(x, tb['hfm_loss_db'], color=S1, marker='o', ms=5, label='HFM')
ax.annotate('LFM: crece sin techo', (x[-1], tb['lfm_loss_db'][-1]),
            xytext=(-6, -16), textcoords='offset points', color=S2,
            fontweight='bold', ha='right')
ax.annotate('HFM: 0,178 dB, plano', (x[-1], tb['hfm_loss_db'][-1]),
            xytext=(-6, 10), textcoords='offset points', color=S1,
            fontweight='bold', ha='right')
ax.axvline(60, color=MUTED, lw=1.0, ls=(0, (4, 3)))
ax.text(63, 12.2, 'canto de Eptesicus\n(TB = 60)', color=MUTED, fontsize=9,
        linespacing=1.3)
ax.set_xscale('log')
ax.minorticks_off()
ax.set_xticks([15, 30, 60, 120, 240, 600])
ax.set_xticklabels(['15', '30', '60', '120', '240', '600'])
ax.set_xlabel('producto tiempo × ancho de banda  (B = 30 kHz fijo, T variable)')
ax.set_title('La tolerancia del hiperbólico no depende de la duración;\n'
             'la del lineal se deshace con ella  (v = 5 m/s)')
tidy(ax, 'pérdida del filtro adaptado (dB)')
ax.legend(frameon=False, loc='upper left')
fig.tight_layout(); fig.savefig(D + 'fig3_perdida_vs_TB.png', dpi=170)
plt.close(fig)

# ── 4 · el sesgo de rango, contra la agudeza del propio murciélago ─────────────
bl = np.abs(np.array(sw['lfm_bias_us']))
bh = np.abs(np.array(sw['hfm_bias_us']))
fig, ax = plt.subplots(figsize=(7.4, 4.8))
ax.axhspan(60, 120, color='#efefec', zorder=0)
ax.axhspan(30, 40, color='#e6e6e2', zorder=0)
ax.text(0.35, 88, 'umbral de discriminación\nde rango (1–2 cm)', ha='left',
        va='center', color=MUTED, fontsize=9, linespacing=1.3)
ax.text(0.35, 34, 'idem, corregido por\nmovimiento de cabeza', ha='left',
        va='center', color=MUTED, fontsize=9, linespacing=1.3)
ax.axhline(1.0, color=MUTED, lw=1.2, ls=(0, (2, 2)))
ax.text(5.0, 1.25, 'agudeza de jitter: ~1 µs', color=MUTED, fontsize=9)
ax.plot(v, bh, color=S1, lw=3.4, label='HFM', zorder=3)
ax.plot(v, bl, color=S2, lw=1.8, ls=(0, (5, 2.5)), label='LFM', zorder=4)
ax.annotate('por encima de ~9,8 m/s el pico\ndel LFM ya no tiene posición:\nsalta de lóbulo',
            (10.6, 58), xytext=(9.4, 5.2), textcoords='data', color=MUTED,
            fontsize=8.5, ha='right', linespacing=1.35,
            arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.9,
                            shrinkA=2, shrinkB=3))
ax.set_yscale('log')
ax.minorticks_off()
ax.set_yticks([1, 10, 100])
ax.set_yticklabels(['1', '10', '100'])
ax.set_xlim(0, 12); ax.set_ylim(0.6, 400)
ax.set_xlabel('velocidad de acercamiento (m/s)')
ax.set_title('Lo que la tolerancia NO arregla: el sesgo de rango\n'
             'las dos familias pagan lo mismo (coinciden al 1,5 % a 5 m/s)')
tidy(ax, 'sesgo de retardo, |µs|   (1 cm ≈ 58 µs)')
ax.legend(frameon=False, loc='lower right')
fig.tight_layout(); fig.savefig(D + 'fig4_sesgo_rango.png', dpi=170)
plt.close(fig)

# ── 5 · tolerancia = ceguera: cuánta velocidad hace falta para perder 1 dB ─────
ceg = ph['ceguera_v_1db']
orden = ['CF-FM Rhinolophus', 'CF tono 2ms', 'LFM 2ms 55-25k', 'HFM 2ms 55-25k']
etiq = ['CF-FM de Rhinolophus\n(83 kHz, 51,5 ms)', 'tono puro CF\n(40 kHz, 2 ms)',
        'LFM\n(55→25 kHz, 2 ms)', 'HFM\n(55→25 kHz, 2 ms)']
val = [ceg[k] for k in orden]
col = [S4, S3, S2, S1]
fig, ax = plt.subplots(figsize=(7.6, 4.2))
y = np.arange(len(val))
# puntos, no barras: el eje es logarítmico y la LONGITUD de una barra en un eje
# log no es proporcional a nada. Un punto codifica posición, que sí vale.
for yi, vi, ci in zip(y, val, col):
    ax.plot([vi], [yi], 'o', ms=11, color=ci, zorder=3)
    ax.text(vi * 1.22, yi, ('%.2f m/s' % vi).replace('.', ','), va='center',
            color=INK, fontsize=10, fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels(etiq, fontsize=9.5)
ax.invert_yaxis()
ax.set_xscale('log'); ax.set_xlim(0.03, 60)
ax.minorticks_off()
ax.set_xticks([0.05, 0.1, 0.5, 1, 5, 10, 30])
ax.set_xticklabels(['0,05', '0,1', '0,5', '1', '5', '10', '30'])
ax.set_xlabel('cambio de velocidad que le cuesta 1 dB  (más a la derecha = más ciega)')
ax.set_title('Tolerancia es ceguera: lo que no te hace perder amplitud,\n'
             'tampoco te deja medirlo')
ax.grid(True, axis='x'); ax.set_axisbelow(True)
ax.spines['left'].set_visible(False)
ax.set_ylim(len(val) - 0.4, -0.6)
fig.tight_layout(); fig.savefig(D + 'fig5_ceguera.png', dpi=170)
plt.close(fig)

print("✓ 5 láminas escritas en", D)
