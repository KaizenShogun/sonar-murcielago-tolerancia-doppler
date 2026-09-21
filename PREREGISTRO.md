# PREREGISTRO V232 — Tolerancia Doppler de los cantos de murciélago (HFM vs LFM vs CF-FM)

Sellado ANTES de escribir la numérica. 2026-09-21.

## Pregunta

Un eco que vuelve de un blanco que se acerca a velocidad `v` no es el canto desplazado en
frecuencia: es el canto **comprimido en el tiempo** por η = (c+v)/(c−v). Con c = 343 m/s y v = 5 m/s
(vuelo normal de un murciélago), η = 1.0296 — un 3 %, enorme para un sonar.

La literatura (Altes & Titlebaum 1970, *Bat signals as optimally Doppler tolerant waveforms*) dice
que los barridos hiperbólicos (HFM) son invariantes al escalado: HFM(ηt) = HFM(t+Δ) salvo fase
constante. Lo que quiero MEDIR yo, a parámetros reales de murciélago, es otra cosa:

1. ¿Cuánto paga de verdad un LFM a parámetros de murciélago? El criterio clásico de tolerancia LFM
   es (2v/c)·T·B ≪ 1. Aquí vale **1.78**. No es ≪ 1. Predigo que el LFM se rompe.
2. Si el HFM es invariante, ¿DÓNDE se va el Doppler? Sospecho que no desaparece: se **traslada** de
   la amplitud al retardo (sesgo de rango). Tolerancia = relocalización, no cancelación.
3. ¿De qué está hecha la pérdida que SÍ tiene el HFM? Hipótesis: entera de borde de banda (el eco
   trae frecuencias ηf₁ que la réplica no contiene), no de fase.

## Instrumento

- Señales analíticas definidas por su fase; el eco se evalúa EXACTAMENTE en el argumento escalado
  `s(η(t−τ))` — sin remuestreo, sin interpolación de señal.
- fs = 2 MHz (rejilla 0.5 µs). Picos por interpolación parabólica sobre |envolvente| → resolución
  efectiva ≈ 0.05 µs, 400× por debajo de las diferencias en juego (~20 µs).
- Envolvente Tukey α=0.2 (principal) y rectangular (sensibilidad).
- Pérdida = 20·log10 de la correlación cruzada normalizada máxima (ambas señales a energía unidad).
  Suelo del instrumento: la auto-correlación a v=0 debe dar < 1e-6 dB.
- Modelo BANDA ANCHA (escalado) vs NARROWBAND (desplazamiento de frecuencia ν = (2v/c)·f_c,
  f_c = centro de banda) como comparación explícita.

## Señales (parámetros publicados)

- **LFM / HFM** tipo *Eptesicus fuscus* fase de aproximación: T = 2.0 ms, 55 → 25 kHz, B = 30 kHz,
  TB = 60. Mismo T, mismo f₁, mismo f₂ para las dos familias: sólo cambia la ley de fase.
- **CF-FM** tipo *Rhinolophus ferrumequinum*: CF 83 kHz durante 50 ms + cola FM 1.5 ms 83 → 68 kHz.
- **Control positivo de catástrofe**: tono puro CF de 2 ms a 40 kHz (debe perder mucho).

## Controles (#44: auditar el instrumento en los dos sentidos)

- NULO: v = 0 ⇒ pérdida 0 y sesgo 0 para las tres familias.
- POSITIVO analítico: ancho −3 dB del LFM ≈ 0.886/B; pendiente de la cresta LFM = T/B.
- POSITIVO de catástrofe: el tono CF de 2 ms pierde > 10 dB a v = 5 m/s (analítico: |sinc(νT)| con
  νT = 2.37).
- Mecanismos rivales declarados: (a) la envolvente podría dominar → mido rect y Tukey;
  (b) el signo del Doppler (acercarse/alejarse) → mido los dos; (c) la rejilla → repito el titular
  a fs = 1, 2 y 4 MHz y exijo que no mueva el resultado más de 0.05 dB.

## Predicciones selladas

**P1 [control+]** LFM, envolvente rect, v=0: ancho −3 dB en retardo dentro de ±15 % de
0.886/B = 29.53 µs → HIT si está en **[25.10, 33.96] µs**.

**P2 [control+]** Pendiente de cresta LFM en el modelo narrowband, ajuste por mínimos cuadrados
sobre ν ∈ [−2000, +2000] Hz en 21 puntos: dentro de ±5 % de T/B = 66.667 ns/Hz → HIT si está en
**[63.33, 70.00] ns/Hz** con R² > 0.99.

**P3 [ARRIESGADA — el criterio]** A v = 5 m/s, Tukey(0.2), modelo banda ancha, las TRES cláusulas:
pérdida LFM **> 1.0 dB**, pérdida HFM **< 0.5 dB**, y diferencia (LFM − HFM) **> 1.0 dB**.
Razón sellada: (2v/c)·T·B = 1.78, así que el LFM NO está en régimen tolerante.

**P4 [ARRIESGADA]** El modelo narrowband SUBESTIMA la pérdida del LFM a v = 5 m/s en más de
**2.0 dB** (pérdida_bandaancha − pérdida_narrowband > 2.0 dB).

**P5 [control de mi álgebra]** Sesgo de retardo a v = 5 m/s (posición del pico del filtro adaptado,
respecto del retardo verdadero): LFM en **[65, 95] µs**, HFM en **[38, 58] µs**, y razón
HFM/LFM en **[0.52, 0.70]**. (Derivado: Δ_HFM = (a/b)(1−η)/η con a/b = T·f₂/(f₁−f₂) = 1.667 ms
⇒ 47.9 µs; Δ_LFM ≈ (η−1)·f_c/k = 78.9 µs con f_c = 40 kHz, k = B/T.)

**P6 [ARRIESGADA — la descomposición]** Con réplica extendida al soporte escalado
(regenerada sobre [Δ, Δ+T/η] con la misma ley), la correlación normalizada del HFM sube a
**≥ 0.999** (pérdida ≤ 0.0044 dB) mientras la del LFM se queda **≤ 0.98** (pérdida ≥ 0.18 dB).
Es decir: toda la pérdida del HFM es borde de banda; la del LFM es fase y no se arregla extendiendo.

**P7 [ARRIESGADA]** CF-FM de *Rhinolophus* sin compensar, v = 5 m/s: pérdida total del canto
completo en **[12, 20] dB**, Y la pérdida de la cola FM sola **< 1.0 dB**.

## Lo que NO estoy midiendo (alcance)

Formas de onda idealizadas a parámetros publicados, no grabaciones. Sin ruido, sin absorción
atmosférica, sin directividad, sin armónicos, sin oído de murciélago (el murciélago no es un filtro
adaptado y nadie ha demostrado que lo sea). El veredicto es sobre **familias de forma de onda en el
régimen de parámetros del murciélago**, no sobre ningún animal concreto.
