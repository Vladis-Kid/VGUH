/* engine.js
 * Client-side port of hplc_sim/models.py + metrics.py, so the dashboard
 * renders instantly without requiring the Python backend to be running.
 * If a backend API URL is configured in the UI, app.js will prefer real
 * server-computed results and fall back to this engine otherwise.
 */
const Engine = (() => {

  // ---- erfc approximation (Abramowitz & Stegun 7.1.26, |error| < 1.5e-7) ----
  function erfc(x) {
    const z = Math.abs(x);
    const t = 1 / (1 + 0.5 * z);
    const tau = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
      t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
      t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? tau : 2 - tau;
  }

  function kIsocratic(k_w, S, phi) {
    const logK = Math.log10(Math.max(k_w, 1e-6)) - S * phi;
    return Math.pow(10, logK);
  }

  function vantHoff(kRef, Tref, T, dH) {
    const R = 8.314462618;
    const lnK = Math.log(Math.max(kRef, 1e-6)) - (dH / R) * (1 / T - 1 / Tref);
    return Math.exp(lnK);
  }

  // gradient_profile: [[t0,phi0],[t1,phi1],...] piecewise-linear
  function phiAt(profile, t) {
    if (t <= profile[0][0]) return profile[0][1];
    for (let i = 0; i < profile.length - 1; i++) {
      const [t1, p1] = profile[i], [t2, p2] = profile[i + 1];
      if (t >= t1 && t <= t2) {
        if (t2 === t1) return p2;
        const frac = (t - t1) / (t2 - t1);
        return p1 + frac * (p2 - p1);
      }
    }
    return profile[profile.length - 1][1];
  }

  function gradientRetentionTime(k_w, S, t0, tDwell, profile) {
    const dt = 0.002;
    let integral = 0, t = 0;
    const maxT = 200;
    while (t < maxT) {
      const tEff = Math.max(0, t - tDwell);
      const phi = phiAt(profile, tEff);
      const k = kIsocratic(k_w, S, phi);
      integral += dt / (t0 * Math.max(k, 1e-6));
      if (integral >= 1) return t + t0;
      t += dt;
    }
    return maxT;
  }

  // EMG peak evaluated on an array of t values (numerically stabilized as in models.py)
  function emgPeak(tArr, area, tR, sigma, tau) {
    if (Math.abs(tau) < 1e-6) tau = tau >= 0 ? 1e-6 : -1e-6;
    const sign = tau > 0 ? 1 : -1;
    const tauAbs = Math.abs(tau);
    const prefactor = area / (2 * tauAbs);
    const out = new Float64Array(tArr.length);
    for (let i = 0; i < tArr.length; i++) {
      const tt = sign > 0 ? tArr[i] : (2 * tR - tArr[i]);
      const z = (tR + sign * (sigma * sigma) / tauAbs - tt) / (sigma * Math.SQRT2);
      const exponent = (sigma * sigma) / (2 * tauAbs * tauAbs) + sign * (tR - tt) / tauAbs;
      const expClamped = Math.max(Math.min(exponent, 700), -700);
      let y = prefactor * Math.exp(expClamped) * erfc(z);
      if (!isFinite(y)) y = 0;
      out[i] = y;
    }
    return out;
  }

  function retentionFactor(tR, t0) { return (tR - t0) / Math.max(t0, 1e-9); }
  function selectivity(k1, k2) { const lo = Math.min(k1, k2), hi = Math.max(k1, k2); return hi / Math.max(lo, 1e-9); }
  function hetp(lengthMm, N) { return (lengthMm * 1000) / Math.max(N, 1e-9); }
  function resolution(tR1, tR2, wh1, wh2) { return 1.18 * Math.abs(tR2 - tR1) / Math.max(wh1 + wh2, 1e-9); }

  // Kozeny-Carman-based backpressure estimate (bar)
  function estimatePressureBar(lengthMm, idMm, particleUm, flowMlMin, porosity, viscosityCp) {
    const L = lengthMm / 1000;          // m
    const r = (idMm / 1000) / 2;        // m
    const dp = particleUm * 1e-6;       // m
    const eps = porosity;
    const eta = viscosityCp * 1e-3;     // Pa*s
    const Q = flowMlMin * 1e-6 / 60;    // m3/s
    const area = Math.PI * r * r;
    const u = Q / (area * eps);         // interstitial linear velocity m/s
    const dP = (180 * eta * L * u * Math.pow(1 - eps, 2)) / (dp * dp * Math.pow(eps, 3));
    return dP / 1e5; // Pa -> bar
  }

  function runSimulation(state) {
    const { column, mobilePhase, method, detector, compounds } = state;
    const T = mobilePhase.temperature_C + 273.15;
    const Tref = 298.15;

    const r_cm = (column.id_mm / 10) / 2;
    const l_cm = column.length_mm / 10;
    const volumeMl = Math.PI * r_cm * r_cm * l_cm;
    const porosity = 0.65;
    const deadVolumeMl = volumeMl * porosity;
    const t0 = deadVolumeMl / Math.max(mobilePhase.flow_ml_min, 1e-6);

    const nPoints = 4000;
    const tEnd = method.run_time_min;
    const t = new Float64Array(nPoints);
    for (let i = 0; i < nPoints; i++) t[i] = (i / (nPoints - 1)) * tEnd;
    const signal = new Float64Array(nPoints);

    let profile;
    if (method.mode === 'isocratic') {
      profile = [[0, method.isocratic_phi], [tEnd, method.isocratic_phi]];
    } else {
      profile = method.gradient_profile.map(p => [p.time_min, p.phi_B]);
    }

    const peaks = [];
    for (const c of compounds) {
      const kwT = vantHoff(c.k_w, Tref, T, c.dH_J_mol);
      let tR;
      if (method.mode === 'isocratic') {
        const k = kIsocratic(kwT, c.S, method.isocratic_phi);
        tR = t0 * (1 + k);
      } else {
        tR = gradientRetentionTime(kwT, c.S, t0, method.dwell_time_min, profile);
      }
      const k = retentionFactor(tR, t0);
      const N_nominal = column.N_per_m_nominal * (column.length_mm / 1000);
      const broadening = 1 / (1 + 0.02 * Math.max(k, 0));
      const N = Math.max(N_nominal * broadening, 100);
      const sigma = tR / Math.sqrt(Math.max(N, 1));
      const tau = c.tau_rel * sigma;
      const area = 1000 * (c.response_factor ?? 1);
      const y = emgPeak(t, area, tR, sigma, tau);
      for (let i = 0; i < nPoints; i++) signal[i] += y[i];

      const wHalf = 2.355 * sigma;
      const wBase = 4.0 * sigma;
      let height = 0;
      for (let i = 0; i < nPoints; i++) if (y[i] > height) height = y[i];

      peaks.push({
        name: c.name, tR, k, N, HETP_um: hetp(column.length_mm, N),
        area, height_mAU: height, w_half_min: wHalf, w_base_min: wBase,
        sigma, tau, start: tR - 3 * sigma, stop: tR + Math.max(4 * sigma, 6 * tau),
      });
    }

    peaks.sort((a, b) => a.tR - b.tR);
    for (let i = 0; i < peaks.length - 1; i++) {
      const p1 = peaks[i], p2 = peaks[i + 1];
      const k1 = p1.k > 0 ? p1.k : 1e-6, k2 = p2.k > 0 ? p2.k : 1e-6;
      p2.alpha_vs_prev = selectivity(k1, k2);
      p2.Rs_vs_prev = resolution(p1.tR, p2.tR, p1.w_half_min, p2.w_half_min);
      // avoid overlapping integration windows with neighbor
      const mid = (p1.tR + p2.tR) / 2;
      if (p1.stop > mid) p1.stop = mid;
      if (p2.start < mid) p2.start = mid;
    }

    // detector noise + baseline drift
    const noiseStd = detector.noise_std, drift = detector.baseline_drift_per_min;
    for (let i = 0; i < nPoints; i++) {
      const g1 = Math.random(), g2 = Math.random();
      const gauss = Math.sqrt(-2 * Math.log(g1 + 1e-12)) * Math.cos(2 * Math.PI * g2);
      signal[i] += gauss * noiseStd + drift * t[i];
    }

    const pressureBar = estimatePressureBar(column.length_mm, column.id_mm, column.particle_um,
      mobilePhase.flow_ml_min, porosity, 0.9);

    const Ns = peaks.map(p => p.N);
    const Rss = peaks.filter(p => p.Rs_vs_prev !== undefined).map(p => p.Rs_vs_prev);

    return {
      time_min: Array.from(t), signal_mAU: Array.from(signal), t0_min: t0,
      peaks, pressure_bar: pressureBar,
      N_mean: Ns.length ? Ns.reduce((a, b) => a + b, 0) / Ns.length : null,
      Rs_min: Rss.length ? Math.min(...Rss) : null,
    };
  }

  return { runSimulation, kIsocratic, vantHoff, gradientRetentionTime, emgPeak, estimatePressureBar };
})();
