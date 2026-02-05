/* =========================================================
   app.js (COMPLETO)
   - Reactivos/Productos ILIMITADOS (añadir/quitar filas)
   - Ajuste de reacción (coeficientes mínimos)
   - Masas molares
   - Cálculo: conocido (masa/moles/volumen disolución/gas) -> objetivo (masa/moles/volumen disolución/gas)
   - Conversión “modo cuaderno” con fracciones + tachado + colores
   ========================================================= */

document.addEventListener("DOMContentLoaded", () => {

  // ===============================
  //  Formato español (coma) + parser que acepta coma en inputs
  // ===============================
  function fmt(num, decimals) {
    return Number(num).toLocaleString("es-ES", {
      minimumFractionDigits: decimals,
      maximumFractionDigits: decimals
    });
  }
  function toNumber(val) {
    // Acepta "2,5" o "2.5"
    if (val === null || val === undefined) return NaN;
    if (typeof val === "number") return val;
    const s = String(val).trim().replace(/\s+/g, "").replace(",", ".");
    const n = Number(s);
    return n;
  }

  // ===============================
  //  Normalizar subíndices Unicode
  // ===============================
  function normalizeFormula(str) {
    return (str || "")
      .replace(/₀/g, "0").replace(/₁/g, "1").replace(/₂/g, "2")
      .replace(/₃/g, "3").replace(/₄/g, "4").replace(/₅/g, "5")
      .replace(/₆/g, "6").replace(/₇/g, "7").replace(/₈/g, "8")
      .replace(/₉/g, "9");
  }

  // ===============================
  //  Tabla periódica (masas atómicas)
  // ===============================
  const ATOMIC_MASSES = {
    H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81, C: 12.011, N: 14.007, O: 15.999,
    F: 18.998, Ne: 20.18, Na: 22.99, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974, S: 32.06,
    Cl: 35.45, Ar: 39.948, K: 39.098, Ca: 40.078, Sc: 44.956, Ti: 47.867, V: 50.942, Cr: 51.996,
    Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38, Ga: 69.723, Ge: 72.63,
    As: 74.922, Se: 78.971, Br: 79.904, Kr: 83.798, Rb: 85.468, Sr: 87.62, Ag: 107.8682,
    Cd: 112.414, Sn: 118.71, Sb: 121.76, I: 126.904, Ba: 137.327, Pt: 195.084, Au: 196.967,
    Hg: 200.592, Pb: 207.2
  };

  // ===============================
  //  Enteros: gcd/lcm + fracciones
  // ===============================
  function gcdInt(a, b) {
    a = Math.abs(Math.trunc(a));
    b = Math.abs(Math.trunc(b));
    while (b !== 0) {
      const t = b;
      b = a % b;
      a = t;
    }
    return a || 1;
  }

  function lcmInt(a, b) {
    a = Math.abs(Math.trunc(a));
    b = Math.abs(Math.trunc(b));
    if (a === 0 || b === 0) return 0;
    return (a / gcdInt(a, b)) * b;
  }

  // Aproxima x como fracción p/q (continued fractions)
  function toFraction(x, maxDen = 4000, eps = 1e-12) {
    if (!isFinite(x)) throw new Error("Número no finito en el ajuste.");
    if (Math.abs(x) < eps) return [0, 1];

    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);

    let a = Math.floor(x);
    let h1 = 1, k1 = 0;
    let h = a, k = 1;

    let frac = x - a;
    while (k <= maxDen && Math.abs(h / k - x) > eps && frac > eps) {
      frac = 1 / frac;
      a = Math.floor(frac);

      const h2 = h1, k2 = k1;
      h1 = h; k1 = k;
      h = a * h1 + h2;
      k = a * k1 + k2;

      frac = frac - a;
    }
    return [sign * h, k];
  }

  // ===============================
  //  Parser de fórmulas (paréntesis)
  // ===============================
  function parseFormula(formula) {
    formula = normalizeFormula(formula);
    const cleaned = formula.replace(/\s+/g, "");
    const stack = [new Map()];
    let i = 0;

    while (i < cleaned.length) {
      const ch = cleaned[i];

      if (ch === "(") {
        stack.push(new Map());
        i++;
        continue;
      }

      if (ch === ")") {
        i++;
        let num = "";
        while (i < cleaned.length && /[0-9]/.test(cleaned[i])) num += cleaned[i++];
        const mult = num ? parseInt(num, 10) : 1;

        const group = stack.pop();
        const top = stack[stack.length - 1];
        for (const [el, c] of group.entries()) {
          top.set(el, (top.get(el) || 0) + c * mult);
        }
        continue;
      }

      if (/[A-Z]/.test(ch)) {
        let sym = ch;
        i++;
        if (i < cleaned.length && /[a-z]/.test(cleaned[i])) sym += cleaned[i++];

        let num = "";
        while (i < cleaned.length && /[0-9]/.test(cleaned[i])) num += cleaned[i++];
        const count = num ? parseInt(num, 10) : 1;

        const top = stack[stack.length - 1];
        top.set(sym, (top.get(sym) || 0) + count);
        continue;
      }

      break; // no soportamos cargas, hidratos, etc.
    }

    return stack[0];
  }

  function molarMassFromFormula(f) {
    const comp = parseFormula(f);
    let mass = 0;
    for (const [el, c] of comp.entries()) {
      if (!ATOMIC_MASSES[el]) throw new Error("Elemento no soportado: " + el);
      mass += ATOMIC_MASSES[el] * c;
    }
    return mass;
  }

  // ===============================
  //  Matriz elementos × especies
  // ===============================
  function buildMatrix(species) {
    const elements = new Set();
    const parsed = species.map(s => ({ ...s, comp: parseFormula(s.formula) }));
    parsed.forEach(s => s.comp.forEach((_, el) => elements.add(el)));

    const elems = [...elements];
    return elems.map(el =>
      parsed.map(s => (s.comp.get(el) || 0) * (s.side === "reactivo" ? -1 : 1))
    );
  }

  // ===============================
  //  Espacio nulo -> coeficientes mínimos
  // ===============================
  function nullspaceVector(A) {
    const m = A.length;
    const n = A[0].length;

    const M = A.map(row => row.slice());
    const pivotCols = [];
    let row = 0;

    // RREF Gauss-Jordan
    for (let col = 0; col < n && row < m; col++) {
      let pivot = -1;
      for (let r = row; r < m; r++) {
        if (Math.abs(M[r][col]) > 1e-12) { pivot = r; break; }
      }
      if (pivot === -1) continue;

      [M[row], M[pivot]] = [M[pivot], M[row]];

      const pv = M[row][col];
      for (let c = col; c < n; c++) M[row][c] /= pv;

      for (let r = 0; r < m; r++) {
        if (r !== row && Math.abs(M[r][col]) > 1e-12) {
          const f = M[r][col];
          for (let c = col; c < n; c++) M[r][c] -= f * M[row][c];
        }
      }

      pivotCols.push(col);
      row++;
    }

    const pivotSet = new Set(pivotCols);
    const freeCols = [];
    for (let c = 0; c < n; c++) if (!pivotSet.has(c)) freeCols.push(c);

    if (freeCols.length === 0) throw new Error("No se puede ajustar la reacción con las especies dadas.");

    const f = freeCols[0];
    const v = Array(n).fill(0);
    v[f] = 1;

    for (let r = 0; r < pivotCols.length; r++) {
      const pc = pivotCols[r];
      let sum = 0;
      for (let c = pc + 1; c < n; c++) sum += M[r][c] * v[c];
      v[pc] = -sum;
    }

    const fracs = v.map(x => toFraction(x, 4000));
    let L = 1;
    for (const [, q] of fracs) L = lcmInt(L, q);

    let ints = fracs.map(([p, q]) => p * (L / q));

    let g = 0;
    for (const x of ints) g = gcdInt(g, Math.abs(x));
    ints = ints.map(x => x / (g || 1));

    if (!ints.some(x => x > 0)) ints = ints.map(x => -x);
    if (ints.some(x => x === 0)) throw new Error("Ajuste degenerado: revisa las especies.");

    return ints;
  }

  // ===============================
  //  Gas ideal
  // ===============================
  const R = 0.082057; // L·atm/(mol·K)
  function celsiusToK(tc) { return tc + 273.15; }
  function pressureToAtm(P, unit) {
    if (unit === "atm") return P;
    if (unit === "kPa") return P / 101.325;
    if (unit === "mmHg") return P / 760.0;
    return P;
  }

  // ===============================
  //  Conversión “modo cuaderno” (tachado + colores)
  // ===============================
  function unitSpan(unitText) {
    const u = (unitText || "").trim().replace(/\s+/g, " ");
    return `<span class="u" data-u="${u}">${u}</span>`;
  }
  function termHTML(valueText, unitText) {
    return `<span class="term"><span class="v">${valueText}</span> ${unitSpan(unitText)}</span>`;
  }
  function fractionHTML(topValue, topUnit, bottomValue, bottomUnit) {
    return `
      <span class="fraction">
        <span class="top">${termHTML(topValue, topUnit)}</span>
        <span class="bar"></span>
        <span class="bottom">${termHTML(bottomValue, bottomUnit)}</span>
      </span>
    `;
  }

  function cancelUnitsIn(containerEl) {
    const left = containerEl.querySelector(".conv-left");
    if (!left) return;

    const tops = Array.from(left.querySelectorAll(".top .u, .factor .u"));
    const bottoms = Array.from(left.querySelectorAll(".bottom .u"));

    const mapTop = new Map();
    const mapBottom = new Map();

    function push(map, el) {
      const key = (el.dataset.u || "").trim();
      if (!key) return;
      if (!map.has(key)) map.set(key, []);
      map.get(key).push(el);
    }

    tops.forEach(el => push(mapTop, el));
    bottoms.forEach(el => push(mapBottom, el));

    let colorIndex = 0;
    for (const [unit, arrTop] of mapTop.entries()) {
      const arrBot = mapBottom.get(unit);
      if (!arrBot || !arrBot.length) continue;

      const k = Math.min(arrTop.length, arrBot.length);
      const cls = `c${colorIndex % 8}`;
      colorIndex++;

      for (let i = 0; i < k; i++) {
        arrTop[i].classList.add("cancel", cls);
        arrBot[i].classList.add("cancel", cls);
      }
    }
  }

  // ===============================
  //  DOM: listas dinámicas de especies
  // ===============================
  const reactantsList = document.getElementById("reactants-list");
  const productsList = document.getElementById("products-list");
  const addReactantBtn = document.getElementById("add-reactant");
  const addProductBtn = document.getElementById("add-product");

  // Tarjeta 1
  const btnBalance = document.getElementById("btn-balance");
  const balanceError = document.getElementById("balance-error");
  const balancedInfo = document.getElementById("balanced-info");
  const balancedEquationText = document.getElementById("balanced-equation-text");
  const molarMassList = document.getElementById("molar-mass-list");
  const resetReaction = document.getElementById("reset-reaction");

  // Tarjeta 2
  const knownSpeciesSelect = document.getElementById("known-species");
  const targetSpeciesSelect = document.getElementById("target-species");

  const knownType = document.getElementById("known-type");
  const knownMassBlock = document.getElementById("known-mass-block");
  const knownMolesBlock = document.getElementById("known-moles-block");
  const knownVolumeBlock = document.getElementById("known-volume-block");

  const knownMassInput = document.getElementById("known-mass");
  const knownMolesInput = document.getElementById("known-moles");

  const knownVolumeKind = document.getElementById("known-volume-kind");
  const knownSolutionFields = document.getElementById("known-solution-fields");
  const knownGasFields = document.getElementById("known-gas-fields");

  const knownSolVolML = document.getElementById("known-sol-vol-ml");
  const knownSolM = document.getElementById("known-sol-molarity");

  const knownGasV = document.getElementById("known-gas-vol-L");
  const knownGasP = document.getElementById("known-gas-P");
  const knownGasPunit = document.getElementById("known-gas-Punit");
  const knownGasTc = document.getElementById("known-gas-Tc");

  const targetType = document.getElementById("target-type");
  const targetVolumeBlock = document.getElementById("target-volume-block");
  const targetVolumeKind = document.getElementById("target-volume-kind");
  const targetSolutionFields = document.getElementById("target-solution-fields");
  const targetGasFields = document.getElementById("target-gas-fields");

  const targetSolM = document.getElementById("target-sol-molarity");
  const targetSolVunit = document.getElementById("target-sol-vunit");

  const targetGasP = document.getElementById("target-gas-P");
  const targetGasPunit = document.getElementById("target-gas-Punit");
  const targetGasTc = document.getElementById("target-gas-Tc");

  const btnCalc = document.getElementById("btn-calc");
  const calcError = document.getElementById("calc-error");
  const resultBox = document.getElementById("result-box");
  const resultText = document.getElementById("result-text");

  const btnShowConversion = document.getElementById("btn-show-conversion");
  const conversionBox = document.getElementById("conversion-box");
  const conversionText = document.getElementById("conversion-text");

  const resetMass = document.getElementById("reset-mass");

  let currentSpecies = [];

  // ===============================
  //  UI adaptativa (disolución/gas)
  // ===============================
  function hide(el) { if (el) el.classList.add("hidden"); }
  function show(el) { if (el) el.classList.remove("hidden"); }

  function updateKnownUI() {
    hide(knownMassBlock);
    hide(knownMolesBlock);
    hide(knownVolumeBlock);

    if (knownType.value === "mass") show(knownMassBlock);
    if (knownType.value === "moles") show(knownMolesBlock);
    if (knownType.value === "volume") {
      show(knownVolumeBlock);
      updateKnownVolumeUI();
    }
  }

  function updateKnownVolumeUI() {
    if (knownType.value !== "volume") return;
    hide(knownSolutionFields);
    hide(knownGasFields);
    if (knownVolumeKind.value === "solution") show(knownSolutionFields);
    if (knownVolumeKind.value === "gas") show(knownGasFields);
  }

  function updateTargetUI() {
    hide(targetVolumeBlock);
    if (targetType.value === "volume") {
      show(targetVolumeBlock);
      updateTargetVolumeUI();
    }
  }

  function updateTargetVolumeUI() {
    if (targetType.value !== "volume") return;
    hide(targetSolutionFields);
    hide(targetGasFields);
    if (targetVolumeKind.value === "solution") show(targetSolutionFields);
    if (targetVolumeKind.value === "gas") show(targetGasFields);
  }

  knownType.addEventListener("change", updateKnownUI);
  knownVolumeKind.addEventListener("change", updateKnownVolumeUI);
  targetType.addEventListener("change", updateTargetUI);
  targetVolumeKind.addEventListener("change", updateTargetVolumeUI);

  updateKnownUI();
  updateTargetUI();

  // ===============================
  //  Añadir / quitar filas (reactivos/productos)
  // ===============================
  const MAX_PER_SIDE = 10; // límite suave para UI

  function attachRemoveHandler(btn) {
    btn.addEventListener("click", (e) => {
      const row = e.currentTarget.closest(".rowline");
      if (!row) return;
      const parent = row.parentElement;
      if (!parent) return;

      // no dejar una columna sin ninguna fila
      const rows = parent.querySelectorAll(".rowline");
      if (rows.length <= 1) {
        const input = row.querySelector("input");
        if (input) input.value = "";
        return;
      }
      row.remove();
    });
  }

  // Para filas iniciales del HTML
  document.querySelectorAll(".remove-species").forEach(attachRemoveHandler);

  function createRow(side) {
    const row = document.createElement("div");
    row.className = "rowline";

    const inp = document.createElement("input");
    inp.type = "text";
    inp.className = "species-input";
    inp.dataset.side = side;
    inp.placeholder = side === "reactivo"
      ? "Ej: H2, Fe2(SO4)3..."
      : "Ej: H2O, CO2...";

    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "icon-btn remove-species";
    btn.title = "Quitar";
    btn.setAttribute("aria-label", "Quitar");
    btn.textContent = "🗑️";
    attachRemoveHandler(btn);

    row.appendChild(inp);
    row.appendChild(btn);
    return row;
  }

  function addRow(side) {
    const parent = side === "reactivo" ? reactantsList : productsList;
    const count = parent.querySelectorAll(".rowline").length;
    if (count >= MAX_PER_SIDE) return;
    parent.appendChild(createRow(side));
  }

  addReactantBtn.addEventListener("click", () => addRow("reactivo"));
  addProductBtn.addEventListener("click", () => addRow("producto"));

  // ===============================
  //  Leer especies (N)
  // ===============================
  function readSpeciesFromUI() {
    const inputs = Array.from(document.querySelectorAll(".species-input"));
    const reactivos = [];
    const productos = [];

    for (const inp of inputs) {
      const side = inp.dataset.side;
      const val = normalizeFormula(inp.value.trim());
      if (!val) continue;
      if (side === "reactivo") reactivos.push(val);
      if (side === "producto") productos.push(val);
    }

    const species = [];
    reactivos.forEach((f, idx) => species.push({ id: `r${idx + 1}`, formula: f, side: "reactivo" }));
    productos.forEach((f, idx) => species.push({ id: `p${idx + 1}`, formula: f, side: "producto" }));

    return { species, reactivosCount: reactivos.length, productosCount: productos.length };
  }

  // ===============================
  //  Ajustar reacción
  // ===============================
  btnBalance.addEventListener("click", () => {
    balanceError.textContent = "";
    balancedInfo.classList.add("hidden");

    knownSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    targetSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    currentSpecies = [];

    const { species, reactivosCount, productosCount } = readSpeciesFromUI();

    if (reactivosCount < 1 || productosCount < 1) {
      balanceError.textContent = "Introduce al menos un reactivo y un producto.";
      return;
    }

    try {
      const A = buildMatrix(species);
      const coeffs = nullspaceVector(A);

      const enriched = species.map((sp, i) => ({
        ...sp,
        coeff: coeffs[i],
        molarMass: molarMassFromFormula(sp.formula)
      }));

      currentSpecies = enriched;

      const left = enriched
        .filter(s => s.side === "reactivo")
        .map(s => `${s.coeff === 1 ? "" : s.coeff + " "}${s.formula}`)
        .join(" + ");

      const right = enriched
        .filter(s => s.side === "producto")
        .map(s => `${s.coeff === 1 ? "" : s.coeff + " "}${s.formula}`)
        .join(" + ");

      balancedEquationText.textContent = left + " → " + right;

      molarMassList.innerHTML = "";
      enriched.forEach(s => {
        const li = document.createElement("li");
        li.textContent = `${s.formula}: ${fmt(s.molarMass, 3)} g/mol`;
        molarMassList.appendChild(li);
      });

      enriched.forEach(s => {
        const o1 = document.createElement("option");
        o1.value = s.id;
        o1.textContent = s.formula;
        knownSpeciesSelect.appendChild(o1);

        const o2 = document.createElement("option");
        o2.value = s.id;
        o2.textContent = s.formula;
        targetSpeciesSelect.appendChild(o2);
      });

      balancedInfo.classList.remove("hidden");
    } catch (err) {
      balanceError.textContent = err.message || "Error al ajustar la reacción.";
    }
  });

  // ===============================
  //  Calcular
  // ===============================
  btnCalc.addEventListener("click", () => {
    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";

    if (!currentSpecies.length) {
      calcError.textContent = "Primero ajusta la reacción.";
      return;
    }

    const knownId = knownSpeciesSelect.value;
    const targetId = targetSpeciesSelect.value;

    if (!knownId || !targetId) {
      calcError.textContent = "Selecciona sustancia conocida y sustancia a calcular.";
      return;
    }
    if (knownId === targetId) {
      calcError.textContent = "Las sustancias deben ser distintas.";
      return;
    }

    const known = currentSpecies.find(s => s.id === knownId);
    const target = currentSpecies.find(s => s.id === targetId);
    if (!known || !target) {
      calcError.textContent = "Error interno: selección inválida.";
      return;
    }

    let nKnown;
    const leftPieces = [];

    try {
      if (knownType.value === "mass") {
        const m = toNumber(knownMassInput.value);
        if (!isFinite(m) || m <= 0) throw new Error("Introduce una masa válida.");

        nKnown = m / known.molarMass;

        leftPieces.push(`<span class="factor">${termHTML(`${fmt(m, 3)}`, `g ${known.formula}`)}</span>`);
        leftPieces.push(`<span class="times">×</span>`);
        leftPieces.push(
          fractionHTML(`1`, `mol ${known.formula}`, `${fmt(known.molarMass, 3)}`, `g ${known.formula}`)
        );
      } else if (knownType.value === "moles") {
        const n = toNumber(knownMolesInput.value);
        if (!isFinite(n) || n <= 0) throw new Error("Introduce moles válidos.");

        nKnown = n;

        leftPieces.push(`<span class="factor">${termHTML(`${fmt(n, 6)}`, `mol ${known.formula}`)}</span>`);
      } else {
        // volumen
        if (knownVolumeKind.value === "solution") {
          const vml = toNumber(knownSolVolML.value);
          const M = toNumber(knownSolM.value);
          if (!isFinite(vml) || vml <= 0) throw new Error("Volumen (mL) inválido.");
          if (!isFinite(M) || M <= 0) throw new Error("Concentración (mol/L) inválida.");

          const vL = vml / 1000;
          nKnown = M * vL;

          leftPieces.push(`<span class="factor">${termHTML(`${fmt(M, 6)}`, `mol/L`)}</span>`);
          leftPieces.push(`<span class="times">×</span>`);
          leftPieces.push(`<span class="factor">${termHTML(`${fmt(vL, 6)}`, `L`)}</span>`);
          leftPieces.push(`<span class="times eq">=</span>`);
          leftPieces.push(`<span class="factor">${termHTML(`${fmt(nKnown, 6)}`, `mol ${known.formula}`)}</span>`);
        } else {
          const V = toNumber(knownGasV.value);
          const Praw = toNumber(knownGasP.value);
          const Tc = toNumber(knownGasTc.value);

          if (!isFinite(V) || V <= 0) throw new Error("Volumen (L) inválido.");
          if (!isFinite(Praw) || Praw <= 0) throw new Error("Presión inválida.");
          if (!isFinite(Tc)) throw new Error("Temperatura (°C) inválida.");

          const P_atm = pressureToAtm(Praw, knownGasPunit.value);
          const T = celsiusToK(Tc);
          if (T <= 0) throw new Error("Temperatura inválida.");

          nKnown = (P_atm * V) / (R * T);

          leftPieces.push(`<span class="factor">${termHTML(`${fmt(Praw, 3)}`, `${knownGasPunit.value}`)}</span>`);
          leftPieces.push(`<span class="times">×</span>`);
          leftPieces.push(`<span class="factor">${termHTML(`${fmt(V, 3)}`, `L`)}</span>`);
          leftPieces.push(`<span class="times">÷</span>`);
          leftPieces.push(`<span class="factor">${termHTML(`${fmt(R, 6)}`, `L·atm/(mol·K)`)}</span>`);
          leftPieces.push(`<span class="times">×</span>`);
          leftPieces.push(`<span class="factor">${termHTML(`${fmt(T, 2)}`, `K`)}</span>`);
          leftPieces.push(`<span class="times eq">=</span>`);
          leftPieces.push(`<span class="factor">${termHTML(`${fmt(nKnown, 6)}`, `mol ${known.formula}`)}</span>`);
        }
      }
    } catch (e) {
      calcError.textContent = e.message || "Dato conocido inválido.";
      return;
    }

    // estequiometría
    const nTarget = nKnown * (target.coeff / known.coeff);

    // salida
    let resultStr = "";
    let rightHTML = "";
    const tailPieces = [];

    try {
      if (targetType.value === "moles") {
        resultStr = `${fmt(nTarget, 6)} mol de ${target.formula}`;
        rightHTML = `<span class="factor result">${termHTML(`${fmt(nTarget, 6)}`, `mol ${target.formula}`)}</span>`;
      } else if (targetType.value === "mass") {
        const mTarget = nTarget * target.molarMass;
        resultStr = `${fmt(mTarget, 3)} g de ${target.formula}`;
        rightHTML = `<span class="factor result">${termHTML(`${fmt(mTarget, 3)}`, `g ${target.formula}`)}</span>`;

        tailPieces.push(`<span class="times">×</span>`);
        tailPieces.push(
          fractionHTML(`${fmt(target.molarMass, 3)}`, `g ${target.formula}`, `1`, `mol ${target.formula}`)
        );
      } else {
        // volumen objetivo
        if (targetVolumeKind.value === "solution") {
          const M = toNumber(targetSolM.value);
          if (!isFinite(M) || M <= 0) throw new Error("Concentración objetivo inválida.");

          const V_L = nTarget / M;
          const unit = targetSolVunit.value;
          const outV = unit === "mL" ? V_L * 1000 : V_L;

          resultStr = `${fmt(outV, 3)} ${unit} de disolución de ${target.formula}`;
          rightHTML = `<span class="factor result">${termHTML(`${fmt(outV, 3)}`, `${unit}`)}</span>`;
        } else {
          const Praw = toNumber(targetGasP.value);
          const Tc = toNumber(targetGasTc.value);

          if (!isFinite(Praw) || Praw <= 0) throw new Error("Presión objetivo inválida.");
          if (!isFinite(Tc)) throw new Error("Temperatura objetivo inválida.");

          const P_atm = pressureToAtm(Praw, targetGasPunit.value);
          const T = celsiusToK(Tc);
          if (T <= 0) throw new Error("Temperatura objetivo inválida.");

          const V_L = (nTarget * R * T) / P_atm;
          resultStr = `${fmt(V_L, 3)} L de ${target.formula} (gas ideal)`;
          rightHTML = `<span class="factor result">${termHTML(`${fmt(V_L, 3)}`, `L`)}</span>`;
        }
      }
    } catch (e) {
      calcError.textContent = e.message || "Datos de salida inválidos.";
      return;
    }

    // Mostrar resultado (solo resultado)
    resultText.textContent = resultStr;
    resultBox.classList.remove("hidden");

    // Construir conversión completa + tachado
    const stoich = fractionHTML(
      `${target.coeff}`, `mol ${target.formula}`,
      `${known.coeff}`, `mol ${known.formula}`
    );

    conversionText.innerHTML = `
      <span class="conversion-line">
        <span class="conv-left">
          ${leftPieces.join("")}
          <span class="times">×</span>
          ${stoich}
          ${tailPieces.join("")}
        </span>
        <span class="times eq">=</span>
        <span class="conv-right">${rightHTML}</span>
      </span>
    `;

    cancelUnitsIn(conversionBox);

    btnShowConversion.textContent = "Mostrar conversión";
    btnShowConversion.classList.remove("hidden");
  });

  // ===============================
  //  Mostrar / ocultar conversión
  // ===============================
  btnShowConversion.addEventListener("click", () => {
    const hidden = conversionBox.classList.contains("hidden");
    if (hidden) {
      conversionBox.classList.remove("hidden");
      btnShowConversion.textContent = "Ocultar conversión";
    } else {
      conversionBox.classList.add("hidden");
      btnShowConversion.textContent = "Mostrar conversión";
    }
  });

  // ===============================
  //  Reinicio tarjeta 1
  // ===============================
  resetReaction.addEventListener("click", () => {
    function resetList(listEl, side) {
      while (listEl.querySelectorAll(".rowline").length > 2) listEl.lastElementChild.remove();
      while (listEl.querySelectorAll(".rowline").length < 2) listEl.appendChild(createRow(side));
      listEl.querySelectorAll("input.species-input").forEach(i => (i.value = ""));
      listEl.querySelectorAll(".remove-species").forEach(btn => {
        if (!btn.dataset.bound) {
          attachRemoveHandler(btn);
          btn.dataset.bound = "1";
        }
      });
    }

    resetList(reactantsList, "reactivo");
    resetList(productsList, "producto");

    balanceError.textContent = "";
    balancedInfo.classList.add("hidden");

    knownSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    targetSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';

    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";

    currentSpecies = [];
  });

  // ===============================
  //  Reinicio tarjeta 2
  // ===============================
  resetMass.addEventListener("click", () => {
    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";

    knownSpeciesSelect.value = "";
    targetSpeciesSelect.value = "";

    knownMassInput.value = "";
    knownMolesInput.value = "";

    knownSolVolML.value = "";
    knownSolM.value = "";

    knownGasV.value = "";
    knownGasP.value = "";
    knownGasTc.value = "";

    targetSolM.value = "";
    targetGasP.value = "";
    targetGasTc.value = "";
  });

  // ===============================
  //  (IMPORTANTE) asegurar handlers
  // ===============================
  document.querySelectorAll(".remove-species").forEach(btn => {
    if (!btn.dataset.bound) {
      attachRemoveHandler(btn);
      btn.dataset.bound = "1";
    }
  });
});
