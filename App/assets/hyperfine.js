(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.ORBITAL_HYPERFINE = api;
})(globalThis, function () {
  "use strict";

  const HARTREE_PER_MEV = 1 / 27.211386245988;
  const parityValue = value => value === "+" ? 1 : value === "-" ? -1 : 0;
  const finiteNumber = (value, label) => {
    if (typeof value !== "number" || !Number.isFinite(value)) {
      throw new Error(`${label} must be a finite JSON number.`);
    }
    return value;
  };

  function twoJ(value, label) {
    const text = String(value).trim();
    const fraction = text.match(/^(\d+)\/2$/);
    if (fraction) {
      const result = Number(fraction[1]);
      if (!Number.isFinite(result) || result > 64) throw new Error(`${label} exceeds the current angular-momentum limit J ≤ 32.`);
      return result;
    }
    const decimal = text.match(/^(\d+)(?:\.(0|5))?$/);
    if (!decimal) throw new Error(`${label} must be a non-negative integer or half-integer.`);
    const result = 2 * Number(decimal[1]) + (decimal[2] === "5" ? 1 : 0);
    if (result > 64) throw new Error(`${label} exceeds the current angular-momentum limit J ≤ 32.`);
    return result;
  }

  function factorial(value) {
    if (!Number.isInteger(value) || value < 0) return 0;
    let result = 1;
    for (let n = 2; n <= value; n++) result *= n;
    return result;
  }

  function phase(exponent) {
    if (!Number.isInteger(exponent)) throw new Error("Angular-momentum phase exponent is not an integer.");
    return Math.abs(exponent) % 2 ? -1 : 1;
  }

  function wigner3j(j1x, j2x, j3x, m1x, m2x, m3x) {
    if (m1x + m2x + m3x !== 0) return 0;
    if ([m1x, m2x, m3x].some((m, i) => Math.abs(m) > [j1x, j2x, j3x][i] || ([j1x, j2x, j3x][i] - m) % 2 !== 0)) return 0;
    if (j3x > j1x + j2x || j3x < Math.abs(j1x - j2x) || (j1x + j2x + j3x) % 2 !== 0) return 0;

    const j1 = j1x / 2, j2 = j2x / 2, j3 = j3x / 2;
    const m1 = m1x / 2, m2 = m2x / 2, m3 = m3x / 2;
    const args = [
      j1 + j2 - j3, j1 - j2 + j3, -j1 + j2 + j3, j1 + j2 + j3 + 1,
      j1 + m1, j1 - m1, j2 + m2, j2 - m2, j3 + m3, j3 - m3,
    ];
    if (args.some(value => !Number.isInteger(value) || value < 0)) return 0;
    const triangle = Math.sqrt(
      factorial(args[0]) * factorial(args[1]) * factorial(args[2]) / factorial(args[3])
    );
    const magnetic = Math.sqrt(args.slice(4).reduce((product, value) => product * factorial(value), 1));
    const low = Math.max(0, Math.ceil(Math.max(j2 - j3 - m1, j1 - j3 + m2)));
    const high = Math.floor(Math.min(j1 + j2 - j3, j1 - m1, j2 + m2));
    let sum = 0;
    for (let z = low; z <= high; z++) {
      const terms = [
        z, j1 + j2 - j3 - z, j1 - m1 - z, j2 + m2 - z,
        j3 - j2 + m1 + z, j3 - j1 - m2 + z,
      ];
      if (terms.some(value => !Number.isInteger(value) || value < 0)) continue;
      const denominator = terms.reduce((product, value) => product * factorial(value), 1);
      sum += phase(z) / denominator;
    }
    return phase((j1x - j2x - m3x) / 2) * triangle * magnetic * sum;
  }

  function clebschGordan(j1, m1, j2, m2, totalJ, totalM) {
    if (m1 + m2 !== totalM) return 0;
    return phase((j1 - j2 + totalM) / 2)
      * Math.sqrt(totalJ + 1)
      * wigner3j(j1, j2, totalJ, m1, m2, -totalM);
  }

  function reducedElement(jx, mx, rank, q, jpx, mpx, reduced) {
    return phase((jx - mx) / 2)
      * wigner3j(jx, 2 * rank, jpx, -mx, 2 * q, mpx)
      * reduced;
  }

  function readStates(states, label, energyKey, angularKey, coefficientKey) {
    if (!Array.isArray(states) || states.length === 0) throw new Error(`${label} must contain at least one state.`);
    const ids = new Set();
    return states.map((state, index) => {
      if (!state || typeof state.id !== "string" || !state.id.trim()) {
        throw new Error(`${label}[${index}] requires a non-empty id.`);
      }
      if (ids.has(state.id)) throw new Error(`${label} contains duplicate state id '${state.id}'.`);
      ids.add(state.id);
      const j = twoJ(state[angularKey], `${label}[${index}] angular momentum`);
      const p = parityValue(state.parity);
      if (!p) throw new Error(`${label}[${index}] parity must be '+' or '-'.`);
      const energy = finiteNumber(state[energyKey], `${label}[${index}].${energyKey}`);
      const level = state.k;
      if (!Number.isInteger(level) || level < 1) throw new Error(`${label}[${index}] requires a positive integer k.`);
      const T = state.T;
      if (label === "nuclear_states" && (typeof T !== "number" || !Number.isFinite(T) || T < 0)) {
        throw new Error(`${label}[${index}] requires a non-negative isospin T.`);
      }
      const coefficients = state[coefficientKey];
      if (!Array.isArray(coefficients) || coefficients.length === 0) {
        throw new Error(`${label}[${index}] requires ${coefficientKey} wavefunction coefficients.`);
      }
      const coefficientIds = new Set();
      const normalized = coefficients.map((item, coefficientIndex) => {
        const id = item?.configuration ?? item?.basis_state;
        const coefficient = item?.coefficient;
        if (typeof id !== "string" || !id.trim() || typeof coefficient !== "number" || !Number.isFinite(coefficient)) {
          throw new Error(`${label}[${index}].${coefficientKey}[${coefficientIndex}] requires a configuration id and finite coefficient.`);
        }
        if (coefficientIds.has(id)) throw new Error(`${label}[${index}] has duplicate wavefunction basis id '${id}'.`);
        coefficientIds.add(id);
        return {configuration:id, coefficient};
      });
      const norm = normalized.reduce((sum, item) => sum + item.coefficient ** 2, 0);
      if (Math.abs(norm - 1) > 1e-6) {
        throw new Error(`${label}[${index}] ${coefficientKey} are not normalized (sum of squared coefficients is ${norm}).`);
      }
      return {...state, twoJ:j, parityValue:p, energy, basis_coefficients:normalized};
    });
  }

  function transitionTable(rows, states, label, rank, signature) {
    if (!Array.isArray(rows)) throw new Error(`${label} reduced matrix elements must be an array.`);
    const stateById = new Map(states.map(state => [state.id, state]));
    const values = new Map();
    for (const [index, row] of rows.entries()) {
      if (!row || !stateById.has(row.bra) || !stateById.has(row.ket)) {
        throw new Error(`${label}[${index}] names an unknown bra or ket state.`);
      }
      const value = finiteNumber(row.value, `${label}[${index}].value`);
      const bra = stateById.get(row.bra), ket = stateById.get(row.ket);
      if (Math.abs(bra.twoJ - ket.twoJ) > 2 * rank || bra.twoJ + ket.twoJ < 2 * rank) {
        throw new Error(`${label}[${index}] violates the rank-${rank} angular-momentum triangle rule.`);
      }
      if (bra.parityValue * ket.parityValue !== signature) {
        throw new Error(`${label}[${index}] violates the declared multipole parity.`);
      }
      const key = `${row.bra}\u0000${row.ket}`;
      if (values.has(key)) throw new Error(`${label} has duplicate reduced matrix element ${row.bra} ← ${row.ket}.`);
      values.set(key, value);
    }
    return (bra, ket) => {
      const key = `${bra.id}\u0000${ket.id}`;
      if (values.has(key)) return values.get(key);
      const reverse = `${ket.id}\u0000${bra.id}`;
      if (values.has(reverse)) return phase((bra.twoJ - ket.twoJ) / 2) * values.get(reverse);
      const allowed = Math.abs(bra.twoJ - ket.twoJ) <= 2 * rank
        && bra.twoJ + ket.twoJ >= 2 * rank
        && bra.parityValue * ket.parityValue === signature;
      if (allowed) throw new Error(`${label} is incomplete: missing signed rank-${rank} element ${bra.id} ← ${ket.id}. Supply explicit zeroes too.`);
      return 0;
    };
  }

  function coupledTensorElement(eBra, nBra, eKet, nKet, F, M, rank, eReduced, nReduced) {
    let total = 0;
    for (let me = -eBra.twoJ; me <= eBra.twoJ; me += 2) {
      const mn = 2 * M - me;
      if (Math.abs(mn) > nBra.twoJ || (nBra.twoJ - mn) % 2 !== 0) continue;
      const cgBra = clebschGordan(eBra.twoJ, me, nBra.twoJ, mn, 2 * F, 2 * M);
      if (!cgBra) continue;
      for (let mePrime = -eKet.twoJ; mePrime <= eKet.twoJ; mePrime += 2) {
        const qx = me - mePrime;
        if (qx % 2 !== 0 || Math.abs(qx) > 2 * rank) continue;
        const q = qx / 2;
        const mnPrime = 2 * M - mePrime;
        if (Math.abs(mnPrime) > nKet.twoJ || (nKet.twoJ - mnPrime) % 2 !== 0) continue;
        const cgKet = clebschGordan(eKet.twoJ, mePrime, nKet.twoJ, mnPrime, 2 * F, 2 * M);
        if (!cgKet) continue;
        const eValue = reducedElement(eBra.twoJ, me, rank, q, eKet.twoJ, mePrime, eReduced);
        const nValue = reducedElement(nBra.twoJ, mn, rank, -q, nKet.twoJ, mnPrime, nReduced);
        total += cgBra * cgKet * phase(q) * eValue * nValue;
      }
    }
    return total;
  }

  function diagonalize(matrix) {
    const n = matrix.length;
    if (!n || n > 160) throw new Error("Coupled basis dimension must be between 1 and 160 for the in-browser dense solver.");
    const a = matrix.map(row => row.slice());
    const vectors = Array.from({length:n}, (_, row) => Array.from({length:n}, (_, col) => row === col ? 1 : 0));
    for (let sweep = 0; sweep < 80; sweep++) {
      let largest = 0;
      for (let p = 0; p < n - 1; p++) {
        for (let q = p + 1; q < n; q++) {
          const apq = a[p][q];
          largest = Math.max(largest, Math.abs(apq));
          if (Math.abs(apq) <= 1e-14 * (1 + Math.abs(a[p][p]) + Math.abs(a[q][q]))) continue;
          const tau = (a[q][q] - a[p][p]) / (2 * apq);
          const t = (tau >= 0 ? 1 : -1) / (Math.abs(tau) + Math.sqrt(1 + tau * tau));
          const c = 1 / Math.sqrt(1 + t * t), s = t * c;
          const app = a[p][p], aqq = a[q][q];
          a[p][p] = app - t * apq;
          a[q][q] = aqq + t * apq;
          a[p][q] = a[q][p] = 0;
          for (let r = 0; r < n; r++) {
            if (r !== p && r !== q) {
              const arp = a[r][p], arq = a[r][q];
              a[r][p] = a[p][r] = c * arp - s * arq;
              a[r][q] = a[q][r] = s * arp + c * arq;
            }
            const vrp = vectors[r][p], vrq = vectors[r][q];
            vectors[r][p] = c * vrp - s * vrq;
            vectors[r][q] = s * vrp + c * vrq;
          }
        }
      }
      if (largest < 1e-13) break;
      if (sweep === 79) throw new Error("Hyperfine diagonalization did not converge.");
    }
    return Array.from({length:n}, (_, col) => ({
      energy: a[col][col],
      vector: vectors.map(row => row[col]),
    })).sort((left, right) => left.energy - right.energy);
  }

  function solveSpinCoupledProduct(input) {
    const electronic = input?.electronic_state;
    const nuclear = input?.nuclear_level;
    if (!electronic || !nuclear) throw new Error("An electronic GRASP state and an evaluated nuclear level are required.");
    if (typeof electronic.id !== "string" || !electronic.id.trim()) throw new Error("The electronic state requires an id.");
    const electronicJ = twoJ(electronic.J, "Electronic J");
    const nuclearJ = twoJ(nuclear.spin, "Evaluated nuclear spin");
    const electronicParity = parityValue(electronic.parity);
    const nuclearParity = parityValue(nuclear.parity);
    if (!electronicParity) throw new Error("Electronic parity must be '+' or '-'.");
    if (!nuclearParity || nuclear.spin_parity_assignment_exact !== true) {
      throw new Error("This ENSDF level has no exact spin/parity assignment and cannot be selected as a definite spin state.");
    }
    if (typeof nuclear.id !== "string" || !nuclear.id.trim()) throw new Error("The evaluated nuclear level requires an id.");
    const excitationMeV = finiteNumber(nuclear.energy_mev, "ENSDF excitation energy");
    if (excitationMeV < 0) throw new Error("ENSDF excitation energy cannot be negative.");
    const totalF = String(input.F ?? "").trim();
    const totalF2 = twoJ(totalF, "Total F");
    if (totalF2 < Math.abs(electronicJ - nuclearJ)
        || totalF2 > electronicJ + nuclearJ
        || (electronicJ + nuclearJ + totalF2) % 2 !== 0) {
      throw new Error("Requested F violates the angular-momentum triangle or integer-step rule for this electronic state and ENSDF level.");
    }
    const expansion = [];
    for (let mJ = -electronicJ; mJ <= electronicJ; mJ += 2) {
      const mI = totalF2 - mJ;
      if (Math.abs(mI) > nuclearJ || (nuclearJ - mI) % 2 !== 0) continue;
      const coefficient = clebschGordan(electronicJ, mJ, nuclearJ, mI, totalF2, totalF2);
      if (Math.abs(coefficient) < 1e-14) continue;
      expansion.push({
        electronic_state:electronic.id,
        nuclear_state:nuclear.id,
        mJ:mJ / 2,
        mI:mI / 2,
        coefficient,
        probability:coefficient ** 2,
      });
    }
    const norm = expansion.reduce((sum, term) => sum + term.coefficient ** 2, 0);
    if (Math.abs(norm - 1) > 1e-10) throw new Error(`Spin-coupled product expansion is not normalized (norm squared ${norm}).`);
    return {
      F:totalF,
      M:totalF,
      parity:electronicParity * nuclearParity > 0 ? "+" : "-",
      energy_offset_mev:excitationMeV,
      energy_reference:"Nuclear excitation energy relative to the selected electronic GRASP state; no hyperfine shift or absolute electronic energy is included.",
      dimension:expansion.length,
      total_wavefunction_expansion:expansion,
      electronic_basis_resolved:false,
      nuclear_basis_resolved:false,
      caveat:"This is the Clebsch–Gordan-coupled product of one selected GRASP electronic level and one ENSDF nuclear spin/parity level. It does not resolve internal nuclear spatial structure, calculate hyperfine shifts, or mix electronic/nuclear levels.",
    };
  }

  function solve(input) {
    if (!input || input.schema !== "orbital-hyperfine-v1") {
      throw new Error("Expected input schema 'orbital-hyperfine-v1'.");
    }
    const electronic = readStates(input.electronic_states, "electronic_states", "energy_hartree", "J", "basis_coefficients");
    const nuclear = readStates(input.nuclear_states, "nuclear_states", "energy_mev", "I", "basis_coefficients");
    const targetF2 = twoJ(input.target?.F, "Target total F");
    const targetJ = twoJ(input.target?.electronic_J, "Target electronic J");
    const targetElectronicParity = parityValue(input.target?.electronic_parity);
    const targetElectronicK = Number(input.target?.electronic_k);
    if (!targetElectronicParity || !Number.isInteger(targetElectronicK) || targetElectronicK < 1
        || typeof input.target?.nuclear_state_id !== "string") {
      throw new Error("Target requires electronic Jπ,k and a selected nuclear state id.");
    }
    const targetElectronicMatches = electronic.filter(state => state.twoJ === targetJ
      && state.parityValue === targetElectronicParity && state.k === targetElectronicK);
    if (targetElectronicMatches.length !== 1) {
      throw new Error(targetElectronicMatches.length
        ? "Electronic state data ambiguously matches the requested Jπ,k; give each state a unique k within its symmetry."
        : "No supplied electronic wavefunction matches the requested Jπ and k.");
    }
    const targetElectronic = targetElectronicMatches[0];
    const targetNuclear = nuclear.find(state => state.id === input.target.nuclear_state_id);
    if (!targetNuclear) throw new Error(`Unknown target nuclear state '${input.target.nuclear_state_id}'.`);
    const targetParity = targetElectronic.parityValue * targetNuclear.parityValue;
    if (!Array.isArray(input.interactions) || input.interactions.length === 0) {
      throw new Error("At least one sourced hyperfine interaction multipole is required.");
    }

    const interactions = input.interactions.map((interaction, index) => {
      const rank = interaction.rank, signature = parityValue(interaction.parity);
      const scale = interaction.scale_hartree;
      if (!Number.isInteger(rank) || rank < 0 || rank > 12 || !signature || !Number.isFinite(scale)) {
        throw new Error(`interactions[${index}] requires integer rank, parity, and finite scale_hartree.`);
      }
      if (typeof interaction.source !== "string" || !interaction.source.trim()) {
        throw new Error(`interactions[${index}] requires provenance (for example GRASP rtransition and KSHELL).`);
      }
      return {
        rank, signature, scale, source:interaction.source,
        electronic:transitionTable(interaction.electronic_reduced, electronic, `interaction ${index} electronic`, rank, signature),
        nuclear:transitionTable(interaction.nuclear_reduced, nuclear, `interaction ${index} nuclear`, rank, signature),
      };
    });

    const basis = [];
    for (const e of electronic) {
      for (const n of nuclear) {
        if (Math.abs(e.twoJ - n.twoJ) <= targetF2 && targetF2 <= e.twoJ + n.twoJ
            && e.parityValue * n.parityValue === targetParity) {
          basis.push({electronic:e, nuclear:n, label:`${e.id} ⊗ ${n.id}; F=${input.target.F}${input.target.parity}`});
        }
      }
    }
    if (!basis.length) throw new Error("No electronic–nuclear product basis states satisfy the requested F and total parity.");
    if (basis.length > 160) throw new Error(`The requested coupled basis has ${basis.length} states; the dense browser solver limit is 160.`);
    const targetBasisIndex = basis.findIndex(item =>
      item.electronic.id === targetElectronic.id && item.nuclear.id === targetNuclear.id
    );
    if (targetBasisIndex < 0) {
      throw new Error("The requested electronic and nuclear states cannot couple to the selected total F and parity.");
    }

    const matrix = basis.map((bra, i) => basis.map((ket, j) => {
      let value = i === j
        ? bra.electronic.energy + bra.nuclear.energy * HARTREE_PER_MEV
        : 0;
      for (const interaction of interactions) {
        const eReduced = interaction.electronic(bra.electronic, ket.electronic);
        const nReduced = interaction.nuclear(bra.nuclear, ket.nuclear);
        if (eReduced && nReduced) {
          value += interaction.scale * coupledTensorElement(
            bra.electronic, bra.nuclear, ket.electronic, ket.nuclear,
            targetF2 / 2, targetF2 / 2, interaction.rank, eReduced, nReduced
          );
        }
      }
      return value;
    }));

    let asymmetry = 0;
    for (let i = 0; i < basis.length; i++) {
      for (let j = 0; j < basis.length; j++) asymmetry = Math.max(asymmetry, Math.abs(matrix[i][j] - matrix[j][i]));
    }
    if (asymmetry > 1e-9) {
      throw new Error(`Interaction matrix is not Hermitian (maximum asymmetry ${asymmetry.toExponential(2)} Ha); check signed reduced elements and phase conventions.`);
    }

    const roots = diagonalize(matrix);
    const selectedRoot = roots.map((root, index) => ({
      ...root, index,
      targetOverlap:root.vector[targetBasisIndex] ** 2,
    })).sort((left, right) => right.targetOverlap - left.targetOverlap)[0];
    const selected = roots[selectedRoot.index];
    const maxComponent = Math.max(...selected.vector.map(Math.abs));
    const sign = selected.vector.find(value => Math.abs(value) === maxComponent) < 0 ? -1 : 1;
    const components = basis.map((entry, index) => ({
      electronic_state:entry.electronic.id,
      nuclear_state:entry.nuclear.id,
      F:input.target.F,
      parity:targetParity > 0 ? "+" : "-",
      coefficient:sign * selected.vector[index],
      probability: selected.vector[index] ** 2,
    })).filter(component => component.probability > 1e-12);
    const expansion = [];
    for (let basisIndex = 0; basisIndex < basis.length; basisIndex++) {
      const outer = sign * selected.vector[basisIndex];
      if (Math.abs(outer) < 1e-12) continue;
      const eState = basis[basisIndex].electronic;
      const nState = basis[basisIndex].nuclear;
      for (const eComponent of eState.basis_coefficients) {
        for (const nComponent of nState.basis_coefficients) {
          for (let me = -eState.twoJ; me <= eState.twoJ; me += 2) {
            const mn = targetF2 - me;
            if (Math.abs(mn) > nState.twoJ || (nState.twoJ - mn) % 2 !== 0) continue;
            const cg = clebschGordan(eState.twoJ, me, nState.twoJ, mn, targetF2, targetF2);
            if (Math.abs(cg) < 1e-14) continue;
            expansion.push({
              electronic_state:eState.id,
              electronic_configuration:eComponent.configuration ?? eComponent.basis_state,
              mJ:me / 2,
              nuclear_state:nState.id,
              nuclear_configuration:nComponent.configuration ?? nComponent.basis_state,
              mI:mn / 2,
              coefficient:outer * eComponent.coefficient * nComponent.coefficient * cg,
            });
            if (expansion.length > 100_000) {
              throw new Error("Total wavefunction expansion exceeds the 100,000-term output limit.");
            }
          }
        }
      }
    }
    const expansionNorm = expansion.reduce((sum, term) => sum + term.coefficient ** 2, 0);
    if (Math.abs(expansionNorm - 1) > 1e-8) {
      throw new Error(`Expanded total wavefunction is not normalized (norm squared ${expansionNorm}). Check state coefficients and angular-momentum conventions.`);
    }
    return {
      energy_hartree:selected.energy,
      energy_mev:selected.energy / HARTREE_PER_MEV,
      energy_reference:"Electronic energies are relative to the supplied GRASP zero and nuclear energies are excitation energies relative to the supplied nuclear reference; nuclear rest-mass energy is excluded.",
      F:input.target.F,
      parity:targetParity > 0 ? "+" : "-",
      hyperfine_root_k:selectedRoot.index + 1,
      target_overlap:selectedRoot.targetOverlap,
      target_identification_warning:selectedRoot.targetOverlap < 0.5
        ? "The selected eigenstate has less than 50% overlap with the requested uncoupled reference; inspect nearby roots and basis convergence."
        : null,
      target_electronic_state:targetElectronic.id,
      target_nuclear_state:targetNuclear.id,
      dimension:basis.length,
      components,
      total_wavefunction_expansion:expansion,
      provenance:interactions.map(item => item.source),
      caveat:"Finite-basis hyperfine eigenvector in supplied GRASP electronic and KSHELL/evaluated nuclear state bases. The accuracy is limited by input states, signed reduced matrix elements, interaction multipoles, and phase conventions.",
    };
  }

  return {solve, solveSpinCoupledProduct, twoJ};
});
