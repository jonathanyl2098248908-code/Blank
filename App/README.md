# Orbital

Orbital is a local atomic-structure explorer. Open `index.html` in a modern browser or serve this directory as static files. Keep the adjacent `assets/` directory; it contains the bundled GRASP WebAssembly executables and third-party notices. GRASP and preview calculations run locally; ENSDF retrieval requires the Python server and an active network connection.

## Calculation modes

- **Browser GRASP MCDHF** runs the GRASP relativistic multiconfiguration Dirac–Hartree–Fock workflow in a WebAssembly worker. It can optionally run GRASP `rci`, and uses GRASP `rdensity` and `rwfntotxt` output to show state-resolved electron density and selected radial large-component signs. You can request electronic J, parity, and level index k; generated starters restrict the J range and reject a reference configuration with the wrong parity, and a browser result is rejected unless GRASP's density output contains the requested Jπ and vector index. Here k identifies a GRASP output vector, not a universal excitation label. **Generate GRASP inputs** uses the selected element, charge, isotope, starter-configuration preset, and requested Jπ to create `rcsf.inp`, `isodata`, and default run transcripts locally. These are starting guesses: one approximate Aufbau/promoted reference configuration, an empirical nuclear radius, mass A u, and zero nuclear moments. Starting orbitals use screened-hydrogenic estimates below Z=13 and Thomas–Fermi estimates at Z=13 and above. Light isotopes use a 0.5 fm nuclear skin to avoid invalid default Fermi parameters. Review and improve them before scientific use; they are not a validated or converged basis.
- **Fast approximate preview** runs a spherical, nonrelativistic Kohn–Sham LDA-SCF model on a logarithmic radial grid. It is an approximate model, not an alternative implementation of GRASP.
- **Native GRASP server** runs an installed GRASP build through a loopback-only Python server and provides native output files and logs.

## Coupled electronic–nuclear hyperfine states

The **Solve hyperfine state** panel diagonalizes a finite-basis coupled Hamiltonian at fixed total `F`, `M=F`, and total parity. It expands the selected eigenvector into products of GRASP electronic CSFs and nuclear basis configurations, with Clebsch–Gordan spin coupling. Nuclear excitation energies are entered in MeV; electronic state energies are entered in Hartree. The output eigenvalue combines the supplied electronic energy reference and nuclear excitation energy, and excludes nuclear rest-mass energy.

The **ENSDF spin-coupled product** route needs no `.snt` interaction file or user-entered matrix elements. It combines a successfully calculated, selected GRASP Jπk level with an evaluated ENSDF level, using Clebsch–Gordan coefficients to couple their J and I to requested F (M=F). The downloaded expansion is in electronic and nuclear magnetic substates. It does not resolve internal nuclear spatial structure, calculate hyperfine shifts, or mix electronic/nuclear levels; its energy offset is only the nuclear excitation relative to the selected electronic level.

The app can retrieve evaluated nuclear levels from the IAEA Nuclear Data Services LiveChart API when served by `server.py`. It requests the selected nuclide's ENSDF level records and displays energies, reported Jπ/isospin, and available static moments with source attribution. This is level-data retrieval, not a nuclear-structure calculation: ENSDF does not provide the many-body nuclear basis eigenvectors or generally complete signed transition matrix elements needed for the full coupled Hamiltonian.

The current hyperfine diagonalizer still accepts a user-prepared JSON dataset with schema `orbital-hyperfine-v1`; `assets/hyperfine-template.json` documents its shape but is marked template-only and cannot be solved. `metadata` must match selected Z/charge/A and identify `electronic_source: "GRASP"`, nuclear source (`"KSHELL"` or `"ENSDF"`), and the source/version for nuclear wavefunctions. KSHELL metadata must also name its interaction and model space; ENSDF metadata must identify the evaluated dataset. The selected electronic `J`, parity, and `k` must match exactly one electronic state; choose a nuclear state by its dataset ID. Every electronic and nuclear state must contain a normalized `basis_coefficients` array (`configuration` plus real coefficient), and nuclear states must include `I`, parity, isospin `T`, level `k`, and excitation energy. `interactions` must provide tensor rank, parity signature, an energy scale in Hartree, provenance, and signed reduced matrix-element tables for both sectors. The app uses the Condon–Shortley Wigner–Eckart convention

`<j m|T(k,q)|j' m'> = (-1)^(j-m) (j k j'; -m q m') <j||T(k)||j'>`

and contracts each multipole as `sum_q (-1)^q T_e(k,q) T_n(k,-q)`. `scale_hartree` converts the product of the supplied electronic and nuclear reduced-matrix-element units to Hartree. Supply every angular-momentum/parity-allowed transition, including explicit zeroes; missing allowed elements, non-normalized states, and non-Hermitian matrices block the calculation. The total expansion downloads as JSON. Its phase is fixed by making the largest basis amplitude positive; signs of source reduced elements still depend on the source phase convention.

The app does **not** yet launch KSHELL, run GRASP `rtransition`, or adapt their native outputs into the coupled-solver dataset. It still requires manually prepared wavefunction and signed matrix-element exports to run the finite-basis hyperfine solve; the ENSDF retrieval does not replace those inputs. KSHELL requires a model-space interaction, and results depend on that interaction and truncation—there is no unique nuclear wavefunction implied by Z/A alone. The dense browser diagonalizer is limited to 160 coupled product states. No coupled result is reported when the requested finite-basis hyperfine Hamiltonian is incomplete.

The GRASP density viewer remains an electronic spherical-average view, not a rendering of the coupled hyperfine wavefunction. The coupled solver computes a finite-basis atom–nucleus expansion only when the user supplies complete normalized state vectors and signed matrix elements. It is not FCI, RCC/CCSDT(Q), or a guaranteed DCB+QED calculation. Breit and QED corrections are only included if selected in the user-supplied GRASP workflow and confirmed in its output; this app does not enable or independently verify them. Neither mode guarantees mathematical exactness. Results depend on the input basis, nuclear model, convergence, state identification, interaction multipoles, and enabled corrections. The viewer ray-marches a full 3D volume from the supplied radial spherical-average density rather than displaying a planar slice; it cannot represent non-spherical density anisotropy. The visualized phase is only the sign of a selected real large radial component with an approximate angular factor, not the phase of a full spinor or many-electron wavefunction.

When using manually supplied inputs, element, isotope, charge, and excitation controls do not rewrite those files. Regenerate the starter files whenever the selected species or state changes, then review them. The isotope stability chart is illustrative, not experimental nuclear data. ENSDF retrieval requires opening the app through the local server and an active network connection. Browser progress reflects workflow stages; it is not a prediction of completion time. Calculations may take substantial time and memory for large CSF bases.

## Run the native GRASP backend locally

The backend uses only Python's standard library. Install GRASP using the [official instructions](https://compas.github.io/grasp/) (Fortran/C/C++ toolchain and CMake required). For example:

```sh
git clone https://github.com/compas/grasp.git
cd grasp
./configure.sh
cmake --build build --target install -j2
```

In another terminal, from this project directory, start the app server with the GRASP executable directory:

```sh
GRASP_BIN=/absolute/path/to/grasp/bin python3 server.py
```

Open `http://127.0.0.1:8765`, choose **Native GRASP server**, and submit the GRASP inputs. The server binds to loopback only. Set `ORBITAL_PORT` to change the port or `ORBITAL_JOB_DIR` to change where per-job files are stored. Use the server-served page for native requests; browsers restrict `file://` pages from calling local HTTP services.

The form accepts a GRASP CSF configuration (`rcsf.inp`), nuclear data (`isodata`), and interactive answer transcripts for `rangular`, `rwfnestimate`, and `rmcdhf`; an optional `rci` transcript invokes RCI. Current GRASP writes its wavefunction output directly, so the app does not use the obsolete `rsave` stage. Prepare inputs for the exact element, ion, isotope, symmetry, and state. Consult the GRASP manual and its test examples. Failed or unsupported calculations are reported, not replaced with preview data.

## Bundled browser engine provenance

The browser executables were built from [GRASP](https://github.com/compas/grasp) revision `9006157730a82ac839f2b4ff4e938bcba63a539e` with Flang/Emscripten using `ghcr.io/r-wasm/flang-wasm:main`. GRASP source compatibility adjustments for the WebAssembly runtime are provided as `assets/grasp-wasm-compat.patch`; bundled dependency licenses and notices are in `assets/THIRD_PARTY_NOTICES.txt`. The standalone browser engine bundle is `assets/grasp-wasm-sources.js`.
