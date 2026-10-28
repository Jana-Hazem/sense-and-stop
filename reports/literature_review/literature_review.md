# Literature Review: CSD-Inspired Waves and Sparse-State Reconstruction

## Scope

This synthesis adapts the team's Neuroflux literature review for the current repository. It focuses on cortical spreading depolarization (CSD), continuum models, and reconstruction from sparse observations. It is a focused background review, not a systematic review. CSD is a proposed physiological substrate of typical migraine aura, but direct electrophysiological recording during a human migraine attack has not been established; human imaging findings are indirect [1]. Accordingly, this repository studies a computational wave-reconstruction framework and makes no clinical or diagnostic claim.

## Biological motivation and limits of the evidence

CSD involves a slowly propagating episode of near-complete neuronal and glial depolarization, large ionic shifts, suppression of normal firing, and a negative direct-current shift. Extracellular potassium accumulation and local cellular interactions contribute to propagation, while pumps and glial buffering support recovery. Metabolic demand, oxygen supply, and perfusion can affect the event and its recovery [2–4]. These features motivate studying a traveling wave with activation and recovery; they do not imply that a single scalar field captures all CSD physiology.

Much quantitative evidence comes from animal or ex vivo preparations. Human imaging is indirect, and invasive recordings from injured cortex are not equivalent to recordings during spontaneous migraine. Species, induction method, tissue condition, and measurement scale can change observed waveforms and vascular responses [1–4]. Experimental measurements can therefore inform modeling assumptions and parameter bounds, but they do not validate a patient-specific reconstruction.

## Model hierarchy and implications for this project

CSD models range from reduced reaction–diffusion equations to detailed ionic electrodiffusion, membrane-current, pump, glial, swelling, and neurovascular models [3–5]. Detailed models can address mechanisms but bring more states and uncertain parameters. A reduced model can be appropriate when the aim is to test reconstruction methods, provided its scope is stated clearly.

The supplied review identifies a limitation in the proposed scalar Fisher–KPP equation,

\[
u_t = D\,u_{xx} + r\,u(1-u).
\]

It is a convenient way to generate a traveling front [6, 7], but its usual monostable dynamics allow growth from arbitrarily small perturbations and settle at an elevated state. It has no explicit excitation threshold, refractory period, or recovery variable, and does not represent ionic gradients, pump-mediated homeostasis, glial buffering, oxygen, or perfusion. Matching its front speed to a reported CSD speed would not resolve these structural omissions.

The review therefore proposes treating Fisher–KPP as a front-reconstruction baseline and considering a two-variable excitable reaction–diffusion model with an activator and a recovery variable as the primary CSD-inspired test system. A slow-fast CSD model provides a relevant physiological-modeling precedent [5]. Such a recovery variable would be phenomenological, not a measured ionic or metabolic quantity. This is a literature-informed recommendation for the team to consider; it does not, by itself, change the model decision recorded elsewhere in the repository. A mechanistic question would require a physiologically grounded model and a larger parameter-identifiability effort [3–5].

## Inverse-problem framing

For observations of the form \(y = H[u] + \epsilon\), the target must be specified. The literature review distinguishes:

1. **Trajectory-state reconstruction:** recover the evolving field \(u(x,t)\) with model coefficients fixed. This is the proposed primary target.
2. **Initial-condition reconstruction:** recover \(u(x,0)\) while holding the model and coefficients fixed. Treat this as a separate experiment.
3. **Parameter estimation:** infer coefficients such as diffusivity or recovery rates. Defer this until identifiability has been assessed; jointly estimating state and parameters is harder.

Reconstruction depends on observability, sensor placement, spatial coverage, and observation timing, not sensor count alone. Noise level and distribution should be explicit. A regularized data-fit objective can stabilize recovery, but regularization introduces a bias–variance trade-off. Measurement noise should also be distinguished from model discrepancy, since the generating wave may not follow the inverse model [8–10].

Synthetic ground truth enables error measurement, but using effectively identical forward and inverse solvers can create an inverse crime and overstate performance. The reviewed work recommends generating data at finer resolution or with an independent solver, and, where feasible, testing model mismatch and realistic sensor averaging [11]. These are methodological recommendations from the literature; project evaluation and results remain reserved for their later phase.

## Research gap and question

The supplied review did not identify, among its examined sources, a CSD-focused study that jointly varies point-sensor number and placement, observation noise, and regularization while testing reconstruction against model or discretization mismatch. This is a bounded gap in the reviewed sources, not a claim that no related work exists. The contribution should be framed as transferring and assessing established inverse-problem methods in a CSD-motivated setting, rather than inventing sparse-wave reconstruction or proposing a clinical diagnostic method.

The review's proposed question is: for a specified two-variable excitable reaction–diffusion model that produces a traveling pulse with recovery, how do point-sensor density and placement, measurement noise, and regularization affect reconstruction error of the full spatiotemporal state, and how robust is recovery when the generating and inverse models differ in resolution or dynamics?

## References

1. O’Hare, L., Asher, J. M., & Hibbard, P. B. (2021). Migraine visual aura and cortical spreading depression—Linking mathematical models to empirical evidence. *Vision, 5*(2), 30. https://doi.org/10.3390/vision5020030
2. Takano, T., et al. (2007). Cortical spreading depression causes and coincides with tissue hypoxia. *Nature Neuroscience, 10*(6), 754–762. https://doi.org/10.1038/nn1902
3. Chang, J. C., et al. (2013). A mathematical model of the metabolic and perfusion effects on cortical spreading depression. *PLoS ONE, 8*(8), e70469. https://doi.org/10.1371/journal.pone.0070469
4. Ellingsrud, A. J., et al. (2022). Validating a computational framework for ionic electrodiffusion with cortical spreading depression as a case study. *eNeuro, 9*(2), ENEURO.0408-21.2022. https://doi.org/10.1523/ENEURO.0408-21.2022
5. Reyner-Parra, D., Bonet, C., Seara, T. M., & Huguet, G. (2023). Traveling waves in a model for cortical spreading depolarization with slow-fast dynamics. *Chaos, 33*, 083154. https://doi.org/10.1063/5.0160509
6. Fisher, R. A. (1937). The wave of advance of advantageous genes. *Annals of Eugenics, 7*, 355–369.
7. Kolmogorov, A. N., Petrovsky, I. G., & Piskunov, N. S. (1937). Study of the diffusion equation with growth of the quantity of matter and its application to a biology problem. *Bulletin of Moscow University, Mathematics and Mechanics, 1*(6), 1–25.
8. Badr, S., Fenton, F. H., & Cherry, E. M. (2023). Reconstructing cardiac voltage using data assimilation: Effects of observation distribution. *Computing in Cardiology, 50*, 407. https://doi.org/10.22489/CinC.2023.407
9. Lukyanenko, D., et al. (2021). Inverse problem of recovering the initial condition for a nonlinear equation of the reaction–diffusion–advection type by data given on the position of a reaction front with a time delay. *Mathematics, 9*(4), 342. https://doi.org/10.3390/math9040342
10. Ahmed, S. E., San, O., & Lakshmivarahan, S. (2023). On the dual advantage of placing observations through forward sensitivity analysis. *Proceedings of the Royal Society A, 479*, 20220815. https://doi.org/10.1098/rspa.2022.0815
11. Kaipio, J., & Somersalo, E. (2007). Statistical inverse problems: Discretization, model reduction and inverse crimes. *Journal of Computational and Applied Mathematics, 198*(2), 493–504. https://doi.org/10.1016/j.cam.2005.09.027
