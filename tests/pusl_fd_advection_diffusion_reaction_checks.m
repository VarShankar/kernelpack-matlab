function pusl_fd_advection_diffusion_reaction_checks()
%PUSL_FD_ADVECTION_DIFFUSION_REACTION_CHECKS Basic checks for PU-SL + FD diffusion + reaction.

common = pusl_advection_diffusion_reaction_common();
domain = common.create_disk_domain(0.16);

solver = kp.solvers.PUSLFDAdvectionDiffusionReactionSolver();
solver.init(domain, 4, 4, 0.01, 0.05, "backward");
solver.enableHomogeneousNeumannMassConservation(true);

Xout = solver.getOutputNodes();
u0 = common.exact_solution(0.0, Xout);
solver.setInitialState(u0);

u1 = solver.bdf1Step(0.01, ...
    @(t, X) common.boundary_zero_swirl(domain, t, X), [], ...
    @(t, state, X) common.reaction(t, state, X), ...
    @(nuValue, t, X) common.forcing(domain, nuValue, t, X), ...
    @(t, X) common.unit_coeff(t, X), ...
    @(t, X) common.zero_coeff(t, X), ...
    @(NeuCoeffs, DirCoeffs, nr, t, X) common.homogeneous_neumann(NeuCoeffs, DirCoeffs, nr, t, X));

uTrue = common.exact_solution(0.01, Xout);
assert(norm(u1 - uTrue) / norm(uTrue) < 5e-1, ...
    'PU-SL + FD diffusion + reaction BDF1 should track the manufactured solution to modest accuracy.');

disp('pusl fd advection diffusion reaction checks passed');
end
