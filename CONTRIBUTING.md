# Contributing

Bug reports, numerical reproductions, documentation corrections, and focused
pull requests are welcome. Please open an issue before beginning a substantial
API or algorithm change so the mathematical and compatibility requirements can
be agreed upon first.

Contributions should:

- preserve the existing `kp` package organization and MATLAB-native APIs;
- include a focused check in `tests/` for behavioral changes;
- avoid changing a published numerical algorithm without documenting and
  validating the mathematical difference; and
- pass `run_public_checks` and `mip_package_checks` before submission.

By contributing, you agree that your contribution is distributed under the
BSD 3-Clause License in this repository.
