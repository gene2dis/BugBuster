## Description

<!-- What does this PR change, and why? -->

## Related issues

<!-- e.g. Closes #12 -->

## How was this tested?

- [ ] `nextflow run main.nf -profile test,docker -stub`
- [ ] `nextflow run main.nf -profile test,docker`
- [ ] `nf-test test`

## Checklist

- [ ] The PR targets `dev` (not `main`)
- [ ] CI passes
- [ ] Commits follow the conventional commit format (`feat(scope): ...`, `fix: ...`)
- [ ] New processes follow the module conventions in [CONTRIBUTING.md](../CONTRIBUTING.md) (stub block, `versions.yml`, container pinned, `meta.yml`)
- [ ] New parameters are declared in `nextflow_schema.json` and documented in `docs/parameters.md` (a new database registry map goes in `defaultIgnoreParams`; see CONTRIBUTING.md, Adding a parameter)
- [ ] New code has a test (CONTRIBUTING.md, Adding a test)
- [ ] `CHANGELOG.md` updated for user-visible changes
