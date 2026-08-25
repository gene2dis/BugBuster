## Description

<!-- What does this PR change, and why? -->

## Related issues

<!-- e.g. Closes #12 -->

## How was this tested?

- [ ] `nextflow run main.nf -profile test,docker -stub`
- [ ] `nextflow run main.nf -profile test,docker`
- [ ] `nf-test test`

## Checklist

- [ ] CI passes
- [ ] Commits follow the conventional commit format (`feat(scope): ...`, `fix: ...`)
- [ ] New processes follow the module conventions in [CONTRIBUTING.md](../CONTRIBUTING.md) (stub block, `versions.yml`, container pinned, `meta.yml`)
- [ ] New parameters are declared in `nextflow_schema.json` and documented in `docs/parameters.md`
- [ ] `CHANGELOG.md` updated for user-visible changes
