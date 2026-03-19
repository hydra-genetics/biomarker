# Release-Please Migration Guide — Snakemake Pipeline Repos

Some important steps:
- Remebmer to run 'git fetch --all' to make sure that all information is up to date
- Make sure that we are on a branch called update-please-release
- please-release-branch should originate from master

This guide covers migrating hydra-genetics Snakemake pipeline repositories to
`googleapis/release-please-action@v4`. These steps differ from the main
`hydra-genetics` CLI repo in two ways: there is no PyPI publish step, and
`release-type: simple` is used since there are no Python package files.

## Key differences from the hydra-genetics CLI repo

| | CLI repo | Pipeline repos |
|---|---|---|
| `release-type` | `python` | `simple` |
| PyPI publish job | yes | no |
| `pyproject.toml` changes | yes | no |

---

## Step 1 — Get current version and bootstrap SHA

```bash
# Latest tag → goes into the manifest
git tag --sort=-v:refname | head -1

# SHA of the commit that tag points to → goes into bootstrap-sha
git rev-list -n 1 $(git tag --sort=-v:refname | head -1)
```

The bootstrap SHA tells release-please where to start scanning commits when
creating the first release PR, preventing it from scanning the entire git
history.

---

## Step 2 — Create `.release-please-manifest.json`

```json
{
  ".": "X.Y.Z"
}
```

Replace `X.Y.Z` with the version from Step 1, without the `v` prefix.

---

## Step 3 — Create `release-please-config.json`

```json
{
  "$schema": "https://raw.githubusercontent.com/googleapis/release-please/main/schemas/config.json",
  "packages": {
    ".": {
      "release-type": "simple",
      "package-name": "<repo-name>",
      "changelog-path": "CHANGELOG.md",
      "include-v-in-tag": true,
      "include-component-in-tag": false,
      "draft": false,
      "prerelease": false
    }
  },
  "bootstrap-sha": "<sha-from-step-1>"
}
```

Replace `<repo-name>` with the repository name (e.g. `prealignment`) and
`<sha-from-step-1>` with the SHA from Step 1.

---

## Step 4 — Replace `.github/workflows/release-please.yaml`

```yaml
name: release-please

on:
  push:
    branches:
      - master

permissions:
  contents: write
  pull-requests: write

jobs:
  release-please:
    runs-on: ubuntu-latest
    steps:
      - name: Generate Token
        uses: actions/create-github-app-token@v2
        id: app-token
        with:
          app-id: ${{ secrets.APP_ID }}
          private-key: ${{ secrets.APP_PRIVATE_KEY }}

      - uses: googleapis/release-please-action@v4
        id: release
        with:
          token: ${{ steps.app-token.outputs.token }}
          target-branch: master
```

---

## Step 5 — Add secrets to the repository

In each repo: **Settings → Secrets and variables → Actions**

| Secret | Value |
|---|---|
| `APP_ID` | GitHub App ID (same as CLI repo) |
| `APP_PRIVATE_KEY` | GitHub App private key (same as CLI repo) |

---

## Step 6 — Confirm the GitHub App is installed on the repo

**GitHub → (App settings) → Install App → select the repository**

Skip this step if the App is already installed organisation-wide.

---

## Checklist per repo

```
[ ] Get current version and bootstrap SHA (Step 1)
[ ] Create .release-please-manifest.json (Step 2)
[ ] Create release-please-config.json (Step 3)
[ ] Replace .github/workflows/release-please.yaml (Step 4)
[ ] Add APP_ID and APP_PRIVATE_KEY secrets (Step 5)
[ ] Confirm GitHub App is installed on the repo (Step 6)
[ ] Push to master and verify the first release PR is created correctly
```
