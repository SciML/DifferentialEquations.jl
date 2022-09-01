name: Invalidations

on:
  pull_request:

concurrency:
  # Skip intermediate builds: always.
  # Cancel intermediate builds: always.
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true

jobs:
  evaluate:
    # Only run on PRs to the default branch.
    # In the PR trigger above branches can be specified only explicitly whereas this check should work for master, main, or any other default branch
    if: github.base_ref == github.event.repository.default_branch
    runs-on: ubuntu-latest
    steps:
    - uses: julia-actions/setup-julia@v1
      with:
        version: '1'
    - uses: actions/checkout@v3
    - uses: julia-actions/julia-buildpkg@v1
    - uses: julia-actions/julia-invalidations@v1
      id: invs_pr

    - uses: actions/checkout@v3
      with:
        ref: ${{ github.event.repository.default_branch }}
    - uses: julia-actions/julia-buildpkg@v1
    - uses: julia-actions/julia-invalidations@v1
      id: invs_default
    
    - name: Report invalidation counts
      run: |
        echo "Invalidations on default branch: ${{ steps.invs_default.outputs.total }} (${{ steps.invs_default.outputs.deps }} via deps)" >> $GITHUB_STEP_SUMMARY
        echo "This branch: ${{ steps.invs_pr.outputs.total }} (${{ steps.invs_pr.outputs.deps }} via deps)" >> $GITHUB_STEP_SUMMARY
    - name: Check if the PR does increase number of invalidations
      if: steps.invs_pr.outputs.total > steps.invs_default.outputs.total
      run: exit 1
