# Endpoint-reference acceptance evidence

The historical campaign runner is no longer part of the active test
architecture. Current search measurements are documented in the
[benchmark guide](../../joint-benchmark.md); reference certification remains
owned by the permanent offline tests.

[Acceptance report](../../joint-component-reference-acceptance.md) ·
[Summary](summary.json) · [Verification](verification.json) ·
[Archive member layout](archives.json) · [Campaign log](campaign.txt)

The six raw campaign archives are retrievable from Git history. Their former
paths, source and retrieval commits, compressed sizes, and SHA-256 values are
recorded in the [archive inventory](archives.json). Retrieve a blob with:

    git show <last_present_commit>:<former_path>

Recreating the report requires its historical source and environment. Current
benchmark runs use the entry point documented in the benchmark guide.
