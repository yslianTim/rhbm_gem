# Endpoint-reference acceptance evidence

Status: Historical. The campaign runner was retired. Current search measurements
are documented in [`../../joint-benchmark.md`](../../joint-benchmark.md), and
reference certification remains owned by the permanent offline tests.

[Acceptance report](../../joint-component-reference-acceptance.md) ·
[Summary](summary.json) · [Verification](verification.json) ·
[Historical file layout](archives.json) · [Campaign log](campaign.txt) ·
[Artifact manifest](../experiment-retirement-baseline/artifact-manifest.json)

The six raw campaign archives were removed from the current worktree in PR 6.
Their compressed sizes, SHA-256 values, source commits and former paths are
recorded as `removed_from_worktree` entries in the artifact manifest. The
retained summary and verification records preserve the campaign conclusions and
provenance; the archive blobs remain retrievable from Git history.

To retrieve an exact historical blob, use its `last_present_commit` and `path`
from the manifest with `git show <commit>:<path>`. The retired campaign source
was last updated at `7f84f931215f7b5a52948f7ecdbfc15fad79d7c5`; recreating the old
report requires that historical source and environment. Current benchmark runs
use the entry point in [`joint-benchmark.md`](../../joint-benchmark.md).
