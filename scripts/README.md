# Scripts Directory Policy

- Keep reusable project scripts under `scripts/`.
- Put machine-specific one-off scripts under `scripts/local/` (ignored by git).
- Local scripts must use environment variables for machine-specific paths.
- Do not hardcode absolute user paths in scripts.
- Use `scripts/run_local_template.sh` as the baseline template.
- Do not place scripts at repository root.
