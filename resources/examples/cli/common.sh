#!/usr/bin/env bash

rhbm_gem_checked_locations=()

rhbm_gem_usage_example() {
    cat <<'EOF'
Examples:
  bash resources/examples/cli/00_quickstart.sh /tmp/rhbm_cli_demo
  bash resources/examples/cli/00_quickstart.sh --workdir /tmp/rhbm_cli_demo --executable ./build-local/bin/RHBM-GEM
  RHBM_GEM_EXECUTABLE="$HOME/.local/bin/RHBM-GEM" bash resources/examples/cli/00_quickstart.sh --workdir /tmp/rhbm_cli_demo
EOF
}

rhbm_gem_canonicalize_existing_path() {
    local target="$1"

    if [ -d "$target" ]; then
        (
            cd "$target" >/dev/null 2>&1 &&
                pwd -P
        )
        return
    fi

    (
        cd "$(dirname "$target")" >/dev/null 2>&1 &&
            printf '%s/%s\n' "$(pwd -P)" "$(basename "$target")"
    )
}

rhbm_gem_append_checked_location() {
    local candidate="$1"

    rhbm_gem_checked_locations+=("$candidate")
}

rhbm_gem_try_candidate() {
    local candidate="$1"
    local normalized

    if [ -z "$candidate" ]; then
        return 1
    fi

    if normalized=$(rhbm_gem_canonicalize_existing_path "$candidate" 2>/dev/null); then
        rhbm_gem_append_checked_location "$normalized"
    else
        rhbm_gem_append_checked_location "$candidate"
    fi

    if [ -x "$candidate" ]; then
        rhbm_gem_canonicalize_existing_path "$candidate"
        return 0
    fi

    return 1
}

rhbm_gem_print_resolution_error() {
    local script_name="$1"
    local note="${2:-}"

    {
        echo "Error: Could not locate an executable RHBM-GEM binary."
        if [ -n "$note" ]; then
            echo "$note"
        fi
        echo "Checked:"
        if [ "${#rhbm_gem_checked_locations[@]}" -eq 0 ]; then
            echo "  - (no candidates were checked)"
        else
            local candidate
            for candidate in "${rhbm_gem_checked_locations[@]}"; do
                echo "  - $candidate"
            done
        fi
        echo "Recovery:"
        echo "  - Pass --executable /path/to/RHBM-GEM"
        echo "  - Or export RHBM_GEM_EXECUTABLE=/path/to/RHBM-GEM"
        echo "  - Or install RHBM-GEM so 'command -v RHBM-GEM' succeeds"
        echo "Usage:"
        echo "  bash $script_name --workdir /tmp/rhbm_cli_demo --executable /path/to/RHBM-GEM"
    } >&2
}

rhbm_gem_resolve_executable() {
    local explicit_executable="${1:-}"
    local script_path="${2:-${BASH_SOURCE[0]}}"
    local script_dir
    local repo_root
    local installed_prefix
    local path_candidate
    local candidate
    local resolved
    local nullglob_state

    rhbm_gem_checked_locations=()
    script_dir="$(
        cd "$(dirname "$script_path")" >/dev/null 2>&1 &&
            pwd -P
    )"

    if [ -n "$explicit_executable" ]; then
        if resolved=$(rhbm_gem_try_candidate "$explicit_executable"); then
            printf '%s\n' "$resolved"
            return 0
        fi
        rhbm_gem_print_resolution_error "$script_path" "The path from --executable is not executable."
        return 1
    fi

    if [ -n "${RHBM_GEM_EXECUTABLE:-}" ]; then
        if resolved=$(rhbm_gem_try_candidate "${RHBM_GEM_EXECUTABLE}"); then
            printf '%s\n' "$resolved"
            return 0
        fi
        rhbm_gem_print_resolution_error "$script_path" "The path from RHBM_GEM_EXECUTABLE is not executable."
        return 1
    fi

    if path_candidate=$(command -v RHBM-GEM 2>/dev/null); then
        if resolved=$(rhbm_gem_try_candidate "$path_candidate"); then
            printf '%s\n' "$resolved"
            return 0
        fi
    else
        rhbm_gem_append_checked_location "PATH:RHBM-GEM"
    fi

    installed_prefix="$(
        cd "${script_dir}/../../../../.." >/dev/null 2>&1 &&
            pwd -P
    )"
    candidate="${installed_prefix}/bin/RHBM-GEM"
    if resolved=$(rhbm_gem_try_candidate "$candidate"); then
        printf '%s\n' "$resolved"
        return 0
    fi

    repo_root="$(
        cd "${script_dir}/../../.." >/dev/null 2>&1 &&
            pwd -P
    )"
    candidate="${repo_root}/build-local/bin/RHBM-GEM"
    if resolved=$(rhbm_gem_try_candidate "$candidate"); then
        printf '%s\n' "$resolved"
        return 0
    fi

    candidate="${repo_root}/build/bin/RHBM-GEM"
    if resolved=$(rhbm_gem_try_candidate "$candidate"); then
        printf '%s\n' "$resolved"
        return 0
    fi

    nullglob_state="$(shopt -p nullglob || true)"
    shopt -s nullglob
    for candidate in "${repo_root}"/build-*/bin/RHBM-GEM; do
        case "$candidate" in
            "${repo_root}/build-local/bin/RHBM-GEM"|"${repo_root}/build/bin/RHBM-GEM")
                continue
                ;;
        esac

        if resolved=$(rhbm_gem_try_candidate "$candidate"); then
            eval "$nullglob_state"
            printf '%s\n' "$resolved"
            return 0
        fi
    done
    eval "$nullglob_state"

    rhbm_gem_print_resolution_error "$script_path"
    return 1
}

rhbm_gem_download_file() {
    local url="$1"

    if command -v curl >/dev/null 2>&1; then
        curl -fL -O "$url"
        return 0
    fi

    if command -v wget >/dev/null 2>&1; then
        wget "$url"
        return 0
    fi

    {
        echo "Error: Neither curl nor wget is available."
        echo "Install one of them, or rerun with --skip-download after preparing the data files manually."
    } >&2
    return 1
}
