#!/usr/bin/env bash
set -euo pipefail

project_root="$1"
tmp_root="$(mktemp -d "${TMPDIR:-/tmp}/rhbm-cli-quickstart.XXXXXX")"
trap 'rm -rf "$tmp_root"' EXIT

assert_file_contains() {
    local path="$1"
    local expected="$2"

    if ! grep -Fq "$expected" "$path"; then
        echo "Expected '$expected' in $path" >&2
        echo "Actual content:" >&2
        cat "$path" >&2
        exit 1
    fi
}

assert_text_contains() {
    local text="$1"
    local expected="$2"

    if [[ "$text" != *"$expected"* ]]; then
        echo "Expected output to contain: $expected" >&2
        echo "Actual output:" >&2
        printf '%s\n' "$text" >&2
        exit 1
    fi
}

copy_cli_sources() {
    local destination_dir="$1"

    mkdir -p "$destination_dir"
    cp "${project_root}/resources/examples/cli/00_quickstart.sh" "$destination_dir/00_quickstart.sh"
    cp "${project_root}/resources/examples/cli/common.sh" "$destination_dir/common.sh"
}

create_stub_executable() {
    local path="$1"
    local marker="$2"

    mkdir -p "$(dirname "$path")"
    cat >"$path" <<EOF
#!/usr/bin/env bash
set -euo pipefail
printf '%s|%s\n' "$marker" "\$*" >> "\${RHBM_GEM_STUB_LOG:?}"
EOF
    chmod +x "$path"
}

prepare_workdir_inputs() {
    local workdir="$1"

    mkdir -p "$workdir/data"
    : >"$workdir/data/6Z6U.cif"
    : >"$workdir/data/emd_11103_additional.map"
}

run_quickstart() {
    local script_path="$1"
    local workdir="$2"
    local log_path="$3"
    local -a env_args=()
    local -a script_args=()
    local -a command=()
    shift 3

    while [ $# -gt 0 ]; do
        case "$1" in
            *=*)
                env_args+=("$1")
                ;;
            *)
                script_args+=("$1")
                ;;
        esac
        shift
    done

    prepare_workdir_inputs "$workdir"
    command=(
        env
        "PATH=/usr/bin:/bin:/usr/sbin:/sbin"
        "RHBM_GEM_STUB_LOG=$log_path"
    )
    if [ "${#env_args[@]}" -gt 0 ]; then
        command+=("${env_args[@]}")
    fi
    command+=(
        bash
        "$script_path"
        --workdir "$workdir"
        --skip-download
    )
    if [ "${#script_args[@]}" -gt 0 ]; then
        command+=("${script_args[@]}")
    fi
    "${command[@]}"
}

scenario_explicit_override_wins() {
    local scenario_root="$tmp_root/explicit"
    local script_path="$scenario_root/repo/resources/examples/cli/00_quickstart.sh"
    local explicit_exec="$scenario_root/explicit/RHBM-GEM"
    local auto_exec="$scenario_root/repo/build-local/bin/RHBM-GEM"
    local log_path="$scenario_root/run.log"

    copy_cli_sources "$scenario_root/repo/resources/examples/cli"
    create_stub_executable "$explicit_exec" "explicit"
    create_stub_executable "$auto_exec" "auto"

    run_quickstart "$script_path" "$scenario_root/workdir" "$log_path" \
        --executable "$explicit_exec"

    assert_file_contains "$log_path" "explicit|potential_analysis"
    assert_file_contains "$log_path" "explicit|result_dump"
    if grep -Fq "auto|" "$log_path"; then
        echo "Expected --executable to take precedence over auto discovery." >&2
        exit 1
    fi
}

scenario_env_override_works() {
    local scenario_root="$tmp_root/env"
    local script_path="$scenario_root/repo/resources/examples/cli/00_quickstart.sh"
    local env_exec="$scenario_root/env/RHBM-GEM"
    local log_path="$scenario_root/run.log"

    copy_cli_sources "$scenario_root/repo/resources/examples/cli"
    create_stub_executable "$env_exec" "env"

    run_quickstart "$script_path" "$scenario_root/workdir" "$log_path" \
        RHBM_GEM_EXECUTABLE="$env_exec"

    assert_file_contains "$log_path" "env|potential_analysis"
    assert_file_contains "$log_path" "env|result_dump"
}

scenario_installed_tree_lookup_works() {
    local scenario_root="$tmp_root/installed"
    local script_path="$scenario_root/prefix/share/RHBM_GEM/resources/examples/cli/00_quickstart.sh"
    local installed_exec="$scenario_root/prefix/bin/RHBM-GEM"
    local log_path="$scenario_root/run.log"

    copy_cli_sources "$scenario_root/prefix/share/RHBM_GEM/resources/examples/cli"
    create_stub_executable "$installed_exec" "installed"

    run_quickstart "$script_path" "$scenario_root/workdir" "$log_path"

    assert_file_contains "$log_path" "installed|potential_analysis"
    assert_file_contains "$log_path" "installed|result_dump"
}

scenario_source_tree_build_local_lookup_works() {
    local scenario_root="$tmp_root/build-local"
    local script_path="$scenario_root/repo/resources/examples/cli/00_quickstart.sh"
    local build_exec="$scenario_root/repo/build-local/bin/RHBM-GEM"
    local log_path="$scenario_root/run.log"

    copy_cli_sources "$scenario_root/repo/resources/examples/cli"
    create_stub_executable "$build_exec" "build-local"

    run_quickstart "$script_path" "$scenario_root/workdir" "$log_path"

    assert_file_contains "$log_path" "build-local|potential_analysis"
    assert_file_contains "$log_path" "build-local|result_dump"
}

scenario_source_tree_build_lookup_works() {
    local scenario_root="$tmp_root/build"
    local script_path="$scenario_root/repo/resources/examples/cli/00_quickstart.sh"
    local build_exec="$scenario_root/repo/build/bin/RHBM-GEM"
    local log_path="$scenario_root/run.log"

    copy_cli_sources "$scenario_root/repo/resources/examples/cli"
    create_stub_executable "$build_exec" "build"

    run_quickstart "$script_path" "$scenario_root/workdir" "$log_path"

    assert_file_contains "$log_path" "build|potential_analysis"
    assert_file_contains "$log_path" "build|result_dump"
}

scenario_failure_output_is_actionable() {
    local scenario_root="$tmp_root/failure"
    local script_path="$scenario_root/repo/resources/examples/cli/00_quickstart.sh"
    local output

    copy_cli_sources "$scenario_root/repo/resources/examples/cli"
    prepare_workdir_inputs "$scenario_root/workdir"

    set +e
    output="$(
        PATH="/usr/bin:/bin:/usr/sbin:/sbin" \
            bash "$script_path" --workdir "$scenario_root/workdir" --skip-download 2>&1
    )"
    local status=$?
    set -e

    if [ "$status" -eq 0 ]; then
        echo "Expected quickstart to fail when no executable is available." >&2
        exit 1
    fi

    assert_text_contains "$output" "Could not locate an executable RHBM-GEM binary."
    assert_text_contains "$output" "Pass --executable /path/to/RHBM-GEM"
    assert_text_contains "$output" "PATH:RHBM-GEM"
}

scenario_explicit_override_wins
scenario_env_override_works
scenario_installed_tree_lookup_works
scenario_source_tree_build_local_lookup_works
scenario_source_tree_build_lookup_works
scenario_failure_output_is_actionable
